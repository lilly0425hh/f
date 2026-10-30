"""NASA cea (https://github.com/nasa/cea, Apache-2.0) 백엔드.

API 근거 (이 세션에서 설치된 cea 3.3.4 의 stub, 샘플(rp1311/example5,8,12,13), 실측으로 확인):
- 사용자 반응물: cea.Reactant(name, formula, molecular_weight, enthalpy, enthalpy_units, temperature)
- 로켓 문제: RocketSolver.solve(soln, weights, pc[bar], pi_p=[Pc/Pe], supar=[Ae/At], iac=True, hc=h/R, n_frz=)
- 스테이션 순서(IAC): 결과 배열 인덱스 0=연소실, 1=목, 2=출구. n_frz 는 1부터 세는 스테이션 번호이며
  1=연소실부터 동결, 2=목부터 동결.
  n_frz=2 는 nasa.github.io/cea Example 12("composition frozen at the throat", n_frz = 2)와 그 출력
  (연소실 열만 조성이 다르고 목 이후 열의 질량분율이 모두 같음)으로 2026-10-30 직접 확인했다.
  n_frz=1 이 연소실 동결이라는 것은 웹 문서에 명시가 없어 샘플 주석과 rocketcea 문서 문자열에서 추론한 값이다.
- 응축상 종은 species_names 에서 기체 종 뒤(인덱스 >= num_gas)에 위치한다.
- calc_property(ENTHALPY) 는 J/kg, nj 는 mol/g 단위.
"""

from __future__ import annotations

import warnings

import numpy as np

from ..mixture import Mixture
from .base import G0, Backend, ConvergenceError, Result, check_physical, validate_run_args

TRACE_CUTOFF = 1e-12  # mol/kg 미만 종은 결과에서 생략
MAX_NAME_LEN = 15     # cea.Reactant 이름 최대 길이


class CEABackend(Backend):
    name = "cea"

    def __init__(self) -> None:
        try:
            import cea  # noqa: F401
        except ImportError as e:  # pragma: no cover
            raise ImportError("NASA cea 패키지가 필요하다: pip install cea (Python>=3.11)") from e
        import cea as _cea

        self._cea = _cea

    def run(self, mixture: Mixture, Pc_MPa: float, Pe_MPa: float | None = None,
            area_ratio: float | None = None, mode: str = "equilibrium",
            frozen_at: str = "chamber") -> Result:
        validate_run_args(Pc_MPa, Pe_MPa, area_ratio, mode)
        if frozen_at not in ("chamber", "throat"):
            raise ValueError("frozen_at 은 'chamber' 또는 'throat'")
        cea = self._cea

        active = [(ing, w) for ing, w in mixture.parts if w > 0]
        # 고체 KNO3 처럼 CEA DB 에 없거나 상(phase)이 다른 반응물은 항상 커스텀 반응물로 전달한다.
                # 반응물 이름은 15자 이하여야 한다 (16자부터 CEA_INVALID_SIZE, 이 세션에서 실측).
        # 인덱스 접두어로 유일성을 보장하고 DB 종 이름과의 충돌도 피한다.
        reactants = [
            cea.Reactant(
                name=f"R{i}_{ing.key}"[:MAX_NAME_LEN], formula={k: float(v) for k, v in ing.formula.items()},
                molecular_weight=ing.M_g_per_mol, enthalpy=ing.dHf_kJ_per_mol,
                enthalpy_units="kJ/mol", temperature=mixture.T_K,
            )
            for i, (ing, _) in enumerate(active)
        ]
        weights = np.array([w for _, w in active], dtype=float) / 100.0
        T0 = np.full(len(active), mixture.T_K)

        reac = cea.Mixture(reactants)
        prod = cea.Mixture(reactants, products_from_reactants=True)
        solver = cea.RocketSolver(prod, reactants=reac)
        sol = cea.RocketSolution(solver)
        hc = reac.calc_property(cea.ENTHALPY, weights, T0) / cea.R

        kwargs: dict = dict(iac=True, hc=hc)
        pc_bar = Pc_MPa * 10.0
        if Pe_MPa is not None:
            kwargs["pi_p"] = [Pc_MPa / Pe_MPa]
        else:
            kwargs["supar"] = [float(area_ratio)]
        if mode == "frozen":
            kwargs["n_frz"] = 1 if frozen_at == "chamber" else 2

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            solver.solve(sol, weights, pc_bar, **kwargs)
        msgs = [str(w.message) for w in caught]

        if not sol.converged or sol.num_pts < 3:
            raise ConvergenceError(
                f"cea 가 수렴하지 못했다 (mode={mode}, last_error={sol.last_error}, num_pts={sol.num_pts}). "
                "동결 유동은 응축상이 있는 계에서 CEA 가 '온도가 응축상 유효 범위를 50 K 이상 벗어남'으로 "
                "중단할 수 있다. " + "; ".join(msgs)
            )

        names = list(prod.species_names)
        condensed = frozenset(names[solver.num_gas:])
        nj_chamber = np.asarray(sol.nj)[0] * 1000.0  # mol/g -> mol/kg
        species = {n: float(x) for n, x in zip(names, nj_chamber) if x > TRACE_CUTOFF}
        mass_fr = sol.mass_fractions
        cond_mass = {n: float(mass_fr[n][0]) for n in condensed if mass_fr[n][0] > 1e-12}
        exit_cond = sum(float(mass_fr[n][-1]) for n in condensed)

        warns = list(msgs) + mixture.provenance_notes()
        res = Result(
            backend=self.name, mode=mode, Pc_MPa=Pc_MPa, Pe_MPa=float(sol.P[-1]) / 10.0,
            area_ratio=float(sol.ae_at[-1]), composition_wt_pct=dict(mixture.wt_pct()),
            Tc_K=float(sol.T[0]), M_eff_g_mol=float(sol.M[0]), MW_g_mol=float(sol.MW[0]),
            gamma_s=float(sol.gamma_s[0]), species_moles_per_kg=species,
            condensed_species=condensed, condensed_mass_fractions=cond_mass,
            c_star_m_s=float(sol.c_star[-1]), Cf=float(sol.coefficient_of_thrust[-1]),
            isp_s=float(sol.Isp[-1]) / G0, isp_vac_s=float(sol.Isp_vacuum[-1]) / G0,
            ve_m_s=float(sol.Isp[-1]), Te_K=float(sol.T[-1]), gamma_e=float(sol.gamma_s[-1]),
            M_e_g_mol=float(sol.M[-1]), Tt_K=float(sol.T[1]),
            exit_condensed_mass_fraction=exit_cond, warnings=warns,
        )
        check_physical(res)
        return res
