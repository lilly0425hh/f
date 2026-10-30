"""RocketCEA (레거시 CEA 래퍼, GPL-3.0) 백엔드 — 교차검증용.

추진제 전체를 하나의 단일 추진제(add_new_propellant)로 정의한다. 카드 형식:
  name <이름> <원소 개수 ...> wt%=<질량%> h,cal=<cal/mol> t(k)=<K>
내부 단위(psia, ft/s, degR)는 SI 로 변환한다. 이 백엔드는 평형 계산만 지원한다
(KNSB 에서는 동결 계산도 0.0 을 돌려주어 의미가 없다).
"""

from __future__ import annotations

import hashlib
import re

from ..mixture import Mixture
from .base import G0, Backend, Result, check_physical, validate_run_args

PSI_PER_PA = 1.0 / 6894.757293168
FT_TO_M = 0.3048
CAL_PER_J = 1.0 / 4.184  # 열화학 칼로리
_CONDENSED = re.compile(r"\((?:L|cr|gr|s|a|b|c|d|I{1,3}|IV)\)$")


def _card_name(ing_key: str, idx: int) -> str:
    base = re.sub(r"[^A-Za-z0-9]", "", ing_key)[:12] or "X"
    return f"{base}{idx}"


class RocketCEABackend(Backend):
    name = "rocketcea"

    def __init__(self) -> None:
        try:
            from rocketcea.cea_obj import CEA_Obj, add_new_propellant
        except ImportError as e:  # pragma: no cover
            raise ImportError("rocketcea 가 필요하다: pip install rocketcea (gfortran 필요)") from e
        self._CEA_Obj, self._add = CEA_Obj, add_new_propellant

    def _build_card(self, mixture: Mixture) -> str:
        lines = []
        for idx, (ing, w) in enumerate(p for p in mixture.parts if p[1] > 0):
            els = " ".join(f"{el} {n:.6g}" for el, n in ing.formula.items())
            if len(ing.formula) > 5:
                raise NotImplementedError("rocketcea 카드는 반응물당 원소 5개까지만 지원한다")
            h_cal = ing.dHf_kJ_per_mol * 1000.0 * CAL_PER_J
            lines.append(f"name {_card_name(ing.key, idx)} {els} wt%={w:.8g} h,cal={h_cal:.6f} t(k)={mixture.T_K:.2f}")
        return "\n" + "\n".join(lines) + "\n"

    def run(self, mixture: Mixture, Pc_MPa: float, Pe_MPa: float | None = None,
            area_ratio: float | None = None, mode: str = "equilibrium") -> Result:
        validate_run_args(Pc_MPa, Pe_MPa, area_ratio, mode)
        if mode != "equilibrium":
            raise NotImplementedError("rocketcea 백엔드는 평형 계산 교차검증 전용이다")

        card = self._build_card(mixture)
        prop_name = "PC" + hashlib.sha1(card.encode()).hexdigest()[:10]
        self._add(prop_name, card)
        cea = self._CEA_Obj(propName=prop_name)

        Pc_psia = Pc_MPa * 1e6 * PSI_PER_PA
        if Pe_MPa is not None:
            eps = float(cea.get_eps_at_PcOvPe(Pc=Pc_psia, MR=1.0, PcOvPe=Pc_MPa / Pe_MPa))
            Pe_out = Pe_MPa
        else:
            eps = float(area_ratio)
            Pe_out = Pc_MPa / float(cea.get_PcOvPe(Pc=Pc_psia, MR=1.0, eps=eps))
        Pe_psia = Pe_out * 1e6 * PSI_PER_PA

        isp_vac, cstar_fps, tc_R, _mw, gamma = cea.get_IvacCstrTc_ChmMwGam(Pc=Pc_psia, MR=1.0, eps=eps)
        cstar = float(cstar_fps) * FT_TO_M
        cf = float(cea.get_PambCf(Pc=Pc_psia, MR=1.0, eps=eps, Pamb=Pe_psia)[0])
        temps = cea.get_Temperatures(Pc=Pc_psia, MR=1.0, eps=eps)  # [연소실, 목, 출구] (degR)
        mw_c, gam_c = cea.get_Chamber_MolWt_gamma(Pc=Pc_psia, MR=1.0, eps=eps)
        mw_e, gam_e = cea.get_exit_MolWt_gamma(Pc=Pc_psia, MR=1.0, eps=eps)

        molw, massfr = cea.get_SpeciesMassFractions(Pc=Pc_psia, MR=1.0, eps=eps)
        species: dict[str, float] = {}
        cond_mass: dict[str, float] = {}
        condensed: set[str] = set()
        n_gas = 0.0           # mol/g
        exit_cond = 0.0
        for raw, fr in massfr.items():
            name = raw.lstrip("*")
            is_cond = bool(_CONDENSED.search(name))
            x_c = float(fr[1])
            if is_cond:
                condensed.add(name)
                exit_cond += float(fr[3])
                if x_c > 1e-12:
                    cond_mass[name] = x_c
            elif x_c > 0:
                n_gas += x_c / float(molw[raw])
            if x_c > 0:
                species[name] = x_c * 1000.0 / float(molw[raw])

        res = Result(
            backend=self.name, mode=mode, Pc_MPa=Pc_MPa, Pe_MPa=Pe_out, area_ratio=eps,
            composition_wt_pct=dict(mixture.wt_pct()),
            Tc_K=float(tc_R) / 1.8, M_eff_g_mol=1.0 / n_gas, MW_g_mol=float(mw_c),
            gamma_s=float(gam_c), species_moles_per_kg=species,
            condensed_species=frozenset(condensed), condensed_mass_fractions=cond_mass,
            c_star_m_s=cstar, Cf=cf, isp_s=cf * cstar / G0, isp_vac_s=float(isp_vac),
            ve_m_s=cf * cstar, Te_K=float(temps[2]) / 1.8, gamma_e=float(gam_e), M_e_g_mol=float(mw_e),
            Tt_K=float(temps[1]) / 1.8, exit_condensed_mass_fraction=exit_cond,
            warnings=mixture.provenance_notes(),
        )
        check_physical(res)
        return res
