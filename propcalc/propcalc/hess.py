"""원소분석 + HHV 로부터 생성엔탈피(ΔfH°)를 역산한다 (헤스의 법칙).

  연료 CcHhOoNnSs + O2 -> c CO2(g) + h/2 H2O(l) + n/2 N2 + s SO2(g)
  ΔfH°(연료) = c·ΔfH°(CO2) + h/2·ΔfH°(H2O,l) + s·ΔfH°(SO2) − ΔcH°
  ΔcH° = −HHV · M    [MJ/kg · g/mol = kJ/mol]

단순화: 폭탄열량계의 정적 연소열(ΔcU)과 ΔcH° 의 차이(Δn_gas·RT, 리그닌 정도에서
약 −3.5 kJ/mol)는 무시한다. 이 크기는 HHV 측정 불확실도보다 훨씬 작다.
"""

from __future__ import annotations

import dataclasses
import warnings
from collections.abc import Mapping
from dataclasses import dataclass, field

from .formula import ATOMIC_WEIGHTS, as_formula_dict, molar_mass
from .mixture import Ingredient

# 생성물 기준값 (kJ/mol, 298.15 K). 출처: CODATA Key Values for Thermodynamics
# (Cox, Wagman, Medvedev 1989). CO2, H2O 는 지시서 기재값과 일치. SO2 는 CODATA 값
# -296.81 을 사용한다 (이 세션에서 원문 미확인, 영향은 황 함유 시료에 한정).
DHF_CO2_G = -393.51
DHF_H2O_L = -285.83
DHF_SO2_G = -296.81

# Channiwala & Parikh 2002 (Fuel 81, 1051) 상관식이 보고한 평균 절대 오차로 알려진 값.
# 원문 미확인이므로 기본 불확실도 추정용으로만 쓴다.
CP_DEFAULT_REL_UNCERTAINTY = 0.0145

HHV_CORRELATION_WARNING = (
    "생성엔탈피는 큰 두 수의 차이므로 HHV 1~2% 오차가 ΔfH°에서 수 % 이상으로 증폭됨. "
    "실측 HHV 사용 권장"
)


def channiwala_parikh_hhv(C: float, H: float, O: float, N: float = 0.0, S: float = 0.0,
                          ash: float = 0.0) -> float:
    """HHV [MJ/kg] = 0.3491C + 1.1783H + 0.1005S − 0.1034O − 0.0151N − 0.0211A (질량%, 건조 기준)."""
    return 0.3491 * C + 1.1783 * H + 0.1005 * S - 0.1034 * O - 0.0151 * N - 0.0211 * ash


@dataclass(frozen=True)
class HessResult:
    formula: dict[str, float]
    M_g_per_mol: float
    hhv_MJ_per_kg: float
    hhv_source: str  # "measured" | "correlation"
    dHf_kJ_per_mol: float
    u_dHf_kJ_per_mol: float | None
    dcH_kJ_per_mol: float
    warnings: tuple[str, ...] = field(default_factory=tuple)


def dHf_from_hhv(
    formula: str | Mapping[str, float],
    hhv_MJ_per_kg: float,
    u_hhv_MJ_per_kg: float | None = None,
    hhv_source: str = "measured",
) -> HessResult:
    f = as_formula_dict(formula)
    unsupported = set(f) - {"C", "H", "O", "N", "S"}
    if unsupported:
        raise NotImplementedError(f"회분/금속 원소는 지원하지 않는다: {sorted(unsupported)}")
    if hhv_MJ_per_kg <= 0:
        raise ValueError("HHV 는 양수여야 한다")
    M = molar_mass(f)
    dcH = -hhv_MJ_per_kg * M
    dHf = (f.get("C", 0.0) * DHF_CO2_G + f.get("H", 0.0) / 2 * DHF_H2O_L
           + f.get("S", 0.0) * DHF_SO2_G - dcH)
    u = None if u_hhv_MJ_per_kg is None else u_hhv_MJ_per_kg * M
    warns = (HHV_CORRELATION_WARNING,) if hhv_source == "correlation" else ()
    for w in warns:
        warnings.warn(w, UserWarning, stacklevel=3)
    return HessResult(f, M, hhv_MJ_per_kg, hhv_source, dHf, u, dcH, warns)


def formula_from_elemental_analysis(C: float, H: float, O: float, N: float = 0.0, S: float = 0.0,
                                    basis_C: float = 10.0) -> dict[str, float]:
    """질량% 원소분석을 탄소 basis_C 개 기준 화학식으로 환산한다."""
    mol = {"C": C / ATOMIC_WEIGHTS["C"], "H": H / ATOMIC_WEIGHTS["H"], "O": O / ATOMIC_WEIGHTS["O"],
           "N": N / ATOMIC_WEIGHTS["N"], "S": S / ATOMIC_WEIGHTS["S"]}
    scale = basis_C / mol["C"]
    return {el: n * scale for el, n in mol.items() if n > 0}


def dHf_from_elemental_analysis(
    C: float, H: float, O: float, N: float = 0.0, S: float = 0.0, ash: float = 0.0,
    hhv_MJ_per_kg: float | None = None,
    u_hhv_MJ_per_kg: float | None = None,
    formula: str | Mapping[str, float] | None = None,
    basis_C: float = 10.0,
) -> HessResult:
    """원소분석(질량%, 건조·무회분 기준) 에서 ΔfH° 를 구한다.

    hhv_MJ_per_kg 를 주면 실측값으로, 주지 않으면 Channiwala–Parikh 상관식을 쓰고 경고한다.
    formula 를 주지 않으면 원소분석에서 C=basis_C 기준 화학식을 만든다.
    """
    total = C + H + O + N + S + ash
    if abs(total - 100.0) > 1.0:
        raise ValueError(f"원소분석 합이 100% 에서 1% 이상 벗어난다: {total:.2f} (회분은 ash 로 지정)")
    f = as_formula_dict(formula) if formula is not None else formula_from_elemental_analysis(
        C, H, O, N, S, basis_C)
    if hhv_MJ_per_kg is None:
        hhv = channiwala_parikh_hhv(C, H, O, N, S, ash)
        u = u_hhv_MJ_per_kg if u_hhv_MJ_per_kg is not None else CP_DEFAULT_REL_UNCERTAINTY * hhv
        return dHf_from_hhv(f, hhv, u, hhv_source="correlation")
    return dHf_from_hhv(f, hhv_MJ_per_kg, u_hhv_MJ_per_kg, hhv_source="measured")


def dHf_from_dry_basis_report(
    C: float, H: float, N: float, S: float, ash: float, hhv_dry_MJ_per_kg: float,
    O: float | None = None,
    u_hhv_dry_MJ_per_kg: float | None = None,
    formula: str | Mapping[str, float] | None = None,
    basis_C: float = 10.0,
    closure_tol: float = 1.5,
) -> HessResult:
    """건조 기준(회분 포함) 분석 보고서를 건조·무회분(daf) 기준으로 환산해 ΔfH° 를 역산한다.

    입력은 모두 같은 건조 시료의 값이어야 한다 (질량%, HHV 는 MJ/kg-건조). O 를 주지 않으면
    O = 100 − (C+H+N+S+회분) 로 차이법 계산하며, 이 경우 다른 원소 분석 오차가 모두 O 에 쌓인다.
    O 를 주면 합이 100±closure_tol % 안이어야 한다.
    회분은 반응물 조성에서 제외된다(이 모델은 C,H,O,N,S 만 지원). 따라서 결과 연료의 질량%는
    무회분 유기물 기준이며, 실제 시료의 회분 질량만큼 배합 질량이 어긋난다.
    """
    if not 0 <= ash < 100:
        raise ValueError(f"회분은 0 이상 100 미만이어야 한다: {ash}")
    warns: list[str] = []
    organic = C + H + N + S + ash
    if O is None:
        O = 100.0 - organic
        if O < 0:
            raise ValueError(f"차이법 산소가 음수다 ({O:.2f}%): C+H+N+S+회분 이 100% 를 넘는다")
        warns.append(f"산소는 차이법({O:.2f}%)이다. 다른 원소 분석 오차가 모두 O 에 쌓이므로 직접 측정을 권장")
    elif abs(organic + O - 100.0) > closure_tol:
        raise ValueError(f"C+H+N+S+O+회분 의 합이 {organic + O:.2f}% 이다 (허용 100±{closure_tol}%). "
                         "분석 항목 누락이나 기준(건량/습량) 혼용을 확인하라")
    k = 1.0 - ash / 100.0
    daf = dict(C=C / k, H=H / k, N=N / k, S=S / k, O=O / k)
    scale = 100.0 / sum(daf.values())  # 합 오차를 무회분 기준에서 정규화
    daf = {el: v * scale for el, v in daf.items()}
    f = as_formula_dict(formula) if formula is not None else formula_from_elemental_analysis(
        daf["C"], daf["H"], daf["O"], daf["N"], daf["S"], basis_C)
    hhv_daf = hhv_dry_MJ_per_kg / k
    u_daf = None if u_hhv_dry_MJ_per_kg is None else u_hhv_dry_MJ_per_kg / k
    r = dHf_from_hhv(f, hhv_daf, u_daf, hhv_source="measured")
    warns.append(f"회분 {ash:g}% 는 반응물에서 제외됨 (무회분 기준 ΔfH°). 배합 질량%는 무회분 유기물 기준")
    return dataclasses.replace(r, warnings=r.warnings + tuple(warns))


def ingredient_from_hhv(key: str, name: str, formula: str | Mapping[str, float], hhv_MJ_per_kg: float,
                        u_hhv_MJ_per_kg: float | None, source: str,
                        hhv_source: str = "measured") -> Ingredient:
    """HHV 로 역산한 ΔfH° 를 가진 Ingredient 를 만든다 (리그닌 등 정의되지 않은 분자용)."""
    r = dHf_from_hhv(formula, hhv_MJ_per_kg, u_hhv_MJ_per_kg, hhv_source)
    return Ingredient(
        key=key, name=name, formula=r.formula, dHf_kJ_per_mol=r.dHf_kJ_per_mol,
        u_dHf_kJ_per_mol=r.u_dHf_kJ_per_mol, M_g_per_mol=r.M_g_per_mol, phase="cr",
        source=f"헤스 역산: HHV={hhv_MJ_per_kg:g} MJ/kg ({hhv_source}); {source}",
        verification="secondary" if hhv_source == "measured" else "doc_only",
        note="원소분석과 HHV 는 같은 시료에서 얻은 값이어야 한다." if hhv_source == "measured" else HHV_CORRELATION_WARNING,
    )
