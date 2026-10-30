"""성능 정리: c*, Cf, 추력, 비추력과 이상기체 해석해 비교.

모든 식은 SI (Pa, m^2, kg/s, m/s).

  c*  = Pc·At / ṁ
  F   = ṁ·ve + (Pe − Pa)·Ae
  Cf  = F/(Pc·At) = ve/c* + (Pe − Pa)·Ae/(Pc·At)
  Isp = F/(ṁ·g0) = Cf·c*/g0
  c*(이상기체) = sqrt(R̄·Tc/M) / [ sqrt(γ)·(2/(γ+1))^((γ+1)/(2(γ−1))) ]
"""

from __future__ import annotations

from dataclasses import dataclass

from .backends.base import G0, Result

R_UNIVERSAL = 8.314462618  # J/(mol K), CODATA 2018 (정확값)

# 아마추어 자료(Nakka)의 KNSB 실측 해면 Isp 137 s / 평형 계산 151.5 s.
# 신뢰도가 낮은 아마추어 자료이므로 참고용 경험 효율로만 쓰고, 사용자가 명시적으로 선택해야 한다.
ETA_ISP_KNSB_AMATEUR = 137.0 / 151.5

UPPER_BOUND_NOTE = (
    "평형 결과는 입자-기체 완전 평형(응축상과 기체의 속도·온도 슬립 없음), 단열, 1D 등엔트로피, "
    "무한 반응속도를 가정하므로 Isp 의 이론 상한이다. 실제 2상 유동 손실은 포함되지 않는다."
)


def _positive(name: str, value: float) -> None:
    if not value > 0:
        raise ValueError(f"{name} 는 양수여야 한다: {value}")


def c_star(Pc_Pa: float, At_m2: float, mdot_kg_s: float) -> float:
    _positive("Pc", Pc_Pa), _positive("At", At_m2), _positive("mdot", mdot_kg_s)
    return Pc_Pa * At_m2 / mdot_kg_s


def thrust(mdot_kg_s: float, ve_m_s: float, Pe_Pa: float, Pa_Pa: float, Ae_m2: float) -> float:
    return mdot_kg_s * ve_m_s + (Pe_Pa - Pa_Pa) * Ae_m2


def thrust_coefficient(ve_m_s: float, c_star_m_s: float, Pe_Pa: float, Pa_Pa: float,
                       Ae_m2: float, Pc_Pa: float, At_m2: float) -> float:
    return ve_m_s / c_star_m_s + (Pe_Pa - Pa_Pa) * Ae_m2 / (Pc_Pa * At_m2)


def isp_from_thrust(F_N: float, mdot_kg_s: float) -> float:
    _positive("mdot", mdot_kg_s)
    return F_N / (mdot_kg_s * G0)


def ideal_gas_c_star(Tc_K: float, M_g_per_mol: float, gamma: float) -> float:
    """이상기체 해석해 c* (m/s). M 은 g/mol."""
    _positive("Tc", Tc_K), _positive("M", M_g_per_mol)
    if not gamma > 1:
        raise ValueError(f"γ 는 1 보다 커야 한다: {gamma}")
    R = R_UNIVERSAL / (M_g_per_mol * 1e-3)
    gam_fn = gamma ** 0.5 * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    return (R * Tc_K) ** 0.5 / gam_fn


@dataclass(frozen=True)
class PerformanceReport:
    assumptions: str
    Pa_MPa: float
    c_star_equilibrium: float
    c_star_analytic_Meff: float
    c_star_analytic_MW: float
    c_star_analytic_Meff_rel_diff: float   # (해석해 − 평형)/평형, M_eff 사용
    c_star_analytic_MW_rel_diff: float     # (해석해 − 평형)/평형, MW 사용
    Cf_at_Pa: float
    isp_at_Pa_s: float
    ve_m_s: float
    note: str
    upper_bound_note: str
    isp_empirical_s: float | None = None
    empirical_note: str = ""

    def summary(self) -> str:
        lines = [
            f"[가정] {self.assumptions}",
            f"[대기압] Pa = {self.Pa_MPa:g} MPa -> Cf = {self.Cf_at_Pa:.4f}, Isp = {self.isp_at_Pa_s:.2f} s",
            f"[c*] 평형(CEA) {self.c_star_equilibrium:.1f} m/s | 해석해(M_eff) {self.c_star_analytic_Meff:.1f} "
            f"({self.c_star_analytic_Meff_rel_diff:+.2%}) | 해석해(MW) {self.c_star_analytic_MW:.1f} "
            f"({self.c_star_analytic_MW_rel_diff:+.2%})",
            f"[주석] {self.note}",
            f"[상한] {self.upper_bound_note}",
        ]
        if self.isp_empirical_s is not None:
            lines.append(f"[경험 효율 적용] Isp = {self.isp_empirical_s:.2f} s ({self.empirical_note})")
        return "\n".join(lines)


def report(result: Result, Pa_MPa: float | None = None, eta_isp: float | None = None) -> PerformanceReport:
    """Result 로부터 대기압 Pa 에서의 Cf, Isp 와 이상기체 해석해 비교를 만든다.

    Pa 가 None 이면 설계 출구압 Pe(이상팽창)를 쓴다. eta_isp 를 주면 경험 효율을 적용한
    Isp 를 별도로 표기한다 (기본은 적용하지 않음).
    """
    if eta_isp is not None and not 0 < eta_isp <= 1:
        raise ValueError("eta_isp 는 (0, 1] 이어야 한다")
    Pc, Pe = result.Pc_MPa * 1e6, result.Pe_MPa * 1e6
    Pa = Pe if Pa_MPa is None else Pa_MPa * 1e6
    eps = result.area_ratio
    # Ae/ṁ = ε·At/ṁ = ε·c*/Pc 이므로 At, ṁ 없이 무차원으로 환산된다
    cf = result.ve_m_s / result.c_star_m_s + (Pe - Pa) * eps / Pc
    isp = cf * result.c_star_m_s / G0

    cs_meff = ideal_gas_c_star(result.Tc_K, result.M_eff_g_mol, result.gamma_s)
    cs_mw = ideal_gas_c_star(result.Tc_K, result.MW_g_mol, result.gamma_s)
    d_meff = (cs_meff - result.c_star_m_s) / result.c_star_m_s
    d_mw = (cs_mw - result.c_star_m_s) / result.c_star_m_s

    xc = result.condensed_mass_fraction
    if xc > 0.01:
        note = (
            f"응축상 질량분율이 {xc:.1%} 이다. 해석해에 CEA 의 M_eff(시스템 질량/기체 몰수)와 γ_s 를 쓰면 "
            f"평형 c* 와 {d_meff:+.2%} 차이로 일치하는데, 이는 CEA 가 응축상을 무슬립·열평형인 유효 기체로 "
            f"다루는 것과 같은 가정이므로 독립 검증은 아니다. 응축상을 한 분자로 세는 MW 를 쓰면 "
            f"c* 가 {d_mw:+.2%} 달라져, 응축상 질량이 유효 분자량에 미치는 영향의 크기를 보여 준다."
        )
    else:
        note = f"응축상이 거의 없다. 해석해(M_eff) 차이 {d_meff:+.2%} 는 γ 가 노즐 내에서 변하기 때문이다."

    emp, emp_note = None, ""
    if eta_isp is not None:
        emp = isp * eta_isp
        emp_note = (f"η={eta_isp:.3f}. KNSB 참고값 137/151.5 는 아마추어 자료(Nakka)이며 신뢰도가 낮으므로 "
                    "논문에서는 실측 정적 시험값으로 대체해야 한다")
    return PerformanceReport(
        assumptions=result.assumptions(), Pa_MPa=Pa / 1e6, c_star_equilibrium=result.c_star_m_s,
        c_star_analytic_Meff=cs_meff, c_star_analytic_MW=cs_mw,
        c_star_analytic_Meff_rel_diff=d_meff, c_star_analytic_MW_rel_diff=d_mw,
        Cf_at_Pa=cf, isp_at_Pa_s=isp, ve_m_s=result.ve_m_s, note=note,
        upper_bound_note=UPPER_BOUND_NOTE, isp_empirical_s=emp, empirical_note=emp_note,
    )


@dataclass(frozen=True)
class MotorPoint:
    mdot_kg_s: float
    thrust_N: float
    isp_s: float
    Cf: float


def motor_point(result: Result, At_m2: float, Pa_MPa: float | None = None) -> MotorPoint:
    """목 면적 At 이 주어졌을 때의 정상상태 질량유량과 추력 (열역학적 스케일링만)."""
    _positive("At", At_m2)
    rep = report(result, Pa_MPa)
    Pc = result.Pc_MPa * 1e6
    mdot = Pc * At_m2 / result.c_star_m_s
    F = rep.Cf_at_Pa * Pc * At_m2
    return MotorPoint(mdot, F, F / (mdot * G0), rep.Cf_at_Pa)
