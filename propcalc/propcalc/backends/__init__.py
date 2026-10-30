from .base import Backend, ConvergenceError, Result

__all__ = ["Backend", "ConvergenceError", "Result", "get_backend"]


def get_backend(name: str = "cea") -> Backend:
    """이름으로 백엔드를 만든다. 'cea' (1순위) 또는 'rocketcea' (교차검증)."""
    if name == "cea":
        from .cea_backend import CEABackend

        return CEABackend()
    if name == "rocketcea":
        from .rocketcea_backend import RocketCEABackend

        return RocketCEABackend()
    raise ValueError(f"알 수 없는 백엔드: {name!r}")


def compare_backends(mixture, Pc_MPa, Pe_MPa=None, area_ratio=None, names=("cea", "rocketcea")):
    """같은 입력을 두 백엔드에 넣고 결과와 상대 차이를 돌려준다 (평형 계산)."""
    results = {n: get_backend(n).run(mixture, Pc_MPa, Pe_MPa=Pe_MPa, area_ratio=area_ratio) for n in names}
    a, b = (results[n] for n in names)
    diff = {
        "dTc_K": a.Tc_K - b.Tc_K,
        "dc_star_rel": (a.c_star_m_s - b.c_star_m_s) / b.c_star_m_s,
        "dIsp_rel": (a.isp_s - b.isp_s) / b.isp_s,
        "dM_eff_rel": (a.M_eff_g_mol - b.M_eff_g_mol) / b.M_eff_g_mol,
        "dcond_mass_frac": a.condensed_mass_fraction - b.condensed_mass_fraction,
    }
    return results, diff
