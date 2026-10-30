import numpy as np
import pytest

from propcalc import performance as perf
from propcalc.backends import get_backend
from propcalc.mixture import Mixture, load_ingredients

pytest.importorskip("cea")


@pytest.fixture(scope="module")
def knsb_result():
    m = Mixture.from_wt_pct({"KNO3": 65, "SORBITOL": 35}, load_ingredients())
    return get_backend("cea").run(m, 6.895, Pe_MPa=0.101325)


def test_identities_for_random_inputs():
    rng = np.random.default_rng(12345)
    for _ in range(200):
        Pc = rng.uniform(1e6, 1e7)
        At = rng.uniform(1e-5, 1e-3)
        mdot = rng.uniform(0.01, 5.0)
        ve = rng.uniform(1000, 2500)
        Pe = rng.uniform(5e4, 3e5)
        Pa = rng.uniform(0, 2e5)
        Ae = At * rng.uniform(2, 15)

        cs = perf.c_star(Pc, At, mdot)
        F = perf.thrust(mdot, ve, Pe, Pa, Ae)
        cf = perf.thrust_coefficient(ve, cs, Pe, Pa, Ae, Pc, At)
        isp = perf.isp_from_thrust(F, mdot)

        assert cs == pytest.approx(Pc * At / mdot)
        assert cf == pytest.approx(F / (Pc * At), rel=1e-12)       # Cf = F/(Pc At)
        assert isp == pytest.approx(cf * cs / perf.G0, rel=1e-12)  # Isp = Cf c*/g0


def test_ideal_gas_c_star_air_reference():
    # γ=1.4, M=28.97 g/mol, T=300 K: c* ≈ 428.5 m/s
    assert perf.ideal_gas_c_star(300.0, 28.97, 1.4) == pytest.approx(428.5, abs=0.5)


def test_ideal_gas_c_star_scalings():
    base = perf.ideal_gas_c_star(1600.0, 40.0, 1.13)
    assert perf.ideal_gas_c_star(4 * 1600.0, 40.0, 1.13) == pytest.approx(2 * base)
    assert perf.ideal_gas_c_star(1600.0, 4 * 40.0, 1.13) == pytest.approx(base / 2)


def test_knsb_analytic_vs_equilibrium(knsb_result, capsys):
    rep = perf.report(knsb_result)
    with capsys.disabled():
        print("\n" + rep.summary())
    # CEA 의 M_eff(시스템 질량/기체 몰수)와 γ_s 를 쓰면 해석해가 평형 c* 와 일치한다
    assert abs(rep.c_star_analytic_Meff_rel_diff) < 0.005
    # 응축상을 한 분자로 세는 MW 를 쓰면 c* 가 수 % 과대평가된다
    assert 0.03 < rep.c_star_analytic_MW_rel_diff < 0.09
    assert "응축상" in rep.note


def test_pressure_consistency_with_cea_outputs(knsb_result):
    r = knsb_result
    rep = perf.report(r, Pa_MPa=r.Pe_MPa)
    assert rep.Cf_at_Pa == pytest.approx(r.Cf, rel=1e-4)
    assert rep.isp_at_Pa_s == pytest.approx(r.isp_s, rel=1e-4)
    vac = perf.report(r, Pa_MPa=0.0)
    assert vac.isp_at_Pa_s == pytest.approx(r.isp_vac_s, rel=1e-4)
    sl = perf.report(r, Pa_MPa=0.101325)
    assert sl.isp_at_Pa_s == pytest.approx(r.isp_s, rel=1e-3)  # Pe 가 정확히 1 atm 인 설계


def test_overexpanded_nozzle_loses_isp(knsb_result):
    hi = perf.report(knsb_result, Pa_MPa=0.2).isp_at_Pa_s
    lo = perf.report(knsb_result, Pa_MPa=0.05).isp_at_Pa_s
    assert lo > hi


def test_thrust_scaling_from_throat_area(knsb_result):
    At = 2.0e-5
    pt = perf.motor_point(knsb_result, At_m2=At, Pa_MPa=0.101325)
    assert pt.mdot_kg_s == pytest.approx(6.895e6 * At / knsb_result.c_star_m_s)
    assert pt.thrust_N == pytest.approx(pt.mdot_kg_s * pt.isp_s * perf.G0, rel=1e-12)


def test_isp_is_flagged_as_upper_bound(knsb_result):
    rep = perf.report(knsb_result)
    assert "상한" in rep.upper_bound_note and "입자-기체" in rep.upper_bound_note


def test_empirical_efficiency_requires_explicit_opt_in(knsb_result):
    rep = perf.report(knsb_result)
    assert rep.isp_empirical_s is None
    rep2 = perf.report(knsb_result, eta_isp=perf.ETA_ISP_KNSB_AMATEUR)
    assert rep2.isp_empirical_s == pytest.approx(rep2.isp_at_Pa_s * perf.ETA_ISP_KNSB_AMATEUR, rel=1e-12)
    assert "아마추어" in rep2.empirical_note


@pytest.mark.parametrize("eta", [0.0, -0.1, 1.2])
def test_efficiency_bounds(knsb_result, eta):
    with pytest.raises(ValueError):
        perf.report(knsb_result, eta_isp=eta)


def test_invalid_physical_inputs():
    with pytest.raises(ValueError):
        perf.c_star(1e6, 0.0, 1.0)
    with pytest.raises(ValueError):
        perf.ideal_gas_c_star(1500, 40, 1.0)  # γ <= 1
