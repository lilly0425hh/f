import pytest

from propcalc.formula import molar_mass
from propcalc.hess import (
    HHV_CORRELATION_WARNING,
    channiwala_parikh_hhv,
    dHf_from_elemental_analysis,
    dHf_from_hhv,
    ingredient_from_hhv,
)

CELLULOSE_EA = dict(C=44.45, H=6.22, O=49.34)
LIGNIN_EA = dict(C=67.4, H=6.5, O=26.1)


def test_cellulose_measured_hhv():
    r = dHf_from_hhv("C6H10O5", 17.43)
    assert r.dHf_kJ_per_mol == pytest.approx(-964.0, abs=3)
    assert r.hhv_source == "measured"
    assert not r.warnings


def test_cellulose_correlation_warns():
    with pytest.warns(UserWarning, match="HHV 1~2% 오차"):
        r = dHf_from_elemental_analysis(**CELLULOSE_EA, formula="C6H10O5")
    assert r.hhv_MJ_per_kg == pytest.approx(17.75, abs=0.05)
    assert r.dHf_kJ_per_mol == pytest.approx(-912.0, abs=5)
    assert HHV_CORRELATION_WARNING in r.warnings
    assert r.hhv_source == "correlation"


def test_channiwala_parikh_formula_terms():
    assert channiwala_parikh_hhv(C=100, H=0, O=0) == pytest.approx(34.91)
    assert channiwala_parikh_hhv(C=50, H=6, O=40, N=1, S=1, ash=2) == pytest.approx(
        0.3491 * 50 + 1.1783 * 6 + 0.1005 * 1 - 0.1034 * 40 - 0.0151 * 1 - 0.0211 * 2
    )


def test_lignin_correlation_and_measured_with_given_formula():
    with pytest.warns(UserWarning):
        r_corr = dHf_from_elemental_analysis(**LIGNIN_EA, formula="C10H11.5O2.9")
    assert r_corr.dHf_kJ_per_mol == pytest.approx(-503.0, abs=10)
    r_meas = dHf_from_elemental_analysis(**LIGNIN_EA, hhv_MJ_per_kg=26.8, formula="C10H11.5O2.9")
    assert r_meas.dHf_kJ_per_mol == pytest.approx(-806.0, abs=10)


def test_lignin_formula_derived_from_elemental_analysis():
    r = dHf_from_elemental_analysis(**LIGNIN_EA, hhv_MJ_per_kg=26.8)  # C10 기준 formula 자동 생성
    assert r.formula["C"] == pytest.approx(10.0)
    assert r.formula["H"] == pytest.approx(11.5, abs=0.1)
    assert r.formula["O"] == pytest.approx(2.9, abs=0.05)
    assert r.dHf_kJ_per_mol == pytest.approx(-806.0, abs=10)


def test_error_amplification_is_large():
    """HHV 6% 차이(26.8 vs 28.5)가 ΔfH° 에서 약 300 kJ/mol 차이를 만든다."""
    with pytest.warns(UserWarning):
        corr = dHf_from_elemental_analysis(**LIGNIN_EA, formula="C10H11.5O2.9")
    meas = dHf_from_elemental_analysis(**LIGNIN_EA, hhv_MJ_per_kg=26.8, formula="C10H11.5O2.9")
    assert abs(corr.dHf_kJ_per_mol - meas.dHf_kJ_per_mol) > 250


def test_uncertainty_propagation_from_hhv():
    r = dHf_from_hhv("C6H10O5", 17.43, u_hhv_MJ_per_kg=0.1)
    assert r.u_dHf_kJ_per_mol == pytest.approx(0.1 * 162.14, rel=1e-3)


def test_sulfur_and_nitrogen_products():
    # 황은 SO2(g), 질소는 N2 로 가정. 질소는 생성물 항에 기여하지 않고 M 만 키운다.
    with_n = dHf_from_hhv({"C": 6, "H": 10, "O": 5, "N": 1}, 17.43)
    assert with_n.dHf_kJ_per_mol == pytest.approx(
        6 * -393.51 + 5 * -285.83 + 17.43 * molar_mass({"C": 6, "H": 10, "O": 5, "N": 1}), abs=1e-6)
    with_s = dHf_from_hhv({"C": 6, "H": 10, "O": 5, "S": 1}, 17.43)
    assert with_s.dHf_kJ_per_mol == pytest.approx(
        6 * -393.51 + 5 * -285.83 - 296.81 + 17.43 * molar_mass({"C": 6, "H": 10, "O": 5, "S": 1}), abs=1e-6)


def test_dHf_per_kg_is_independent_of_formula_basis():
    """kJ/mol 값은 화학식 단위에 따라 달라지지만 kJ/kg 은 불변이다."""
    a = dHf_from_hhv("C6H10O5", 17.43)
    b = dHf_from_hhv("C12H20O10", 17.43)
    assert b.dHf_kJ_per_mol == pytest.approx(2 * a.dHf_kJ_per_mol, rel=1e-9)
    assert b.dHf_kJ_per_mol / b.M_g_per_mol == pytest.approx(a.dHf_kJ_per_mol / a.M_g_per_mol, rel=1e-9)


def test_input_validation():
    with pytest.raises(ValueError):
        dHf_from_hhv("C6H10O5", -1.0)
    with pytest.raises(ValueError):
        dHf_from_elemental_analysis(C=44, H=6, O=20)  # 합이 100 이 아님 (회분 미지정)


def test_ingredient_from_hhv_carries_provenance():
    ing = ingredient_from_hhv("LIGNIN_TEST", "lignin (test)", "C10H11.5O2.9", 26.8,
                              u_hhv_MJ_per_kg=0.5, source="test fixture")
    assert ing.dHf_kJ_per_mol == pytest.approx(-806.0, abs=10)
    assert ing.u_dHf_kJ_per_mol == pytest.approx(0.5 * 178.1, rel=0.01)
    assert "HHV" in ing.source and "test fixture" in ing.source


# --- 건조 기준(회분 포함) 보고서 접수 -------------------------------------------------
# 아래 숫자는 실제 코코피트 데이터가 아니라 변환 로직 검증용 합성 데이터(셀룰로오스 조성)이다.
from propcalc.hess import dHf_from_dry_basis_report

_DAF = dict(C=44.44, H=6.22, N=0.0, S=0.0, O=49.34)  # 건조·무회분 기준, 합 100.0


def _dry(ash_pct, hhv_daf=17.43):
    k = 1 - ash_pct / 100
    return dict(C=_DAF["C"] * k, H=_DAF["H"] * k, N=0.0, S=0.0, O=_DAF["O"] * k,
                ash=ash_pct, hhv_dry_MJ_per_kg=hhv_daf * k)


@pytest.mark.parametrize("ash", [0.0, 4.9, 9.1, 15.0])
def test_dry_basis_report_reproduces_daf_result(ash):
    direct = dHf_from_hhv("C6H10O5", 17.43)
    r = dHf_from_dry_basis_report(**_dry(ash), formula="C6H10O5")
    assert r.dHf_kJ_per_mol == pytest.approx(direct.dHf_kJ_per_mol, abs=1e-6)
    assert r.hhv_MJ_per_kg == pytest.approx(17.43, abs=1e-9)  # 무회분 기준 HHV 로 환산됨


def test_oxygen_by_difference_matches_direct_oxygen_and_is_flagged():
    d = _dry(9.1)
    with_o = dHf_from_dry_basis_report(**d, formula="C6H10O5")
    d2 = {k: v for k, v in d.items() if k != "O"}
    by_diff = dHf_from_dry_basis_report(**d2, formula="C6H10O5")
    assert by_diff.dHf_kJ_per_mol == pytest.approx(with_o.dHf_kJ_per_mol, abs=1e-6)
    assert any("차이법" in w for w in by_diff.warnings)
    assert not any("차이법" in w for w in with_o.warnings)


def test_ash_exclusion_is_always_reported():
    r = dHf_from_dry_basis_report(**_dry(9.1), formula="C6H10O5")
    assert any("회분" in w and "9.1" in w for w in r.warnings)


def test_uncertainty_is_scaled_to_ash_free_basis():
    r = dHf_from_dry_basis_report(**_dry(10.0), formula="C6H10O5", u_hhv_dry_MJ_per_kg=0.18)
    assert r.u_dHf_kJ_per_mol == pytest.approx(0.18 / 0.9 * 162.14, rel=1e-3)


def test_closure_failure_is_rejected():
    d = _dry(9.1)
    d["O"] -= 6.0  # 합이 94% 정도 -> 분석 누락/오류 의심
    with pytest.raises(ValueError, match="합"):
        dHf_from_dry_basis_report(**d, formula="C6H10O5")


def test_negative_oxygen_by_difference_is_rejected():
    d = {k: v for k, v in _dry(9.1).items() if k != "O"}
    d["C"] += 50.0  # C+H+N+S+회분 = 105% -> 차이법 산소 음수
    with pytest.raises(ValueError, match="산소"):
        dHf_from_dry_basis_report(**d, formula="C6H10O5")


@pytest.mark.parametrize("ash", [-1.0, 100.0, 120.0])
def test_invalid_ash(ash):
    d = _dry(9.1)
    d["ash"] = ash
    with pytest.raises(ValueError):
        dHf_from_dry_basis_report(**d)
