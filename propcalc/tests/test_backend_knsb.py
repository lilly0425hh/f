"""3단계 합격 기준: KNSB 65/35, Pc=6.895 MPa, Pe=0.101325 MPa, 평형 유동.

문헌값은 아마추어 자료(Nakka)와 학술지 값이 섞여 있으므로 범위 비교만 한다.
"""

import json
from pathlib import Path

import pytest

from propcalc.backends import ConvergenceError, compare_backends, get_backend
from propcalc.mixture import Mixture, load_ingredients
from propcalc.reference import compare_to_nakka

pytest.importorskip("cea")

PC, PE = 6.895, 0.101325
REF = Path(__file__).parent / "reference" / "knsb_65_35.json"


@pytest.fixture(scope="module")
def knsb():
    return Mixture.from_wt_pct({"KNO3": 65, "SORBITOL": 35}, load_ingredients())


@pytest.fixture(scope="module")
def res(knsb):
    return get_backend("cea").run(knsb, PC, Pe_MPa=PE)


def test_chamber_temperature(res):
    assert 1550 <= res.Tc_K <= 1650  # 문헌 약 1600 K


def test_c_star(res):
    assert 880 <= res.c_star_m_s <= 950  # 908 (PROPEP/Nakka), 938 (J. Chil. Chem. Soc. 2023)


def test_condensed_mass_fraction(res):
    assert 0.40 <= res.condensed_mass_fraction <= 0.47  # 문헌 0.436
    assert any(k.startswith("K2CO3") for k in res.condensed_mass_fractions)


def test_effective_molecular_weight(res):
    assert 37 <= res.M_eff_g_mol <= 42  # 문헌 39.9


def test_isp(res):
    assert 145 <= res.isp_s <= 170      # 151.5 s, 164 s
    assert 145 <= res.isp_vac_s <= 170


def test_fuel_rich_products_include_co_and_h2(res):
    gas = res.gas_mole_fractions
    assert gas["CO"] > 0.05 and gas["H2"] > 0.05


def test_assumptions_are_reported(res):
    s = res.assumptions()
    for token in ("평형", "Pc=6.895", "백엔드=cea", "KNO3 65"):
        assert token in s


def test_provenance_warning_is_attached(res):
    # KNO3, 소르비톨 ΔfH° 는 1차 문헌 원문 미확인이므로 결과에 경고가 실려야 한다
    assert any("KNO3" in w for w in res.warnings)


def test_nakka_reference_equation_comparison(res, knsb, capsys):
    """어긋나도 실패로 처리하지 않는다 (DB 차이는 기록만). 표만 출력한다."""
    table = compare_to_nakka(res, knsb)
    with capsys.disabled():
        print("\n[Nakka 평형식(6.8 MPa) 대비 소르비톨 1 mol 당 생성물 mol]\n" + table)
    assert "K2CO3" in table


def test_builtin_gas_phase_kno3_would_be_wrong(knsb):
    """회귀 방지: CEA 내장 'KNO3' 는 기체상(약 -316 kJ/mol)이라 쓰면 Tc 가 400 K 가량 높게 나온다."""
    import cea
    import numpy as np

    gas = cea.Mixture(["KNO3"]).calc_property(cea.ENTHALPY, np.array([1.0]), np.array([298.15]))
    per_mol = gas * 101.1032 / 1e6  # J/kg -> kJ/mol
    assert -330 < per_mol < -300  # 고체(-494.6)가 아님
    assert knsb.ingredient("KNO3").dHf_kJ_per_mol < -490


@pytest.mark.needs_rocketcea
def test_backends_agree(knsb):
    pytest.importorskip("rocketcea")
    results, diff = compare_backends(knsb, PC, Pe_MPa=PE)
    assert abs(diff["dTc_K"]) <= 5
    assert abs(diff["dc_star_rel"]) <= 0.01
    assert abs(diff["dIsp_rel"]) <= 0.01
    assert abs(diff["dcond_mass_frac"]) <= 0.01


def test_area_ratio_input_matches_pressure_input(knsb, res):
    r2 = get_backend("cea").run(knsb, PC, area_ratio=res.area_ratio)
    assert r2.Pe_MPa == pytest.approx(PE, rel=1e-4)
    assert r2.isp_s == pytest.approx(res.isp_s, rel=1e-4)


def test_frozen_flow_with_condensed_phase_is_refused(knsb):
    """CEA 동결 유동은 응축상이 있는 KNSB 에서 중단된다. 조용히 넘어가면 안 된다."""
    with pytest.raises(ConvergenceError, match="동결"):
        get_backend("cea").run(knsb, PC, Pe_MPa=PE, mode="frozen")


@pytest.mark.parametrize("kwargs", [
    dict(Pe_MPa=0.1, area_ratio=8.0),   # 둘 다 지정
    dict(),                             # 둘 다 없음
    dict(Pe_MPa=7.0),                   # Pe > Pc
    dict(area_ratio=0.5),
])
def test_invalid_expansion_arguments(knsb, kwargs):
    with pytest.raises(ValueError):
        get_backend("cea").run(knsb, PC, **kwargs)


def test_invalid_mode(knsb):
    with pytest.raises(ValueError):
        get_backend("cea").run(knsb, PC, Pe_MPa=PE, mode="nonsense")


def test_regression_against_stored_cea_result(res):
    """저장된 CEA 결과(회귀 기준값) 와 일치해야 한다. 기준 파일은 범위 검사 통과 후 생성했다."""
    ref = json.loads(REF.read_text(encoding="utf-8"))
    out = ref["outputs"]
    assert res.Tc_K == pytest.approx(out["Tc_K"], rel=1e-4)
    assert res.c_star_m_s == pytest.approx(out["c_star_m_s"], rel=1e-4)
    assert res.isp_s == pytest.approx(out["isp_s"], rel=1e-4)
    assert res.isp_vac_s == pytest.approx(out["isp_vac_s"], rel=1e-4)
    assert res.M_eff_g_mol == pytest.approx(out["M_eff_g_mol"], rel=1e-4)
    assert res.condensed_mass_fraction == pytest.approx(out["condensed_mass_fraction"], rel=1e-3)
    for sp, val in out["species_moles_per_kg_top"].items():
        assert res.species_moles_per_kg[sp] == pytest.approx(val, rel=1e-3)
