"""어린 코코넛 외피(RMUTP 2025) DB 항목이 원문 표의 입력값에서 재현되는지 검사한다.

코코피트가 아니라 외피이며, 환산 가정(수분 포함 기준, 수분의 H·O 제거)은 이 저장소의 것이다.
"""

import math

import pytest

from propcalc.hess import dHf_from_dry_basis_report
from propcalc.mixture import load_ingredients

M_H2O = 18.01528
W_H, W_O = 2 * 1.00794 / M_H2O, 15.9994 / M_H2O

# 논문 Table 1(수분, 회분, HHV)과 Table 2(원소, 질량%), 같은 산지 행
TABLE = {
    "YOUNG_COCONUT_HUSK_RB": dict(C=45.58, H=6.13, N=0.34, S=0.07, O=41.83, ash=6.05, M=3.51, hhv=17.52, u=0.26),
    "YOUNG_COCONUT_HUSK_SS": dict(C=45.93, H=5.89, N=0.36, S=0.04, O=42.88, ash=4.90, M=4.73, hhv=18.39, u=0.13),
}


def _report(r, wet_corrected: bool):
    if wet_corrected:
        f = 1 - r["M"] / 100
        return dHf_from_dry_basis_report(
            C=r["C"] / f, N=r["N"] / f, S=r["S"] / f,
            H=(r["H"] - W_H * r["M"]) / f, O=(r["O"] - W_O * r["M"]) / f,
            ash=r["ash"] / f, hhv_dry_MJ_per_kg=r["hhv"] / f, u_hhv_dry_MJ_per_kg=r["u"] / f)
    return dHf_from_dry_basis_report(
        C=r["C"], H=r["H"], N=r["N"], S=r["S"], O=r["O"], ash=r["ash"],
        hhv_dry_MJ_per_kg=r["hhv"], u_hhv_dry_MJ_per_kg=r["u"])


@pytest.mark.parametrize("key", TABLE)
def test_db_entry_matches_table_inputs(key):
    db = load_ingredients()
    ing, r = db[key], TABLE[key]
    c, a = _report(r, True), _report(r, False)
    assert ing.dHf_kJ_per_mol == pytest.approx(c.dHf_kJ_per_mol, abs=0.1)
    half = abs(a.dHf_kJ_per_mol - c.dHf_kJ_per_mol) / 2
    assert ing.u_dHf_kJ_per_mol == pytest.approx(math.hypot(c.u_dHf_kJ_per_mol, half), abs=0.1)
    assert ing.M_g_per_mol == pytest.approx(c.M_g_per_mol, abs=0.01)
    for el, v in c.formula.items():
        assert ing.formula[el] == pytest.approx(v, abs=1e-3)
    assert ing.verification == "secondary"
    assert "코코피트가 아니다" in ing.note


@pytest.mark.parametrize("key", TABLE)
def test_table_closure_with_ash(key):
    # 표의 O 는 회분을 뺀 차이법: C+H+N+S+O+회분 = 100 (식 (2)에는 회분이 없으나 표 값은 그렇다)
    r = TABLE[key]
    assert r["C"] + r["H"] + r["N"] + r["S"] + r["O"] + r["ash"] == pytest.approx(100.0, abs=0.01)


def test_uniform_division_by_moisture_is_rejected():
    # O 가 수분의 산소를 포함하므로 전체를 (1-수분)으로 나누면 합이 100% 를 넘는다. 닫힘 검사가 거부해야 한다.
    r = TABLE["YOUNG_COCONUT_HUSK_RB"]
    f = 1 - r["M"] / 100
    with pytest.raises(ValueError, match="합"):
        dHf_from_dry_basis_report(
            C=r["C"] / f, H=r["H"] / f, N=r["N"] / f, S=r["S"] / f, O=r["O"] / f,
            ash=r["ash"] / f, hhv_dry_MJ_per_kg=r["hhv"] / f)


def test_lot_difference_exceeds_stated_uncertainty():
    db = load_ingredients()
    rb, ss = db["YOUNG_COCONUT_HUSK_RB"], db["YOUNG_COCONUT_HUSK_SS"]
    gap = abs(rb.dHf_kJ_per_mol - ss.dHf_kJ_per_mol)
    assert gap > max(rb.u_dHf_kJ_per_mol, ss.u_dHf_kJ_per_mol)
