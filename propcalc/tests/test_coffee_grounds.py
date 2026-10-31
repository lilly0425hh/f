"""커피박 DB 항목(ACS Omega 2021, PMC8756566)이 원문 표 입력값에서 재현되는지 검사한다."""

import math
import statistics

import pytest

from propcalc.hess import dHf_from_dry_basis_report
from propcalc.mixture import load_ingredients

M_H2O = 18.01528
W_H, W_O = 2 * 1.00794 / M_H2O, 15.9994 / M_H2O
TABLE = dict(C=50.0, H=6.7, N=2.3, S=0.9, O=39.0, ash=0.9, M=3.3, hhv=22.3)
# PMC9061234 Table 1 의 SCG 행 16개 (수치만): 문헌 간 산포의 근거
LIT_SCG_HHV = [19.61, 23.10, 23.20, 19.30, 24.00, 23.40, 20.20, 21.1, 24.39, 22.83, 22.37, 19.55, 21.14, 24.91, 22.74, 21.15]


def _stated_dry(u):
    r = TABLE
    return dHf_from_dry_basis_report(C=r["C"], H=r["H"], N=r["N"], S=r["S"], O=r["O"], ash=r["ash"],
                                     hhv_dry_MJ_per_kg=r["hhv"], u_hhv_dry_MJ_per_kg=u)


def _moisture_corrected(u):
    r, f = TABLE, 1 - TABLE["M"] / 100
    return dHf_from_dry_basis_report(
        C=r["C"] / f, N=r["N"] / f, S=r["S"] / f, H=(r["H"] - W_H * r["M"]) / f, O=(r["O"] - W_O * r["M"]) / f,
        ash=r["ash"] / f, hhv_dry_MJ_per_kg=r["hhv"] / f, u_hhv_dry_MJ_per_kg=u / f)


def test_literature_spread_values():
    assert len(LIT_SCG_HHV) == 16
    assert statistics.mean(LIT_SCG_HHV) == pytest.approx(22.06, abs=0.01)
    assert statistics.stdev(LIT_SCG_HHV) == pytest.approx(1.80, abs=0.01)


def test_db_entry_matches_table_inputs():
    ing = load_ingredients()["SPENT_COFFEE_GROUNDS"]
    sd = statistics.stdev(LIT_SCG_HHV)
    a, c = _stated_dry(sd), _moisture_corrected(sd)
    assert ing.dHf_kJ_per_mol == pytest.approx(a.dHf_kJ_per_mol, abs=0.1)
    half = abs(a.dHf_kJ_per_mol - c.dHf_kJ_per_mol) / 2
    assert ing.u_dHf_kJ_per_mol == pytest.approx(math.hypot(a.u_dHf_kJ_per_mol, half), abs=0.5)
    assert ing.M_g_per_mol == pytest.approx(a.M_g_per_mol, abs=0.01)
    for el, v in a.formula.items():
        assert ing.formula[el] == pytest.approx(v, abs=1e-3)
    assert ing.verification == "secondary"


def test_uncertainty_is_dominated_by_literature_spread_not_basis():
    ing = load_ingredients()["SPENT_COFFEE_GROUNDS"]
    hhv_part = statistics.stdev(LIT_SCG_HHV) * ing.M_g_per_mol
    assert hhv_part > 0.9 * ing.u_dHf_kJ_per_mol
