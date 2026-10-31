import pytest

from propcalc import sweep
from propcalc.hess import dHf_from_hhv

FORMULA = {"C": 10.0, "H": 14.45, "O": 6.35}


def test_dHf_column_matches_hess_inversion():
    df = sweep.sweep_fuel_hhv(FORMULA, [18.0, 20.0], [65.0])
    for hhv, g in df.groupby("hhv_MJ_per_kg"):
        assert g.dHf_kJ_per_mol.iloc[0] == pytest.approx(dHf_from_hhv(FORMULA, hhv).dHf_kJ_per_mol)


def test_higher_hhv_gives_higher_isp_and_tc_at_fixed_ratio():
    # 일관 구간(상관식 HHV 20.6 ± 2) 안에서 연료 HHV 가 클수록 성능이 커야 한다
    df = sweep.sweep_fuel_hhv(FORMULA, [19.0, 20.0, 21.0, 22.0], [65.0])
    ok = df[df.ok].sort_values("hhv_MJ_per_kg")
    assert len(ok) == 4
    assert ok.isp_s.is_monotonic_increasing
    assert ok.Tc_K.is_monotonic_increasing


def test_rows_carry_assumptions_and_hypothetical_marking():
    df = sweep.sweep_fuel_hhv(FORMULA, [19.0], [60.0, 65.0])
    assert {"hhv_MJ_per_kg", "dHf_kJ_per_mol", "ox_wt_pct", "ok", "flag"} <= set(df.columns)
    text = df.attrs["assumptions"]
    assert "가정형 지도" in text and "회분 제외" in text and "평형 유동" in text
    assert df.fuel.iloc[0].startswith("HYPO")


def test_failed_points_are_kept_not_dropped():
    df = sweep.sweep_fuel_hhv(FORMULA, [16.0], range(40, 86))
    assert len(df) == 46
    bad = df[~df.ok]
    assert len(bad) > 0
    assert bad.isp_s.isna().all()
    assert (bad.error != "").all()
