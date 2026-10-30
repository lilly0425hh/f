"""5단계: 스윕, 연료 대체 비교, 민감도, 그래프 가정 표기."""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from propcalc import sweep
from propcalc.hess import ingredient_from_hhv
from propcalc.mixture import load_ingredients, stoichiometry

pytest.importorskip("cea")

DB = load_ingredients()
PC, PE = 6.895, 0.101325

# 리그닌 고정값은 DB 에 없다. 아래는 지시서의 참고 조성 + HHV 26.8 로 만든 계산 검증 전용 가상 시료.
LIGNIN = ingredient_from_hhv(
    "LIGNIN_TEST", "연목 크라프트 리그닌 (계산 검증 전용)", "C10H11.5O2.9", 26.8,
    u_hhv_MJ_per_kg=0.5, source="test fixture: 지시서 4절 참고 조성(BioResources 2020), 실제 시료 아님",
)


@pytest.fixture(scope="module")
def knsb_sweep():
    return sweep.sweep_oxidizer("SORBITOL", ox_wt_pct=range(50, 86), Pc_MPa=PC, Pe_MPa=PE, db=DB)


def test_sweep_shape_and_convergence(knsb_sweep):
    df = knsb_sweep
    assert len(df) == 36 and df["ok"].all()
    assert df["ox_wt_pct"].is_monotonic_increasing
    assert {"Tc_K", "c_star_m_s", "isp_s", "cond_mass_frac", "C_gr_mass_frac", "has_C_gr"} <= set(df.columns)


def test_sweep_carries_assumptions(knsb_sweep):
    a = knsb_sweep.attrs["assumptions"]
    for token in ("평형", "Pc=6.895", "백엔드=cea"):
        assert token in a


def test_isp_maximum_is_on_fuel_rich_side(knsb_sweep, capsys):
    """기대: Isp 최대점이 화학량론(74.3%)보다 연료 과잉 쪽. 다르면 원인을 기록한다 (실패 아님)."""
    opt = sweep.isp_optimum(knsb_sweep)
    s = stoichiometry(DB["KNO3"], DB["SORBITOL"])
    with capsys.disabled():
        print(f"\n[KNSB 스윕] Isp 최대 {opt['isp_s']:.2f} s @ KNO3 {opt['ox_wt_pct']:g}% "
              f"(화학량론 {s.mass_pct_oxidizer:.1f}%, 오프셋 {opt['ox_wt_pct'] - s.mass_pct_oxidizer:+.1f}%p), "
              f"c* 최대 @ {knsb_sweep.loc[knsb_sweep.c_star_m_s.idxmax(), 'ox_wt_pct']:g}%, "
              f"Tc 최대 @ {knsb_sweep.loc[knsb_sweep.Tc_K.idxmax(), 'ox_wt_pct']:g}%")
    assert opt["ox_wt_pct"] < s.mass_pct_oxidizer
    assert opt["stoich_ox_wt_pct"] == pytest.approx(s.mass_pct_oxidizer)


def test_carbon_not_present_for_sorbitol_in_range(knsb_sweep):
    assert not knsb_sweep["has_C_gr"].any()


def test_melting_plateau_is_flagged(knsb_sweep):
    """K2CO3 융점 구간에서 γ_s<1 같은 평형 모델 아티팩트가 있으면 표시되어야 한다."""
    bad = knsb_sweep[knsb_sweep["gamma_s"] < 1.0]
    assert len(bad) >= 1
    assert bad["flag"].str.contains("γ_s<1").all()
    # 표시된 점은 연료 과잉 쪽 plateau 통과 구간(50~60%)과 연료 부족 쪽 plateau 통과 구간(84~88%)에만 있어야 한다
    flagged = knsb_sweep.loc[knsb_sweep["flag"] != "", "ox_wt_pct"]
    assert (flagged.between(50, 60) | flagged.between(84, 88)).all()


def test_plateau_neighborhood_is_flagged_where_cstar_kinks(knsb_sweep):
    """c* 가 56~58% 에서 꺾이는 것은 목 온도가 K2CO3 상전이 plateau(≈1173 K) 근처에 있기 때문이다.
    상 공존이 관측된 온도 ±20 K 안의 점은 연소실·목 어느 쪽이든 표시해야 한다."""
    df = knsb_sweep.set_index("ox_wt_pct")
    for ox in (56, 57, 58):
        assert "plateau" in df.loc[ox, "flag"], ox
    for ox in (62, 65, 68, 72):
        assert df.loc[ox, "flag"] == "", ox


def test_lignin_fuel_rich_forms_graphite():
    df = sweep.sweep_oxidizer(LIGNIN, ox_wt_pct=[35, 40, 50, 60, 70, 80, 84, 90], Pc_MPa=PC, Pe_MPa=PE, db=DB)
    assert df["ok"].all()
    rich = df[df["ox_wt_pct"] <= 50]
    assert rich["has_C_gr"].all() and (rich["C_gr_mass_frac"] > 0.01).all()
    assert not df.loc[df["ox_wt_pct"] >= 84, "has_C_gr"].any()
    # 연료 과잉일수록 흑연이 늘어난다
    assert df.loc[df.ox_wt_pct == 35, "C_gr_mass_frac"].iloc[0] > df.loc[df.ox_wt_pct == 50, "C_gr_mass_frac"].iloc[0]


def test_rejected_points_are_recorded_not_dropped():
    """CEA 가 쓰레기 값을 내는 조성은 행을 지우지 않고 ok=False 와 사유를 남기고, 최대점 계산에서 제외한다."""
    df = sweep.sweep_oxidizer("SORBITOL", ox_wt_pct=[72, 87, 90], Pc_MPa=PC, Pe_MPa=PE, db=DB)
    assert list(df["ok"]) == [True, False, False]
    assert df.loc[1, "error"].startswith("비물리적") and np.isnan(df.loc[1, "isp_s"])
    assert sweep.isp_optimum(df)["ox_wt_pct"] == 72


def test_neighbors_of_rejected_points_are_flagged():
    """소르비톨 87% 는 거부되고 88% 는 값이 매끄러워 보이지만 불안정 구간 안이므로 표시해야 한다.
    '이웃'은 스윕 격자 순서상의 바로 앞뒤 점이다."""
    df = sweep.sweep_oxidizer("SORBITOL", ox_wt_pct=[83, 84, 85, 86, 87, 88], Pc_MPa=PC, Pe_MPa=PE, db=DB)
    assert list(df["ok"]) == [True, True, True, True, False, True]
    assert "인접 점 거부" in df.loc[3, "flag"] and "인접 점 거부" in df.loc[5, "flag"]
    assert df.loc[0, "flag"] == "" and df.loc[1, "flag"] == ""


def test_compare_fuels_long_table():
    df = sweep.compare_fuels({"소르비톨": "SORBITOL", "셀룰로오스": "CELLULOSE_COTTON", "리그닌(가상)": LIGNIN},
                             ox_wt_pct=[60, 65, 70], Pc_MPa=PC, Pe_MPa=PE, db=DB)
    assert set(df["fuel"]) == {"소르비톨", "셀룰로오스", "리그닌(가상)"}
    assert len(df) == 9 and df["ok"].all()
    assert (df["stoich_ox_wt_pct"] > 70).all()
    best = sweep.best_by_isp(df)
    assert len(best) == 3


def test_fuel_blend_endpoints_match_pure_fuels():
    blend = sweep.sweep_fuel_blend("SORBITOL", "CELLULOSE_COTTON", fraction_b_pct=[0, 25, 50, 100],
                                   ox_wt_pct=65.0, Pc_MPa=PC, Pe_MPa=PE, db=DB)
    pure_a = sweep.sweep_oxidizer("SORBITOL", ox_wt_pct=[65], Pc_MPa=PC, Pe_MPa=PE, db=DB)
    pure_b = sweep.sweep_oxidizer("CELLULOSE_COTTON", ox_wt_pct=[65], Pc_MPa=PC, Pe_MPa=PE, db=DB)
    assert blend["isp_s"].iloc[0] == pytest.approx(pure_a["isp_s"].iloc[0], rel=1e-9)
    assert blend["isp_s"].iloc[-1] == pytest.approx(pure_b["isp_s"].iloc[0], rel=1e-9)
    assert len(blend) == 4 and blend["ok"].all()


def test_sensitivity_bounds_bracket_nominal():
    from propcalc.mixture import Mixture

    m = Mixture.from_wt_pct({"KNO3": 65, "CELLULOSE_COTTON": 35}, DB)
    s = sweep.sensitivity_dHf(m, "CELLULOSE_COTTON", Pc_MPa=PC, Pe_MPa=PE)  # u = ±30 kJ/mol (DB)
    t = s.table
    assert t.loc["low", "Tc_K"] < t.loc["nominal", "Tc_K"] < t.loc["high", "Tc_K"]
    assert t.loc["low", "isp_s"] < t.loc["nominal", "isp_s"] < t.loc["high", "isp_s"]
    assert s.u_kJ_per_mol == 30.0
    # ±30 kJ/mol 의 영향이 0 이 아닌 의미 있는 크기여야 한다
    assert t.loc["high", "Tc_K"] - t.loc["low", "Tc_K"] > 10


def test_sensitivity_monte_carlo_is_reproducible_and_consistent_with_bounds():
    from propcalc.mixture import Mixture

    m = Mixture.from_wt_pct({"KNO3": 65, "CELLULOSE_COTTON": 35}, DB)
    a = sweep.sensitivity_dHf(m, "CELLULOSE_COTTON", method="mc", n=60, seed=7, Pc_MPa=PC, Pe_MPa=PE)
    b = sweep.sensitivity_dHf(m, "CELLULOSE_COTTON", method="mc", n=60, seed=7, Pc_MPa=PC, Pe_MPa=PE)
    assert a.table.equals(b.table)
    bd = sweep.sensitivity_dHf(m, "CELLULOSE_COTTON", Pc_MPa=PC, Pe_MPa=PE).table
    half_range = (bd.loc["high", "Tc_K"] - bd.loc["low", "Tc_K"]) / 2
    # u 를 1σ 로 해석한 MC 표준편차 ≈ 선형 응답 × u = 경계 반폭
    assert a.table.loc["std", "Tc_K"] == pytest.approx(half_range, rel=0.35)
    assert a.n_failed == 0


def test_sensitivity_requires_known_uncertainty():
    from propcalc.mixture import Mixture

    m = Mixture.from_wt_pct({"KNO3": 65, "SORBITOL": 35}, DB)
    with pytest.raises(ValueError, match="불확실도"):
        sweep.sensitivity_dHf(m, "KNO3", Pc_MPa=PC, Pe_MPa=PE)  # KNO3 u 는 미확인(null)
    ok = sweep.sensitivity_dHf(m, "KNO3", u_kJ_per_mol=0.6, Pc_MPa=PC, Pe_MPa=PE)  # 명시하면 가능
    assert ok.table.loc["high", "Tc_K"] > ok.table.loc["low", "Tc_K"]


def test_plot_contains_assumptions(knsb_sweep, tmp_path):
    fig = sweep.plot_sweep(knsb_sweep, path=tmp_path / "knsb.png")
    texts = " ".join(t.get_text() for t in fig.texts)
    assert "Pc=6.895" in texts and "평형" in texts and "cea" in texts
    assert (tmp_path / "knsb.png").stat().st_size > 10_000
    assert (tmp_path / "knsb.csv").exists()  # 표 보기(접근성)용
