"""5단계 분석: KNSB 스윕, 연료 대체 비교, ΔfH° 민감도. 그림은 notebooks/figures/ 에 저장한다.

    python notebooks/stage5_analysis.py

리그닌은 정의된 분자가 아니므로 DB 에 고정값이 없다. 여기서는 지시서의 참고 조성(C10H11.5O2.9)과
가정 HHV 26.8±0.5 MJ/kg 로 만든 '가상 시료'를 쓴다. 실제 시료의 원소분석·HHV 로 교체해야 의미가 있다.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from propcalc import sweep
from propcalc.hess import ingredient_from_hhv
from propcalc.mixture import Mixture, load_ingredients

OUT = Path(__file__).parent / "figures"
PC, PE = 6.895, 0.101325
DB = load_ingredients()
LIGNIN = ingredient_from_hhv(
    "LIGNIN_TEST", "리그닌(가상 시료)", "C10H11.5O2.9", 26.8, u_hhv_MJ_per_kg=0.5,
    source="가상 시료: 지시서 4절 참고 조성(BioResources 2020) + 가정 HHV. 실제 시료로 교체할 것")
FUELS = {"소르비톨": "SORBITOL", "셀룰로오스(면)": "CELLULOSE_COTTON", "리그닌(가상 시료)": LIGNIN}


def fig_knsb():
    df = sweep.sweep_oxidizer("SORBITOL", ox_wt_pct=np.arange(50, 86, 1.0), Pc_MPa=PC, Pe_MPa=PE, db=DB)
    sweep.plot_sweep(df, OUT / "knsb_sweep.png", title="KNO3/소르비톨 산화제 비율 스윕")
    o = sweep.isp_optimum(df)
    print(f"[KNSB] Isp 최대 {o['isp_s']:.2f} s @ {o['ox_wt_pct']:g}% (화학량론 {o['stoich_ox_wt_pct']:.2f}%), "
          f"Tc={o['Tc_K']:.0f} K, M_eff={o['M_eff_g_mol']:.1f}, 응축상 {o['cond_mass_frac']:.3f}")
    for ox in (o["ox_wt_pct"], round(o["stoich_ox_wt_pct"])):
        r = df[df.ox_wt_pct == ox].iloc[0]
        print(f"   ox {ox:g}%: Tc={r.Tc_K:.1f} M_eff={r.M_eff_g_mol:.2f} cond={r.cond_mass_frac:.4f} Isp={r.isp_s:.3f}")
    print("   flag 구간:", df.loc[df.flag != "", ["ox_wt_pct", "Tc_K", "Tt_K", "gamma_s", "flag"]].round(3).to_dict("records"))
    return df


def fig_fuels():
    grid = np.arange(30, 91, 1.0)
    df = sweep.compare_fuels(FUELS, grid, Pc_MPa=PC, Pe_MPa=PE, db=DB)
    df.to_csv(OUT / "fuel_comparison.csv", index=False)
    colors = dict(zip(FUELS, sweep.SERIES_COLORS))
    markers = dict(zip(FUELS, ("o", "s", "^")))
    sweep.apply_font()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5), facecolor=sweep.SURFACE)
    for ax, (col, lab) in zip(axes, (("isp_s", "비추력 Isp (s, Pa=Pe)"), ("Tc_K", "연소실 온도 Tc (K)"))):
        sweep._style(ax)
        for name in FUELS:
            d = df[(df.fuel == name) & df.ok]
            ax.plot(d.ox_wt_pct, d[col], color=colors[name], linewidth=2, label=name, marker=markers[name],
                    markevery=10, markersize=7, markeredgecolor=sweep.SURFACE, markeredgewidth=1.5)
            st = d.stoich_ox_wt_pct.iloc[0]
            ax.axvline(st, color=colors[name], linewidth=1, linestyle=":")
            fl = d[d.flag.fillna("") != ""]
            ax.plot(fl.ox_wt_pct, fl[col], linestyle="none", marker="o", markersize=4, markerfacecolor=sweep.SURFACE,
                    markeredgecolor=sweep.INK_2, markeredgewidth=0.9, zorder=3)
            if col == "isp_s":
                b = d.loc[d.isp_s.idxmax()]
                ax.plot(b.ox_wt_pct, b.isp_s, marker=markers[name], markersize=10, color=colors[name],
                        markeredgecolor=sweep.INK, markeredgewidth=1.5)
        ax.set_xlabel("KNO3 질량% (나머지 연료)", fontsize=9)
        ax.set_ylabel(lab, fontsize=10)
        ax.legend(frameon=False, fontsize=9, labelcolor=sweep.INK, loc="lower left" if col == "isp_s" else "upper left")
    fig.suptitle("연료 대체 비교 (점선 = 각 연료의 화학량론 비율, 큰 마커 = Isp 최대)", x=0.05, ha="left", fontsize=12)
    rej = df[~df.ok].groupby("fuel").ox_wt_pct.apply(lambda s: f"{int(s.min())}~{int(s.max())}%").to_dict()
    note = "제외(CEA 비물리적 결과): " + ", ".join(f"{k} KNO3 {v}" for k, v in rej.items()) if rej else ""
    fig.text(0.05, 0.012, "가정: " + df.attrs["assumptions"] + " | 리그닌은 가상 시료(검증 전용)\n○ = 상전이 plateau/상 공존 의심 점. " + note,
             fontsize=8.5, color=sweep.INK_2, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    fig.savefig(OUT / "fuel_comparison.png", dpi=150, facecolor=sweep.SURFACE)
    best = sweep.best_by_isp(df)
    print("[연료 비교] Isp 최대점:")
    print(best[["fuel", "ox_wt_pct", "stoich_ox_wt_pct", "isp_s", "Tc_K", "c_star_m_s", "cond_mass_frac",
                "C_gr_mass_frac"]].round(3).to_string(index=False))
    cg = df[(df.fuel == "리그닌(가상 시료)") & df.has_C_gr]
    print(f"[리그닌] C(gr) 생성 구간: KNO3 {cg.ox_wt_pct.min():g}~{cg.ox_wt_pct.max():g}% (최대 질량분율 {cg.C_gr_mass_frac.max():.3f})")
    for n in ("소르비톨", "셀룰로오스(면)"):
        d = df[(df.fuel == n) & df.has_C_gr]
        print(f"[{n}] C(gr) 생성 구간:", "없음" if d.empty else f"{d.ox_wt_pct.min():g}~{d.ox_wt_pct.max():g}%")
    for n in FUELS:
        bad = df[(df.fuel == n) & (df.flag != "")]
        if len(bad):
            print(f"[{n}] 모델 아티팩트 의심 점: {bad.ox_wt_pct.tolist()}")


def fig_sensitivity():
    rows = []
    cases = [("셀룰로오스(면)", "CELLULOSE_COTTON", None, 65.0), ("리그닌(가상 시료)", LIGNIN, None, 65.0)]
    for label, fuel, u, ox in cases:
        ing = LIGNIN if fuel is LIGNIN else DB[fuel]
        m = Mixture.from_wt_pct({"KNO3": ox, ing.key: 100 - ox}, DB, extra=[ing])
        s = sweep.sensitivity_dHf(m, ing.key, u_kJ_per_mol=u, Pc_MPa=PC, Pe_MPa=PE)
        mc = sweep.sensitivity_dHf(m, ing.key, u_kJ_per_mol=u, method="mc", n=300, seed=1, Pc_MPa=PC, Pe_MPa=PE)
        t = s.table
        print(f"[민감도] {label} ΔfH°={s.dHf_kJ_per_mol:.1f}±{s.u_kJ_per_mol:.1f} kJ/mol, KNO3 {ox:g}%: "
              f"Tc {t.loc['low','Tc_K']:.0f}~{t.loc['high','Tc_K']:.0f} K (공칭 {t.loc['nominal','Tc_K']:.0f}), "
              f"Isp {t.loc['low','isp_s']:.2f}~{t.loc['high','isp_s']:.2f} s (공칭 {t.loc['nominal','isp_s']:.2f}); "
              f"MC σ(Tc)={mc.table.loc['std','Tc_K']:.1f} K, σ(Isp)={mc.table.loc['std','isp_s']:.3f} s")
        rows.append((label, s))
    # 소르비톨 대비 비교: 불확실도가 작은 소르비톨(±1.4)
    m = Mixture.from_wt_pct({"KNO3": 65, "SORBITOL": 35}, DB)
    s = sweep.sensitivity_dHf(m, "SORBITOL", Pc_MPa=PC, Pe_MPa=PE)
    t = s.table
    print(f"[민감도] 소르비톨 ΔfH°=±{s.u_kJ_per_mol:g}: Tc {t.loc['low','Tc_K']:.1f}~{t.loc['high','Tc_K']:.1f} K, "
          f"Isp {t.loc['low','isp_s']:.3f}~{t.loc['high','isp_s']:.3f} s")
    rows.insert(0, ("소르비톨", s))
    sweep.apply_font()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4), facecolor=sweep.SURFACE, sharey=True)
    for ax, (col, lab, unit) in zip(axes, (("Tc_K", "Tc 변화 (K)", "K"), ("isp_s", "Isp 변화 (s)", "s"))):
        sweep._style(ax)
        ax.grid(axis="x", color=sweep.GRID, linewidth=0.8)
        ax.grid(axis="y", visible=False)
        ax.axvline(0, color=sweep.INK_2, linewidth=1)
        for i, (label, s) in enumerate(rows):
            nom = s.table.loc["nominal", col]
            lo, hi = s.table.loc["low", col] - nom, s.table.loc["high", col] - nom
            ax.plot([lo, hi], [i, i], color=sweep.SERIES_COLORS[i], linewidth=4, solid_capstyle="round")
            ax.plot([lo, hi], [i, i], linestyle="none", marker="|", markersize=12, color=sweep.SERIES_COLORS[i])
            ax.annotate(f"공칭 {nom:.{0 if unit == 'K' else 1}f} {unit}, ΔfH° ±{s.u_kJ_per_mol:.3g} kJ/mol",
                        (lo, i), xytext=(0, 9), textcoords="offset points", va="bottom", ha="left", fontsize=9, color=sweep.INK,
                        bbox=dict(fc=sweep.SURFACE, ec="none", pad=1.5))
        ax.set_yticks(range(len(rows)), [r[0] for r in rows], color=sweep.INK, fontsize=10)
        ax.set_xlabel(lab + " (공칭 대비, ΔfH° 상·하한)", fontsize=9)
        ax.set_xlim(-1.35 * max(abs(s.table.loc[k, col] - s.table.loc["nominal", col]) for _, s in rows for k in ("low", "high")),
                    1.9 * max(abs(s.table.loc[k, col] - s.table.loc["nominal", col]) for _, s in rows for k in ("low", "high")))
        ax.set_ylim(len(rows) - 0.5, -0.5)
    fig.suptitle("연료 ΔfH° 불확실도가 Tc·Isp 에 미치는 영향 (공칭 대비 편차)", x=0.04, ha="left", fontsize=12)
    fig.text(0.04, 0.012, "가정: " + sweep._assumption_text("cea", "equilibrium", PC, PE, "KNO3 65% + 연료 35%, ΔfH° 상·하한 2점 계산"),
             fontsize=8.5, color=sweep.INK_2, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    fig.savefig(OUT / "sensitivity.png", dpi=150, facecolor=sweep.SURFACE)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_knsb()
    fig_fuels()
    fig_sensitivity()
    print("figures ->", OUT)
