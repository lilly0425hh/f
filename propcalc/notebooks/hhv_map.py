"""가정형 HHV 지도: 연료 HHV(건조·무회분 기준)가 Tc, c*, Isp 를 얼마나 움직이는가.

    python notebooks/hhv_map.py

코코피트에 대한 결과가 아니다. 코코피트의 같은 시료 원소분석·실측 HHV 를 원문에서 확인한 문헌이 없어,
조성을 고정한 가상 연료의 HHV 만 격자로 바꾼다. 실제 코코피트(또는 어떤 연료)의 HHV 가 확보되면
HHV_daf = HHV_건조 / (1 − 회분 분율) 로 환산해 이 지도의 가로축에서 읽는다.
참고 점(셀룰로오스, 어린 코코넛 외피 두 시료)은 코코피트가 아니며 위치 비교용이다.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from propcalc import sweep
from propcalc.formula import molar_mass
from propcalc.hess import DHF_CO2_G, DHF_H2O_L, DHF_SO2_G, channiwala_parikh_hhv
from propcalc.mixture import load_ingredients

OUT = Path(__file__).parent / "figures"
PC, PE = 6.895, 0.101325
DB = load_ingredients()
HHV_GRID = np.arange(16.0, 22.51, 0.5)
OX_GRID = np.arange(40.0, 86.1, 1.0)
OX_REF = 65.0
# 조성 고정: 어린 코코넛 외피 두 시료의 건조·무회분 조성 평균(CHO만). 코코피트 조성이 아님.
FORMULA = {"C": 10.0, "H": 14.45, "O": 6.35}
# 조성 민감도 확인용: 셀룰로오스형(H/C 1.67, O/C 0.83). 리그닌형(상관식 HHV 28.5)은 이 HHV 범위와 맞지 않아 뺐다.
ALT_FORMULAS = {"외피 평균 C10H14.45O6.35": FORMULA, "셀룰로오스형 C10H16.67O8.33": {"C": 10, "H": 16.67, "O": 8.33}}
# 고정 조성과 일관된 HHV 구간 = 상관식(Channiwala-Parikh) HHV ± BAND. 폭 2 MJ/kg 은 이 스크립트의 가정이다.
BAND = 2.0
# 순차(크기) 인코딩: 파랑 한 색상, 밝음→어두움 (dataviz 기준 팔레트 100→700)
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CMAP = LinearSegmentedColormap.from_list("blue_seq", SEQ)
REFS = {"셀룰로오스(면)": "CELLULOSE_COTTON", "외피 RB": "YOUNG_COCONUT_HUSK_RB", "외피 SS": "YOUNG_COCONUT_HUSK_SS"}


def hhv_daf_of(key: str) -> float:
    """DB 항목의 ΔfH° 에서 헤스 관계로 되돌린 HHV(daf, MJ/kg). dHf = A − dcH, dcH = −HHV·M."""
    i = DB[key]
    a = i.formula.get("C", 0) * DHF_CO2_G + i.formula.get("H", 0) / 2 * DHF_H2O_L + i.formula.get("S", 0) * DHF_SO2_G
    return (i.dHf_kJ_per_mol - a) / i.M_g_per_mol


def correlation_hhv(f: dict) -> float:
    """조성(원자비)에서 질량%를 구해 Channiwala-Parikh 상관식 HHV(daf, MJ/kg)를 계산한다."""
    m = molar_mass(f)
    w = {e: 100.0 * n * a / m for (e, n), a in zip(f.items(), (12.011, 1.008, 15.999))}
    return channiwala_parikh_hhv(w["C"], w["H"], w["O"])


def formula_check(hhv: float = 19.0) -> None:
    """같은 HHV 에서 조성만 바꿨을 때의 차이. 산소수지가 다르므로 고정 비율과 각자의 Isp 최대 비율 둘 다 본다."""
    print(f"[조성 민감도] HHV(daf) {hhv:g} MJ/kg 고정")
    res = {}
    for name, f in ALT_FORMULAS.items():
        d = sweep.sweep_fuel_hhv(f, [hhv], OX_GRID, Pc_MPa=PC, Pe_MPa=PE, db=DB)
        r65 = d[d.ox_wt_pct == OX_REF].iloc[0]
        ok = d[d.ok]
        b = ok.loc[ok.isp_s.idxmax()]
        res[name] = (r65, b)
        print(f"   {name:28s} ΔfH°={r65.dHf_kJ_per_mol:8.1f} kJ/mol | KNO3 {OX_REF:g}%: Tc={r65.Tc_K:7.1f} K Isp={r65.isp_s:6.2f} s"
              f" | Isp 최대: {b.isp_s:6.2f} s @ {b.ox_wt_pct:g}% (Tc={b.Tc_K:.0f} K, flag='{b.flag or ''}')")
    (a65, ab), (b65, bb) = res.values()
    print(f"   → 조성만 바꿨을 때 KNO3 {OX_REF:g}% 고정에서 Tc 폭 {abs(a65.Tc_K - b65.Tc_K):.0f} K, Isp 폭 {abs(a65.isp_s - b65.isp_s):.1f} s; "
          f"각자의 Isp 최대점에서는 Isp 폭 {abs(ab.isp_s - bb.isp_s):.1f} s, 최적 비율 차이 {abs(ab.ox_wt_pct - bb.ox_wt_pct):g}%p")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = sweep.sweep_fuel_hhv(FORMULA, HHV_GRID, OX_GRID, Pc_MPa=PC, Pe_MPa=PE, db=DB)
    df.to_csv(OUT / "hhv_map.csv", index=False)
    ok = df[df.ok]
    print(f"[가정형 HHV 지도] {df.attrs['assumptions']}")
    print(f"   계산 {len(df)}점 중 수렴 거부 {int((~df.ok).sum())}점, 모델 아티팩트 의심(flag) {int((df.flag.fillna('') != '').sum())}점 "
          "(거부·의심 점은 그림에 표시, CSV에 보존)")

    ref65 = ok[ok.ox_wt_pct == OX_REF].set_index("hhv_MJ_per_kg")
    best = ok.loc[ok.groupby("hhv_MJ_per_kg").isp_s.idxmax()].set_index("hhv_MJ_per_kg")
    print(f"\n[KNO3 {OX_REF:g}%] HHV(daf)별 성능")
    print(ref65[["dHf_kJ_per_mol", "Tc_K", "c_star_m_s", "isp_s"]].round(1).to_string())
    print("\n[HHV별 Isp 최대점]")
    print(best[["ox_wt_pct", "isp_s", "Tc_K", "cond_mass_frac"]].round(2).to_string())
    mid = ref65.loc[18.0:21.0]
    slope_T = np.polyfit(mid.index, mid.Tc_K, 1)[0]
    slope_I = np.polyfit(mid.index, mid.isp_s, 1)[0]
    print(f"\n[기울기, HHV 18~21, KNO3 {OX_REF:g}%] dTc/dHHV = {slope_T:.0f} K per MJ/kg, dIsp/dHHV = {slope_I:.1f} s per MJ/kg "
          "(README의 선형 외삽 약 150 K/(MJ/kg) 과 비교)")
    formula_check()
    hc = correlation_hhv(FORMULA)
    print(f"\n[조성-HHV 일관성] 고정 조성 {FORMULA}의 상관식 HHV = {hc:.1f} MJ/kg. 이 조성에서 일관된 HHV 구간을 "
          f"{hc - BAND:.1f}~{hc + BAND:.1f} MJ/kg(±{BAND:g}, 가정)으로 표시한다. 구간 밖은 조성과 HHV 가 서로 맞지 않는 가상 점이다.")

    refs = {k: hhv_daf_of(v) for k, v in REFS.items()}
    print("\n[참고 점, HHV(daf)] " + ", ".join(f"{k} {v:.2f}" for k, v in refs.items()) + " (코코피트 아님)")
    draw(df, ref65, best, refs)


def draw(df, ref65, best, refs) -> None:
    hc = correlation_hhv(FORMULA)
    sweep.apply_font()
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), facecolor=sweep.SURFACE, gridspec_kw={"width_ratios": [1.15, 1, 1]})
    ax = axes[0]
    ax.set_facecolor(sweep.SURFACE)
    piv = df.pivot(index="hhv_MJ_per_kg", columns="ox_wt_pct", values="isp_s")
    bad = df.assign(bad=(~df.ok) | (df.flag.fillna("") != "")).pivot(index="hhv_MJ_per_kg", columns="ox_wt_pct", values="bad")
    xs, ys = piv.columns.to_numpy(), piv.index.to_numpy()
    dx, dy = np.diff(xs).mean(), np.diff(ys).mean()
    im = ax.pcolormesh(np.append(xs - dx / 2, xs[-1] + dx / 2), np.append(ys - dy / 2, ys[-1] + dy / 2),
                       np.ma.masked_invalid(piv.to_numpy()), cmap=CMAP, shading="flat")
    for (i, h), (j, x) in ((a, b) for a in enumerate(ys) for b in enumerate(xs)):
        if bad.iloc[i, j]:
            ax.plot(x, h, marker="x" if np.isnan(piv.iloc[i, j]) else "o", markersize=5, markerfacecolor="none",
                    markeredgecolor=sweep.SERIES_COLORS[1], markeredgewidth=1.4, linestyle="none")
    ax.plot(best.ox_wt_pct, best.index, color=sweep.SERIES_COLORS[1], linewidth=2, marker="o", markersize=5,
            markeredgecolor=sweep.SURFACE, markeredgewidth=1.5, label="HHV별 Isp 최대 비율")
    ax.legend(loc="upper left", fontsize=9, frameon=True, facecolor=sweep.SURFACE, edgecolor=sweep.GRID, framealpha=1,
              labelcolor=sweep.INK)
    ax.axvline(OX_REF, color=sweep.INK_2, linewidth=1, linestyle=":")
    ax.axhline(hc - BAND, color=sweep.INK_2, linewidth=0.8, linestyle="--")
    ax.axhline(hc + BAND, color=sweep.INK_2, linewidth=0.8, linestyle="--")
    ax.set_xlabel("KNO3 질량% (나머지 가상 연료)", fontsize=9)
    ax.set_ylabel("연료 HHV, 건조·무회분 기준 (MJ/kg)", fontsize=10)
    ax.set_title("Isp (s)", loc="left", fontsize=11, color=sweep.INK)
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.outline.set_edgecolor(sweep.GRID)
    cb.ax.tick_params(colors=sweep.INK_2, labelsize=9)
    ax.tick_params(colors=sweep.INK_2, labelsize=9)
    for s in ax.spines.values():
        s.set_color(sweep.GRID)

    for ax, col, ttl, unit in ((axes[1], "isp_s", "Isp (s)", "s"), (axes[2], "Tc_K", "Tc (K)", "K")):
        sweep._style(ax)
        ax.plot(ref65.index, ref65[col], color=sweep.SERIES_COLORS[0], linewidth=2, label=f"KNO3 {OX_REF:g}% 고정")
        ax.plot(best.index, best[col], color=sweep.SERIES_COLORS[1], linewidth=2, label="HHV별 Isp 최대 비율")
        for ser in (ref65, best):
            fl = ser[ser.flag.fillna("") != ""]
            ax.plot(fl.index, fl[col], linestyle="none", marker="o", markersize=8, markerfacecolor=sweep.SURFACE,
                    markeredgecolor=sweep.INK_2, markeredgewidth=1.6, zorder=4)
        ax.axvspan(hc - BAND, hc + BAND, color=sweep.GRID, alpha=0.55, linewidth=0, zorder=0)
        top = ax.get_ylim()[1]
        for name, x in refs.items():
            ax.axvline(x, color=sweep.INK_2, linewidth=1, linestyle=":")
            ax.text(x, top, name, rotation=90, va="top", ha="right", fontsize=8.5, color=sweep.INK_2)
        ax.set_xlabel("연료 HHV, 건조·무회분 기준 (MJ/kg)", fontsize=9)
        ax.set_ylabel(ttl, fontsize=10)
        ax.set_title(ttl + " vs 연료 HHV", loc="left", fontsize=11, color=sweep.INK)
        ax.legend(frameon=False, fontsize=9, labelcolor=sweep.INK, loc="lower right")
    fig.suptitle("가정형 HHV 지도 (코코피트 결과 아님): 연료 HHV 가 성능을 움직이는 크기", x=0.04, ha="left", fontsize=13, color=sweep.INK)
    fig.text(0.04, 0.012, "가정: " + df.attrs["assumptions"] + "\n점선 = 참고 시료의 HHV(daf), 코코피트 아님. 회색 음영·파선 = 고정 조성과 일관된 HHV 구간(상관식 ±2, 가정). 히트맵의 주황 ○/× = 상전이 plateau 의심 점 / CEA 비물리적 결과로 제외된 점, 선 그래프의 회색 ○ = 의심 점."
             " 표: hhv_map.csv", fontsize=8.5, color=sweep.INK_2, ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.07, 1, 0.94))
    fig.savefig(OUT / "hhv_map.png", dpi=150, facecolor=sweep.SURFACE)


if __name__ == "__main__":
    main()
