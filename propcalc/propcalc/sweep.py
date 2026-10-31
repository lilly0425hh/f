"""배합비 스윕, 연료 대체 비교, ΔfH° 민감도 분석, 그래프.

모든 표(DataFrame)는 attrs["assumptions"] 에 계산 가정(평형/동결, Pc, Pe, 백엔드)을 가진다.
수렴하지 못한 점은 ok=False 와 error 문구로 남기고 값은 NaN 으로 둔다 (조용히 건너뛰지 않는다).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .backends import ConvergenceError, get_backend
from .backends.base import Result
from .hess import dHf_from_hhv
from .mixture import Ingredient, Mixture, load_ingredients, stoichiometry

DEFAULT_PC_MPA = 6.895
DEFAULT_PE_MPA = 0.101325
_METRICS = ("Tc_K", "c_star_m_s", "isp_s", "isp_vac_s", "M_eff_g_mol", "cond_mass_frac")


def _base_name(species: str) -> str:
    return species.split("(")[0]


def _flag(r: Result) -> str:
    """평형 모델이 상변화 구간에서 만드는 아티팩트가 의심되는 점을 표시한다."""
    flags = []
    if r.gamma_s < 1.0:
        flags.append("γ_s<1")
    names = [_base_name(k) for k in r.condensed_mass_fractions]
    if len(names) != len(set(names)):
        flags.append("연소실 상 공존")
    return ", ".join(flags)


def _row(r: Result | None, error: str = "") -> dict:
    if r is None:
        return {"ok": False, "error": error, "flag": "", **{m: np.nan for m in _METRICS},
                "Tt_K": np.nan, "gamma_s": np.nan, "C_gr_mass_frac": np.nan, "has_C_gr": False}
    cgr = r.condensed_mass_fractions.get("C(gr)", 0.0)
    return {
        "ok": True, "error": "", "Tc_K": r.Tc_K, "Tt_K": r.Tt_K, "c_star_m_s": r.c_star_m_s,
        "isp_s": r.isp_s, "isp_vac_s": r.isp_vac_s, "M_eff_g_mol": r.M_eff_g_mol, "gamma_s": r.gamma_s,
        "cond_mass_frac": r.condensed_mass_fraction, "C_gr_mass_frac": cgr, "has_C_gr": cgr > 1e-6,
        "flag": _flag(r),
    }


def _resolve(item: str | Ingredient, db: Mapping[str, Ingredient]) -> Ingredient:
    return item if isinstance(item, Ingredient) else db[item]


def _assumption_text(backend: str, mode: str, Pc: float, Pe: float, note: str = "") -> str:
    flow = {"equilibrium": "평형 유동", "frozen": "동결 유동"}[mode]
    return f"{flow}, Pc={Pc:g} MPa, Pe={Pe:.4g} MPa, 백엔드={backend}" + (f" | {note}" if note else "")


def _run_safe(be, mixture: Mixture, Pc: float, Pe: float, mode: str) -> tuple[Result | None, str]:
    try:
        return be.run(mixture, Pc, Pe_MPa=Pe, mode=mode), ""
    except ConvergenceError as e:
        return None, str(e)[:200]


def sweep_oxidizer(
    fuel: str | Ingredient,
    oxidizer: str | Ingredient = "KNO3",
    ox_wt_pct: Iterable[float] = range(50, 86),
    Pc_MPa: float = DEFAULT_PC_MPA,
    Pe_MPa: float = DEFAULT_PE_MPA,
    backend: str = "cea",
    mode: str = "equilibrium",
    db: Mapping[str, Ingredient] | None = None,
) -> pd.DataFrame:
    """산화제 질량% 를 바꿔 가며 Tc, c*, Isp, 응축상, C(gr) 를 계산한다."""
    db = db if db is not None else load_ingredients()
    ox, fu = _resolve(oxidizer, db), _resolve(fuel, db)
    be = get_backend(backend)
    rows = []
    for x in ox_wt_pct:
        m = Mixture.from_wt_pct({ox.key: float(x), fu.key: 100.0 - float(x)}, db, extra=[ox, fu])
        r, err = _run_safe(be, m, Pc_MPa, Pe_MPa, mode)
        rows.append({"fuel": fu.key, "ox_wt_pct": float(x), "fuel_wt_pct": 100.0 - float(x), **_row(r, err)})
    df = pd.DataFrame(rows)
    # 거부된 점 바로 옆의 점은 값이 매끄러워 보여도 불안정 구간 안이므로 표시한다
    near_bad = (~df["ok"]).shift(1, fill_value=False) | (~df["ok"]).shift(-1, fill_value=False)
    mark = df["ok"] & near_bad
    df.loc[mark, "flag"] = df.loc[mark, "flag"].map(lambda f: (f + ", " if f else "") + "인접 점 거부")
    _flag_plateau_neighborhood(df)
    try:
        st = stoichiometry(ox, fu).mass_pct_oxidizer
    except (NotImplementedError, ValueError):
        st = float("nan")
    df["stoich_ox_wt_pct"] = st
    df.attrs["assumptions"] = _assumption_text(backend, mode, Pc_MPa, Pe_MPa, f"{ox.key}/{fu.key}")
    df.attrs["stoich_ox_wt_pct"] = st
    df.attrs["oxidizer"], df.attrs["fuel"] = ox.key, fu.key
    return df


def _flag_plateau_neighborhood(df: pd.DataFrame, window_K: float = 20.0) -> None:
    """상 공존이 관측된 온도(상전이 plateau)를 데이터에서 찾아, 연소실 또는 목 온도가 그 ±window_K 안인
    점을 표시한다. 평형 모델은 plateau 에서 γ_s, c* 가 불연속적으로 변해 곡선이 꺾인다."""
    co = df["ok"] & df["flag"].str.contains("상 공존")
    for T in sorted(df.loc[co, "Tc_K"].round(0).unique()):
        near = df["ok"] & (((df["Tc_K"] - T).abs() <= window_K) | ((df["Tt_K"] - T).abs() <= window_K))
        tag = f"상전이 plateau(≈{T:.0f} K) 근처"
        df.loc[near, "flag"] = df.loc[near, "flag"].map(lambda f: f if tag in f else (f + ", " if f else "") + tag)


def isp_optimum(df: pd.DataFrame) -> dict:
    """Isp 최대점과 화학량론 대비 위치."""
    ok = df[df["ok"]]
    i = ok["isp_s"].idxmax()
    st = float(df["stoich_ox_wt_pct"].iloc[0])
    return {"ox_wt_pct": float(ok.loc[i, "ox_wt_pct"]), "isp_s": float(ok.loc[i, "isp_s"]),
            "Tc_K": float(ok.loc[i, "Tc_K"]), "M_eff_g_mol": float(ok.loc[i, "M_eff_g_mol"]),
            "cond_mass_frac": float(ok.loc[i, "cond_mass_frac"]), "stoich_ox_wt_pct": st,
            "fuel_rich": bool(ok.loc[i, "ox_wt_pct"] < st)}


def compare_fuels(
    fuels: Mapping[str, str | Ingredient],
    ox_wt_pct: Sequence[float],
    oxidizer: str | Ingredient = "KNO3",
    Pc_MPa: float = DEFAULT_PC_MPA,
    Pe_MPa: float = DEFAULT_PE_MPA,
    backend: str = "cea",
    mode: str = "equilibrium",
    db: Mapping[str, Ingredient] | None = None,
) -> pd.DataFrame:
    """연료별로 같은 산화제 비율에서 성능을 비교한 긴 형식 표 (fuel 열 = 라벨)."""
    db = db if db is not None else load_ingredients()
    frames = []
    for label, fu in fuels.items():
        d = sweep_oxidizer(fu, oxidizer, ox_wt_pct, Pc_MPa, Pe_MPa, backend, mode, db)
        d["fuel"] = label
        frames.append(d)
    out = pd.concat(frames, ignore_index=True)
    out.attrs["assumptions"] = frames[0].attrs["assumptions"].rsplit(" | ", 1)[0]
    return out


def sweep_fuel_hhv(
    formula: str | Mapping[str, float],
    hhv_MJ_per_kg: Iterable[float],
    ox_wt_pct: Iterable[float] = range(50, 81, 2),
    oxidizer: str | Ingredient = "KNO3",
    Pc_MPa: float = DEFAULT_PC_MPA,
    Pe_MPa: float = DEFAULT_PE_MPA,
    backend: str = "cea",
    mode: str = "equilibrium",
    db: Mapping[str, Ingredient] | None = None,
) -> pd.DataFrame:
    """가정형 지도: 조성(formula)을 고정한 가상 연료의 HHV 만 바꿔 가며 산화제 스윕을 반복한다.

    HHV 는 건조·무회분(daf) 기준이며, 측정값이 아니라 격자 값이다. 가상 연료는 verification=doc_only 이고
    어떤 실제 시료도 대표하지 않는다. 실제 시료의 값이 확보되면 이 지도 위에 점으로 얹어 해석하는 용도다.
    긴 형식 표를 반환한다 (hhv_MJ_per_kg, dHf_kJ_per_mol 열 추가).
    """
    db = db if db is not None else load_ingredients()
    frames = []
    for hhv in hhv_MJ_per_kg:
        r = dHf_from_hhv(formula, float(hhv), None, hhv_source="measured")
        fuel = Ingredient.custom(
            f"HYPO{int(round(hhv * 10))}", f"가상 연료 HHV={hhv:g} MJ/kg(daf)", r.formula, r.dHf_kJ_per_mol, None,
            source=f"가정형 HHV 지도 격자 값 (HHV={hhv:g} MJ/kg, 조성 고정). 측정값이 아님",
            verification="doc_only", note="어떤 실제 시료도 대표하지 않는다")
        d = sweep_oxidizer(fuel, oxidizer, ox_wt_pct, Pc_MPa, Pe_MPa, backend, mode, db)
        d.insert(0, "dHf_kJ_per_mol", r.dHf_kJ_per_mol)
        d.insert(0, "hhv_MJ_per_kg", float(hhv))
        frames.append(d)
    out = pd.concat(frames, ignore_index=True)
    f_txt = " ".join(f"{el}{v:g}" for el, v in r.formula.items())
    out.attrs["assumptions"] = frames[0].attrs["assumptions"].rsplit(" | ", 1)[0] + \
        f" | 가정형 지도: 연료 조성 {f_txt}(고정), HHV(daf)는 격자 값, 회분 제외"
    out.attrs["oxidizer"] = frames[0].attrs["oxidizer"]
    return out


def best_by_isp(df: pd.DataFrame) -> pd.DataFrame:
    """연료별 Isp 최대 행."""
    ok = df[df["ok"]]
    return ok.loc[ok.groupby("fuel")["isp_s"].idxmax()].reset_index(drop=True)


def sweep_fuel_blend(
    fuel_a: str | Ingredient,
    fuel_b: str | Ingredient,
    fraction_b_pct: Iterable[float] = range(0, 101, 10),
    ox_wt_pct: float = 65.0,
    oxidizer: str | Ingredient = "KNO3",
    Pc_MPa: float = DEFAULT_PC_MPA,
    Pe_MPa: float = DEFAULT_PE_MPA,
    backend: str = "cea",
    mode: str = "equilibrium",
    db: Mapping[str, Ingredient] | None = None,
) -> pd.DataFrame:
    """산화제 비율을 고정하고 연료 부분 중 fuel_b 의 질량% 를 바꾼다 (혼합 연료)."""
    db = db if db is not None else load_ingredients()
    ox, fa, fb = _resolve(oxidizer, db), _resolve(fuel_a, db), _resolve(fuel_b, db)
    be = get_backend(backend)
    fuel_part = 100.0 - ox_wt_pct
    rows = []
    for xb in fraction_b_pct:
        wt = {ox.key: float(ox_wt_pct), fa.key: fuel_part * (1 - xb / 100.0), fb.key: fuel_part * xb / 100.0}
        wt = {k: v for k, v in wt.items() if v > 0}
        m = Mixture.from_wt_pct(wt, db, extra=[ox, fa, fb])
        r, err = _run_safe(be, m, Pc_MPa, Pe_MPa, mode)
        rows.append({"fuel_b_pct_of_fuel": float(xb), "ox_wt_pct": ox_wt_pct, **_row(r, err)})
    df = pd.DataFrame(rows)
    df.attrs["assumptions"] = _assumption_text(backend, mode, Pc_MPa, Pe_MPa,
                                               f"{ox.key} {ox_wt_pct:g}% + 연료({fa.key}/{fb.key})")
    return df


@dataclass(frozen=True)
class SensitivityResult:
    method: str
    key: str
    dHf_kJ_per_mol: float
    u_kJ_per_mol: float
    table: pd.DataFrame          # 행: bounds -> low/nominal/high, mc -> mean/std/p2.5/p97.5/nominal
    samples: pd.DataFrame | None
    n_failed: int
    assumptions: str


def sensitivity_dHf(
    mixture: Mixture,
    key: str,
    u_kJ_per_mol: float | None = None,
    method: str = "bounds",
    n: int = 200,
    seed: int = 0,
    Pc_MPa: float = DEFAULT_PC_MPA,
    Pe_MPa: float = DEFAULT_PE_MPA,
    backend: str = "cea",
) -> SensitivityResult:
    """한 성분의 ΔfH° 불확실도가 Tc, c*, Isp 에 미치는 영향.

    bounds: ΔfH° ± u 를 한계값으로 보고 3점 계산.
    mc: ΔfH° ~ Normal(평균, σ=u) 로 n 번 추출 (u 를 1σ 로 해석; 출처의 ± 가 1σ 인지 한계인지는
        출처마다 다르므로 해석 선택을 결과 해석에 명시할 것).
    u 가 None 이면 DB 의 u_dHf 를 쓰며, 그것도 없으면(미확인) ValueError.
    """
    ing = mixture.ingredient(key)
    u = u_kJ_per_mol if u_kJ_per_mol is not None else ing.u_dHf_kJ_per_mol
    if u is None:
        raise ValueError(f"{key}: ΔfH° 불확실도가 미확인이다. u_kJ_per_mol 을 명시하라")
    if u < 0:
        raise ValueError("u 는 0 이상이어야 한다")
    be = get_backend(backend)
    nominal = ing.dHf_kJ_per_mol

    def metrics(dHf: float) -> dict | None:
        try:
            r = be.run(mixture.replace_dHf(key, dHf), Pc_MPa, Pe_MPa=Pe_MPa)
        except ConvergenceError:
            return None
        return {"Tc_K": r.Tc_K, "c_star_m_s": r.c_star_m_s, "isp_s": r.isp_s, "isp_vac_s": r.isp_vac_s}

    assumptions = _assumption_text(backend, "equilibrium", Pc_MPa, Pe_MPa,
                                   f"{key} ΔfH°={nominal:g}±{u:g} kJ/mol ({method})")
    if method == "bounds":
        pts = {"low": metrics(nominal - u), "nominal": metrics(nominal), "high": metrics(nominal + u)}
        if any(v is None for v in pts.values()):
            raise ConvergenceError("민감도 경계점 계산이 수렴하지 못했다")
        return SensitivityResult("bounds", key, nominal, u, pd.DataFrame(pts).T, None, 0, assumptions)
    if method == "mc":
        rng = np.random.default_rng(seed)
        draws = rng.normal(nominal, u, size=n)
        got = [(d, metrics(d)) for d in draws]
        good = [(d, m) for d, m in got if m is not None]
        samples = pd.DataFrame([{"dHf_kJ_per_mol": d, **m} for d, m in good])
        cols = ["Tc_K", "c_star_m_s", "isp_s", "isp_vac_s"]
        table = pd.DataFrame({
            "mean": samples[cols].mean(), "std": samples[cols].std(ddof=1),
            "p2.5": samples[cols].quantile(0.025), "p97.5": samples[cols].quantile(0.975),
        }).T
        table.loc["nominal"] = pd.Series(metrics(nominal))
        return SensitivityResult("mc", key, nominal, u, table, samples, len(got) - len(good), assumptions)
    raise ValueError("method 는 'bounds' 또는 'mc'")


# --- 그래프 ---------------------------------------------------------------------

# dataviz 기준 팔레트 (밝은 표면). 범주 슬롯 순서 고정: 파랑, 주황, 아쿠아, 노랑.
# validate_palette.js --mode light 에서 4색 모두 PASS; 아쿠아·노랑은 대비 3:1 미만이라
# 직접 라벨과 표(CSV) 제공이 필수 (relief rule).
SERIES_COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e3e2de", "#fcfcfb"
_FONT_CANDIDATES = ["NanumGothic", "Noto Sans CJK KR", "Malgun Gothic", "AppleGothic", "WenQuanYi Zen Hei"]


def apply_font() -> str:
    """한글을 그릴 수 있는 설치된 폰트를 하나 골라 matplotlib 에 적용한다 (없으면 DejaVu Sans)."""
    import logging

    from matplotlib import font_manager, pyplot as plt

    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    installed = {f.name for f in font_manager.fontManager.ttflist}
    name = next((c for c in _FONT_CANDIDATES if c in installed), "DejaVu Sans")
    plt.rcParams["font.family"] = name
    plt.rcParams["axes.unicode_minus"] = False
    return name


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=INK_2, labelsize=9, length=0)
    ax.xaxis.label.set_color(INK_2)
    ax.yaxis.label.set_color(INK_2)


def plot_sweep(df: pd.DataFrame, path: str | Path | None = None, title: str | None = None):
    """산화제 스윕 4분할 그래프 (Tc, c*, Isp, 응축상/C(gr)). 가정 조건을 그림 하단에 표기한다.

    축은 패널마다 하나이며 이중축은 쓰지 않는다. 화학량론 위치는 점선, Isp 최대점은 강조 마커,
    평형 모델 아티팩트 의심 점(flag)은 속이 빈 주황 마커로 표시한다.
    """
    import matplotlib.pyplot as plt

    apply_font()
    ok = df[df["ok"]]
    x = ok["ox_wt_pct"].to_numpy()
    st = df.attrs.get("stoich_ox_wt_pct", float("nan"))
    ox_name, fu_name = df.attrs.get("oxidizer", "산화제"), df.attrs.get("fuel", "연료")
    flagged = ok[ok["flag"] != ""]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2), facecolor=SURFACE)
    panels = [("Tc_K", "연소실 온도 Tc (K)"), ("c_star_m_s", "특성속도 c* (m/s)"),
              ("isp_s", "비추력 Isp (s, Pa=Pe)"), (None, "응축상 질량분율")]
    for ax, (col, ylabel) in zip(axes.ravel(), panels):
        _style(ax)
        if col:
            ax.plot(x, ok[col], color=SERIES_COLORS[0], linewidth=2, solid_capstyle="round")
            if len(flagged):
                ax.plot(flagged["ox_wt_pct"], flagged[col], linestyle="none", marker="o", markersize=8,
                        markerfacecolor=SURFACE, markeredgecolor=SERIES_COLORS[1], markeredgewidth=2)
            if col == "isp_s":
                i = ok[col].idxmax()
                ax.plot(ok.loc[i, "ox_wt_pct"], ok.loc[i, col], marker="o", markersize=9,
                        color=SERIES_COLORS[0], markeredgecolor=SURFACE, markeredgewidth=2)
                ax.annotate(f"최대 {ok.loc[i, col]:.1f} s\n@ {ok.loc[i, 'ox_wt_pct']:g}%",
                            (ok.loc[i, "ox_wt_pct"], ok.loc[i, col]), xytext=(-70, -38),
                            textcoords="offset points", fontsize=9, color=INK)
        else:
            ax.plot(x, ok["cond_mass_frac"], color=SERIES_COLORS[0], linewidth=2)
            ax.plot(x, ok["C_gr_mass_frac"], color=SERIES_COLORS[1], linewidth=2, linestyle=(0, (5, 2)))
            ax.text(x[-1], ok["cond_mass_frac"].iloc[-1], " 전체 응축상", color=INK, fontsize=9, va="center")
            ax.text(x[-1], max(ok["C_gr_mass_frac"].iloc[-1], 0.0), " C(gr)", color=INK, fontsize=9, va="center")
        if st == st:
            ax.axvline(st, color=INK_2, linewidth=1, linestyle=":")
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_xlabel(f"{ox_name} 질량% (나머지 {fu_name})", fontsize=9)
    if st == st:
        axes[0, 0].annotate(f"화학량론 {st:.1f}%", (st, axes[0, 0].get_ylim()[0]), xytext=(-4, 6),
                            textcoords="offset points", ha="right", fontsize=9, color=INK_2)
    if len(flagged):
        axes[0, 1].text(0.02, 0.04, "○ 상변화 구간 의심 (γ_s<1 / 상 공존)", transform=axes[0, 1].transAxes,
                        fontsize=9, color=INK_2)
    fig.suptitle(title or f"{ox_name}/{fu_name} 산화제 비율 스윕", x=0.06, ha="left", fontsize=13, color=INK)
    fig.text(0.06, 0.012, "가정: " + df.attrs.get("assumptions", "(미기재)"), fontsize=8.5, color=INK_2, ha="left")
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    if path is not None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150, facecolor=SURFACE)
        df.to_csv(path.with_suffix(".csv"), index=False)   # 표 보기 (색 대비 보완)
    return fig
