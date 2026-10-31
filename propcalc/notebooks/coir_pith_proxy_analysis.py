"""코코피트 대용(proxy) 분석: 어린 코코넛 외피(RMUTP 2025)의 ΔfH° 범위가 성능에 미치는 영향.

    python notebooks/coir_pith_proxy_analysis.py

이 분석은 코코피트 자체의 계산이 아니다. DB 에는 코코피트 항목이 없고(같은 시료의 원소분석·회분·수분·실측 HHV 를
원문에서 확인한 코코피트 문헌이 없음), 같은 논문의 어린 코코넛 외피 두 시료만 있다. 코코피트가 이 외피와 비슷하다는
것은 검증되지 않은 가정이며, 아래 범위는 코코피트의 실제 범위가 아니라 "이 두 외피 시료가 보여 주는 규모"다.
회분(건조 환산 5~6%)은 반응물에서 제외되어 있다.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from propcalc import sweep
from propcalc.mixture import Mixture, load_ingredients

OUT = Path(__file__).parent / "figures"
PC, PE = 6.895, 0.101325
OX_PCT = 65.0
DB = load_ingredients()
HUSKS = {"외피 RB": "YOUNG_COCONUT_HUSK_RB", "외피 SS": "YOUNG_COCONUT_HUSK_SS"}
ASSUMPTION = (f"평형 유동, Pc={PC} MPa, Pe={PE} MPa, cea 백엔드, KNO3 {OX_PCT:g} wt% + 연료 {100 - OX_PCT:g} wt%, "
              "연료는 어린 코코넛 외피(코코피트 아님), 회분 제외, ΔfH°는 수분 포함 기준 가정의 환산값")


def bounds_table() -> pd.DataFrame:
    rows = []
    for label, key in HUSKS.items():
        m = Mixture.from_wt_pct({"KNO3": OX_PCT, key: 100 - OX_PCT}, DB)
        s = sweep.sensitivity_dHf(m, key, Pc_MPa=PC, Pe_MPa=PE)
        for point in ("low", "nominal", "high"):
            t = s.table.loc[point]
            dHf = s.dHf_kJ_per_mol + {"low": -s.u_kJ_per_mol, "nominal": 0.0, "high": s.u_kJ_per_mol}[point]
            rows.append({"fuel": label, "point": point, "dHf_kJ_per_mol": round(dHf, 1),
                         "Tc_K": t.Tc_K, "c_star_m_s": t.c_star_m_s, "isp_s": t.isp_s})
    return pd.DataFrame(rows)


def optimum_table() -> pd.DataFrame:
    fuels = {"소르비톨(참고)": "SORBITOL", **HUSKS}
    df = sweep.compare_fuels(fuels, np.arange(40, 86, 1.0), Pc_MPa=PC, Pe_MPa=PE, db=DB)
    df.to_csv(OUT / "coir_pith_proxy_sweep.csv", index=False)
    best = sweep.best_by_isp(df)
    best["fail_rows"] = [int((~df[df.fuel == f].ok).sum()) for f in best.fuel]
    return best[["fuel", "ox_wt_pct", "isp_s", "Tc_K", "cond_mass_frac", "fail_rows"]]


def main() -> None:
    OUT.mkdir(exist_ok=True)
    b = bounds_table()
    b.to_csv(OUT / "coir_pith_proxy_bounds.csv", index=False)
    print(f"[ΔfH° 한계값, KNO3 {OX_PCT:g}%] {ASSUMPTION}")
    print(b.round(1).to_string(index=False))
    lo, hi = b.Tc_K.min(), b.Tc_K.max()
    print(f"\n두 외피 시료의 한계값 전체 범위: Tc {lo:.0f}~{hi:.0f} K (폭 {hi - lo:.0f} K), "
          f"Isp {b.isp_s.min():.1f}~{b.isp_s.max():.1f} s (폭 {b.isp_s.max() - b.isp_s.min():.1f} s)")
    print("\n[산화제 스윕 최대 Isp]")
    print(optimum_table().round(3).to_string(index=False))


if __name__ == "__main__":
    main()
