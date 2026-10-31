"""커피박(SCG) 성능 계산: 한 논문의 한 시료 값과 문헌 간 산포.

    python notebooks/coffee_grounds_analysis.py

입력은 ACS Omega 2021(PMC8756566)의 단일 시료이고, ΔfH° 불확실도(±436 kJ/mol)는 측정 SD 가 아니라 16개 SCG HHV 값의
문헌 간 산포(SD 1.80 MJ/kg)다. 따라서 아래 범위는 "커피박 일반이 이 범위 안에 있다"는 보증이 아니라, 문헌에 보고된
HHV 산포가 성능에 주는 규모다. 회분(0.9%)은 반응물에서 제외된다.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from propcalc import sweep
from propcalc.mixture import Mixture, load_ingredients

OUT = Path(__file__).parent / "figures"
PC, PE = 6.895, 0.101325
DB = load_ingredients()
KEY = "SPENT_COFFEE_GROUNDS"
FUELS = {"소르비톨(참고)": "SORBITOL", "셀룰로오스(참고)": "CELLULOSE_COTTON", "커피박": KEY}
ASSUMPTION = (f"평형 유동, Pc={PC} MPa, Pe={PE} MPa, cea 백엔드, KNO3 65 wt% + 연료 35 wt%, 회분 제외, "
              "커피박 ΔfH° 불확실도는 문헌 간 HHV 산포(SD 1.80 MJ/kg)")


def bounds() -> pd.DataFrame:
    rows = []
    m = Mixture.from_wt_pct({"KNO3": 65.0, KEY: 35.0}, DB)
    ing = DB[KEY]
    for u_label, u in (("문헌 SD 1개 폭(±436)", ing.u_dHf_kJ_per_mol), ("문헌 SD 0.5개 폭(±218)", ing.u_dHf_kJ_per_mol / 2)):
        s = sweep.sensitivity_dHf(m, KEY, u_kJ_per_mol=u, Pc_MPa=PC, Pe_MPa=PE)
        for point, sign in (("low", -1), ("nominal", 0), ("high", 1)):
            t = s.table.loc[point]
            rows.append({"범위": u_label, "point": point, "dHf_kJ_per_mol": round(ing.dHf_kJ_per_mol + sign * u, 1),
                         "Tc_K": t.Tc_K, "c_star_m_s": t.c_star_m_s, "isp_s": t.isp_s})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    print(f"[커피박 ΔfH° 한계값, KNO3 65%] {ASSUMPTION}")
    try:
        b = bounds()
        b.to_csv(OUT / "coffee_grounds_bounds.csv", index=False)
        print(b.round(1).to_string(index=False))
    except Exception as e:  # 한계값이 CEA 수렴 영역을 벗어나면 숨기지 않고 알린다
        print("한계값 계산 실패:", type(e).__name__, str(e)[:200])
    df = sweep.compare_fuels(FUELS, np.arange(40, 86, 1.0), Pc_MPa=PC, Pe_MPa=PE, db=DB)
    df.to_csv(OUT / "coffee_grounds_sweep.csv", index=False)
    best = sweep.best_by_isp(df)
    best["거부 행"] = [int((~df[df.fuel == f].ok).sum()) for f in best.fuel]
    best["의심 행"] = [int(((df.fuel == f) & df.ok & (df.flag.fillna("") != "")).sum()) for f in best.fuel]
    print("\n[산화제 스윕 최대 Isp]")
    print(best[["fuel", "ox_wt_pct", "isp_s", "Tc_K", "cond_mass_frac", "거부 행", "의심 행"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
