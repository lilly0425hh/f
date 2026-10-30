"""문헌 참고값 비교 (범위/대략 비교 전용, 실패 판정에 쓰지 않는다).

Nakka(아마추어 자료, 신뢰도 낮음) 와 PROPEP 은 칼륨 화합물 종 DB 가 CEA 와 다르다.
"""

from __future__ import annotations

from .backends.base import Result
from .mixture import Mixture

# Nakka/PROPEP KNSB 평형식, 6.8 MPa: C6H14O6 + 3.345 KNO3 -> ... (소르비톨 1 mol 당 mol)
# 출처: https://www.nakka-rocketry.net/sorbchem.html (아마추어 자료)
NAKKA_KNSB_6_8MPA: dict[str, float] = {
    "CO2": 1.870, "CO": 2.490, "H2O": 4.828, "H2": 2.145, "N2": 1.672, "K2CO3": 1.644, "KOH": 0.057,
}


def species_per_mol_fuel(result: Result, mixture: Mixture, fuel_key: str) -> dict[str, float]:
    """연소실 종별 몰수를 연료 1 mol 기준으로 환산한다. K2CO3 는 상(L/고체) 구분 없이 합친다."""
    n_fuel = mixture.moles_per_kg(fuel_key)
    out: dict[str, float] = {}
    for name, mol_per_kg in result.species_moles_per_kg.items():
        base = name.split("(")[0] if name.endswith(")") else name
        out[base] = out.get(base, 0.0) + mol_per_kg / n_fuel
    return out


def compare_to_nakka(result: Result, mixture: Mixture, fuel_key: str = "SORBITOL") -> str:
    """CEA 결과와 Nakka 평형식을 표로 비교한 문자열 (KOH 는 KOH 와 K2O2H2 를 분리해 표시)."""
    ours = species_per_mol_fuel(result, mixture, fuel_key)
    rows = [f"{'종':8s} {'Nakka':>8s} {'이 계산':>8s}"]
    for sp, ref in NAKKA_KNSB_6_8MPA.items():
        rows.append(f"{sp:8s} {ref:8.3f} {ours.get(sp, 0.0):8.3f}")
    others = sorted(((k, v) for k, v in ours.items() if k not in NAKKA_KNSB_6_8MPA and v > 0.01),
                    key=lambda kv: -kv[1])
    for k, v in others:
        rows.append(f"{k:8s} {'-':>8s} {v:8.3f}")
    return "\n".join(rows)
