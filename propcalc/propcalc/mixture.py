"""배합 전처리: 질량% -> 원소 몰수, 혼합물 엔탈피, 산소수지, 화학량론비.

수치 데이터(생성엔탈피 등)는 data/ingredients.yaml 에서만 읽으며,
각 항목은 출처(source)와 검증 수준(verification)을 반드시 가져야 한다.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml

from .formula import ATOMIC_WEIGHTS, as_formula_dict, molar_mass

DEFAULT_DB_PATH = Path(__file__).resolve().parents[1] / "data" / "ingredients.yaml"
VERIFICATION_LEVELS = ("primary", "secondary", "doc_only")
T_REF_K = 298.15
_SUPPORTED_STOICH_ELEMENTS = {"C", "H", "O", "N", "K"}


@dataclass(frozen=True)
class Ingredient:
    key: str
    name: str
    formula: dict[str, float]
    dHf_kJ_per_mol: float
    u_dHf_kJ_per_mol: float | None
    M_g_per_mol: float
    phase: str
    source: str
    verification: str
    note: str = ""
    density_g_cm3: float | None = None

    def __post_init__(self) -> None:
        if not self.source or not str(self.source).strip():
            raise ValueError(f"{self.key}: source 가 비어 있다 (출처 없는 값은 넣지 않는다)")
        if self.verification not in VERIFICATION_LEVELS:
            raise ValueError(f"{self.key}: verification 은 {VERIFICATION_LEVELS} 중 하나여야 한다")

    @classmethod
    def custom(
        cls,
        key: str,
        name: str,
        formula: str | Mapping[str, float],
        dHf_kJ_per_mol: float,
        u_dHf_kJ_per_mol: float | None,
        source: str,
        verification: str = "doc_only",
        phase: str = "cr",
        note: str = "",
    ) -> Ingredient:
        f = as_formula_dict(formula)
        return cls(key, name, f, float(dHf_kJ_per_mol), u_dHf_kJ_per_mol, molar_mass(f),
                   phase, source, verification, note)

    def with_dHf(self, dHf_kJ_per_mol: float) -> Ingredient:
        return dataclasses.replace(self, dHf_kJ_per_mol=float(dHf_kJ_per_mol))


def load_ingredients(path: str | Path | None = None) -> dict[str, Ingredient]:
    """ingredients.yaml 을 읽는다. source/verification 이 없으면 ValueError."""
    p = Path(path) if path else DEFAULT_DB_PATH
    raw = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    out: dict[str, Ingredient] = {}
    for key, e in raw.items():
        for field_name in ("name", "formula", "dHf_kJ_per_mol", "M_g_per_mol", "phase", "source", "verification"):
            if field_name not in e or e[field_name] in (None, ""):
                raise ValueError(f"{key}: 필수 필드 '{field_name}' 가 없다 (source 포함)")
        out[key] = Ingredient(
            key=key,
            name=e["name"],
            formula=as_formula_dict(e["formula"]),
            dHf_kJ_per_mol=float(e["dHf_kJ_per_mol"]),
            u_dHf_kJ_per_mol=None if e.get("u_dHf_kJ_per_mol") is None else float(e["u_dHf_kJ_per_mol"]),
            M_g_per_mol=float(e["M_g_per_mol"]),
            phase=str(e["phase"]),
            source=str(e["source"]).strip(),
            verification=str(e["verification"]),
            note=str(e.get("note", "") or "").strip(),
            density_g_cm3=e.get("density_g_cm3"),
        )
    return out


@dataclass(frozen=True)
class Stoich:
    mol_ox_per_mol_fuel: float
    mass_pct_oxidizer: float
    mass_pct_fuel: float
    products: dict[str, float]


def stoichiometry(oxidizer: Ingredient, fuel: Ingredient) -> Stoich:
    """완전산화 화학량론. 생성물 가정: K 가 있으면 K2CO3, CO2, H2O, N2, 없으면 CO2, H2O, N2.

    원소수지 연립방정식을 직접 풀어 산화제 몰수를 구한다.
    """
    elems = set(oxidizer.formula) | set(fuel.formula)
    if not elems <= _SUPPORTED_STOICH_ELEMENTS:
        raise NotImplementedError(f"지원하지 않는 원소: {sorted(elems - _SUPPORTED_STOICH_ELEMENTS)}")

    product_formulas: dict[str, dict[str, float]] = {}
    if "K" in oxidizer.formula:
        product_formulas["K2CO3"] = {"K": 2, "C": 1, "O": 3}
    product_formulas.update({"CO2": {"C": 1, "O": 2}, "H2O": {"H": 2, "O": 1}, "N2": {"N": 2}})

    eq_elems = [e for e in ("C", "H", "O", "N", "K") if e in elems or e in
                {el for f in product_formulas.values() for el in f}]
    names = list(product_formulas)
    # 미지수: [x(산화제 몰수), 생성물 몰수...]; 식: x*ox[e] - sum(p*comp[e]) = -fuel[e]
    A = np.zeros((len(eq_elems), 1 + len(names)))
    rhs = np.zeros(len(eq_elems))
    for i, e in enumerate(eq_elems):
        A[i, 0] = oxidizer.formula.get(e, 0.0)
        for j, p in enumerate(names):
            A[i, 1 + j] = -product_formulas[p].get(e, 0.0)
        rhs[i] = -fuel.formula.get(e, 0.0)
    if A.shape[0] != A.shape[1]:
        raise NotImplementedError("이 산화제/연료 조합은 정방 원소수지로 풀리지 않는다")
    sol = np.linalg.solve(A, rhs)
    x, prods = float(sol[0]), {n: float(v) for n, v in zip(names, sol[1:])}
    if x <= 0 or any(v < -1e-12 for v in prods.values()):
        raise ValueError(f"화학량론 해가 물리적이지 않다 (x={x}, products={prods})")
    m_ox = x * oxidizer.M_g_per_mol
    m_fuel = fuel.M_g_per_mol
    return Stoich(x, 100 * m_ox / (m_ox + m_fuel), 100 * m_fuel / (m_ox + m_fuel), prods)


@dataclass(frozen=True)
class BalanceReport:
    left: dict[str, float]
    right: dict[str, float]
    mismatches: dict[str, tuple[float, float]]

    @property
    def balanced(self) -> bool:
        return not self.mismatches


def check_balance(
    reactants: Sequence[tuple[float, str | Mapping[str, float]]],
    products: Sequence[tuple[float, str | Mapping[str, float]]],
    tol: float = 1e-9,
) -> BalanceReport:
    """반응식의 원소수지를 검사한다. 불일치 원소는 (좌변, 우변) 값으로 보고한다."""

    def tally(side):
        t: dict[str, float] = {}
        for coef, f in side:
            for el, n in as_formula_dict(f).items():
                t[el] = t.get(el, 0.0) + coef * n
        return t

    left, right = tally(reactants), tally(products)
    mism = {
        el: (left.get(el, 0.0), right.get(el, 0.0))
        for el in sorted(set(left) | set(right))
        if abs(left.get(el, 0.0) - right.get(el, 0.0)) > tol
    }
    return BalanceReport(left, right, mism)


@dataclass
class Mixture:
    """질량% 로 정의된 배합. 기준 온도에서의 반응물 혼합물."""

    parts: list[tuple[Ingredient, float]]
    T_K: float = T_REF_K
    _tol: float = field(default=1e-6, repr=False)

    def __post_init__(self) -> None:
        if not self.parts:
            raise ValueError("성분이 없다")
        if any(w < 0 for _, w in self.parts):
            raise ValueError("질량%는 음수일 수 없다")
        total = sum(w for _, w in self.parts)
        if abs(total - 100.0) > self._tol:
            raise ValueError(f"질량% 합이 100이 아니다: {total}")
        keys = [i.key for i, _ in self.parts]
        if len(set(keys)) != len(keys):
            raise ValueError("같은 성분이 중복되었다")

    @classmethod
    def from_wt_pct(cls, wt: Mapping[str, float], db: Mapping[str, Ingredient] | None = None,
                    extra: Sequence[Ingredient] = (), T_K: float = T_REF_K) -> Mixture:
        pool = dict(db if db is not None else load_ingredients())
        pool.update({i.key: i for i in extra})
        parts = []
        for key, w in wt.items():
            if key not in pool:
                raise KeyError(f"DB 에 없는 성분: {key}")
            parts.append((pool[key], float(w)))
        return cls(parts, T_K)

    # --- 조회 -----------------------------------------------------------
    def wt_pct(self) -> list[tuple[str, float]]:
        return [(i.key, w) for i, w in self.parts]

    def ingredient(self, key: str) -> Ingredient:
        for i, _ in self.parts:
            if i.key == key:
                return i
        raise KeyError(key)

    def moles_per_kg(self, key: str) -> float:
        for i, w in self.parts:
            if i.key == key:
                return 10.0 * w / i.M_g_per_mol  # w% of 1000 g -> 10*w g
        raise KeyError(key)

    def mole_ratio(self, a: str, b: str) -> float:
        return self.moles_per_kg(a) / self.moles_per_kg(b)

    def element_moles_per_kg(self) -> dict[str, float]:
        b: dict[str, float] = {}
        for ing, w in self.parts:
            n = 10.0 * w / ing.M_g_per_mol
            for el, cnt in ing.formula.items():
                b[el] = b.get(el, 0.0) + n * cnt
        return b

    def h0_kJ_per_kg(self) -> float:
        """혼합물 비엔탈피 (kJ/kg), 기준 온도 298.15 K 반응물 가정."""
        if abs(self.T_K - T_REF_K) > 1e-9:
            raise NotImplementedError("298.15 K 이외의 반응물 온도는 아직 지원하지 않는다")
        return sum(10.0 * w / i.M_g_per_mol * i.dHf_kJ_per_mol for i, w in self.parts)

    def oxygen_balance_pct(self) -> float:
        """산소수지(%). 양수 = 산화제 과잉, 음수 = 연료 과잉.

        생성물 가정: K 가 있으면 K2CO3, CO2, H2O, N2, 없으면 CO2, H2O, N2.
        100 g 당 부족/잉여 산소 질량(g)을 % 로 표시한다.
        """
        b = self.element_moles_per_kg()
        extra = set(b) - _SUPPORTED_STOICH_ELEMENTS
        if extra:
            raise NotImplementedError(f"산소수지 미지원 원소: {sorted(extra)}")
        need = 2 * b.get("C", 0.0) + 0.5 * b.get("H", 0.0) + 0.5 * b.get("K", 0.0)
        return (b.get("O", 0.0) - need) * ATOMIC_WEIGHTS["O"] / 10.0

    # --- 변형 -----------------------------------------------------------
    def with_moisture(self, fuel_key: str, moisture_wt_pct: float,
                      db: Mapping[str, Ingredient] | None = None) -> Mixture:
        """연료 성분의 일부를 별도 H2O(l) 성분으로 분리한다 (흡습성 고분자용)."""
        if not 0 <= moisture_wt_pct < 100:
            raise ValueError("수분 함량은 0~100 미만이어야 한다")
        water = (db if db is not None else load_ingredients())["WATER"]
        parts = []
        for ing, w in self.parts:
            if ing.key == fuel_key:
                parts.append((ing, w * (1 - moisture_wt_pct / 100)))
                parts.append((water, w * moisture_wt_pct / 100))
            else:
                parts.append((ing, w))
        merged: dict[str, tuple[Ingredient, float]] = {}
        for ing, w in parts:
            merged[ing.key] = (ing, merged[ing.key][1] + w if ing.key in merged else w)
        return Mixture(list(merged.values()), self.T_K)

    def replace_dHf(self, key: str, dHf_kJ_per_mol: float) -> Mixture:
        """한 성분의 생성엔탈피만 바꾼 새 배합 (민감도 분석용)."""
        return Mixture([(i.with_dHf(dHf_kJ_per_mol) if i.key == key else i, w) for i, w in self.parts], self.T_K)

    def provenance_notes(self) -> list[str]:
        """1차 문헌으로 검증되지 않은 입력값 목록."""
        notes = []
        for i, w in self.parts:
            if w > 0 and i.verification != "primary":
                notes.append(f"{i.key}: ΔfH°={i.dHf_kJ_per_mol:g} kJ/mol 은 '{i.verification}' 수준 ({i.source.splitlines()[0][:80]})")
        return notes
