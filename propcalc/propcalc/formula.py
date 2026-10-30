"""화학식 파서와 몰질량.

소수 아래첨자와 괄호를 지원한다. 예: "C6H14O6", "C9H8.5O1.9(OCH3)1.0".
"""

from __future__ import annotations

import re
from collections.abc import Mapping

# 표준 원자량 (g/mol). 출처: IUPAC 2005/2007 약식 표준 원자량
# (Wieser, Pure Appl. Chem. 78, 2051 (2006); Na는 2007 개정값).
# 이 모듈이 지원하는 원소만 수록한다. 필요하면 출처와 함께 추가할 것.
ATOMIC_WEIGHTS: dict[str, float] = {
    "H": 1.00794,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "Na": 22.98976928,
    "Mg": 24.3050,
    "Al": 26.9815386,
    "Si": 28.0855,
    "P": 30.973762,
    "S": 32.065,
    "Cl": 35.453,
    "K": 39.0983,
    "Ca": 40.078,
    "Fe": 55.845,
}

_TOKEN = re.compile(r"\s*(?:([A-Z][a-z]?)|(\()|(\))|(\d+\.?\d*|\.\d+))")


def parse_formula(formula: str) -> dict[str, float]:
    """화학식 문자열을 {원소: 원자수} dict로 변환한다."""
    if not isinstance(formula, str) or not formula.strip():
        raise ValueError("화학식이 비어 있다")

    tokens: list[tuple[str, str]] = []
    pos = 0
    text = formula.strip()
    while pos < len(text):
        m = _TOKEN.match(text, pos)
        if not m or m.end() == pos:
            raise ValueError(f"화학식 '{formula}'의 {pos}번째 문자 '{text[pos]}'를 해석할 수 없다")
        if m.group(1):
            tokens.append(("el", m.group(1)))
        elif m.group(2):
            tokens.append(("(", "("))
        elif m.group(3):
            tokens.append((")", ")"))
        else:
            tokens.append(("num", m.group(4)))
        pos = m.end()

    idx = 0

    def parse_group() -> dict[str, float]:
        nonlocal idx
        counts: dict[str, float] = {}
        while idx < len(tokens) and tokens[idx][0] != ")":
            kind, value = tokens[idx]
            if kind == "el":
                if value not in ATOMIC_WEIGHTS:
                    raise ValueError(f"지원하지 않는 원소 '{value}' (화학식 '{formula}')")
                idx += 1
                sub = {value: 1.0}
            elif kind == "(":
                idx += 1
                sub = parse_group()
                if idx >= len(tokens) or tokens[idx][0] != ")":
                    raise ValueError(f"괄호가 닫히지 않았다: '{formula}'")
                idx += 1
            else:
                raise ValueError(f"숫자가 원소나 괄호 뒤가 아닌 위치에 있다: '{formula}'")
            mult = 1.0
            if idx < len(tokens) and tokens[idx][0] == "num":
                mult = float(tokens[idx][1])
                idx += 1
            for el, n in sub.items():
                counts[el] = counts.get(el, 0.0) + n * mult
        return counts

    result = parse_group()
    if idx != len(tokens):
        raise ValueError(f"짝이 맞지 않는 ')': '{formula}'")
    if not result or any(n <= 0 for n in result.values()):
        raise ValueError(f"원자수가 0 이하인 원소가 있다: '{formula}'")
    return result


def as_formula_dict(formula: str | Mapping[str, float]) -> dict[str, float]:
    """문자열 또는 dict를 검증된 원소 dict로 통일한다."""
    if isinstance(formula, str):
        return parse_formula(formula)
    out = {str(k): float(v) for k, v in formula.items()}
    for el, n in out.items():
        if el not in ATOMIC_WEIGHTS:
            raise ValueError(f"지원하지 않는 원소 '{el}'")
        if n < 0:
            raise ValueError(f"원자수가 음수다: {el}={n}")
    return {k: v for k, v in out.items() if v != 0}


def molar_mass(formula: str | Mapping[str, float]) -> float:
    """몰질량 (g/mol)."""
    return sum(ATOMIC_WEIGHTS[el] * n for el, n in as_formula_dict(formula).items())


def format_formula(formula: Mapping[str, float]) -> str:
    """원소 dict를 읽기 쉬운 문자열로 되돌린다 (C, H 우선, 나머지는 알파벳순)."""
    order = [e for e in ("C", "H") if e in formula] + sorted(e for e in formula if e not in ("C", "H"))
    parts = []
    for el in order:
        n = formula[el]
        parts.append(el if n == 1 else f"{el}{n:g}")
    return "".join(parts)
