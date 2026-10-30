import pytest

from propcalc.formula import format_formula, molar_mass, parse_formula


def test_simple():
    assert parse_formula("C6H14O6") == {"C": 6, "H": 14, "O": 6}


def test_decimal_and_parentheses_lignin():
    f = parse_formula("C9H8.5O1.9(OCH3)1.0")
    assert f["C"] == pytest.approx(10.0)
    assert f["H"] == pytest.approx(11.5)
    assert f["O"] == pytest.approx(2.9)


def test_nested_parentheses():
    assert parse_formula("Ca(OH)2") == {"Ca": 1, "O": 2, "H": 2}
    assert parse_formula("K4(Fe(CN)6)") == {"K": 4, "Fe": 1, "C": 6, "N": 6}


def test_lignin_molar_mass_matches_literature():
    # BioResources 2020 보고값 M = 178.1
    assert molar_mass("C10H11.5O2.9") == pytest.approx(178.1, abs=0.2)


@pytest.mark.parametrize("bad", ["", "  ", "Xx2", "C6H(14", "C6H14)", "6C", "C0H2", "C-1"])
def test_invalid(bad):
    with pytest.raises(ValueError):
        parse_formula(bad)


def test_format_roundtrip():
    assert format_formula({"O": 6, "C": 6, "H": 14}) == "C6H14O6"
