import pytest

from propcalc.formula import molar_mass
from propcalc.mixture import (
    Ingredient,
    Mixture,
    check_balance,
    load_ingredients,
    stoichiometry,
)

DB = load_ingredients()


def test_db_every_entry_has_source():
    for ing in DB.values():
        assert ing.source and ing.verification


def test_db_rejects_entry_without_source(tmp_path):
    p = tmp_path / "bad.yaml"
    p.write_text(
        "X:\n  name: x\n  formula: {C: 1}\n  dHf_kJ_per_mol: -1.0\n"
        "  u_dHf_kJ_per_mol: 1.0\n  M_g_per_mol: 12.0\n  phase: cr\n"
    )
    with pytest.raises(ValueError, match="source"):
        load_ingredients(p)


def test_db_molar_mass_consistent_with_formula():
    for ing in DB.values():
        assert ing.M_g_per_mol == pytest.approx(molar_mass(ing.formula), abs=0.005), ing.key


def test_knsb_stoichiometry():
    s = stoichiometry(DB["KNO3"], DB["SORBITOL"])
    assert s.mol_ox_per_mol_fuel == pytest.approx(5.2, abs=0.005)
    assert s.mass_pct_oxidizer == pytest.approx(74.3, abs=0.1)
    assert s.mass_pct_fuel == pytest.approx(25.7, abs=0.1)
    assert s.products["K2CO3"] == pytest.approx(2.6, abs=0.005)
    assert s.products["CO2"] == pytest.approx(3.4, abs=0.005)
    assert s.products["H2O"] == pytest.approx(7.0, abs=0.005)
    assert s.products["N2"] == pytest.approx(2.6, abs=0.005)


def test_knsb_65_35_mole_ratio():
    m = Mixture.from_wt_pct({"KNO3": 65, "SORBITOL": 35}, DB)
    ratio = m.mole_ratio("KNO3", "SORBITOL")
    assert ratio == pytest.approx(3.346, abs=0.005)


def test_cellulose_stoichiometry():
    s = stoichiometry(DB["KNO3"], DB["CELLULOSE_COTTON"])
    assert s.mol_ox_per_mol_fuel == pytest.approx(4.8, abs=0.01)
    assert s.mass_pct_oxidizer == pytest.approx(75.0, abs=0.2)


def test_softwood_kraft_lignin_stoichiometry():
    lignin = Ingredient.custom(
        key="LIGNIN_REF", name="softwood kraft lignin (계산 검증 전용)",
        formula="C10H11.5O2.9", dHf_kJ_per_mol=-503.0, u_dHf_kJ_per_mol=None,
        source="test fixture: 지시서 4절 참고 조성 (BioResources 2020), ΔfH°는 검증 대상 아님",
    )
    s = stoichiometry(DB["KNO3"], lignin)
    assert s.mol_ox_per_mol_fuel == pytest.approx(9.1, abs=0.1)
    assert s.mass_pct_oxidizer == pytest.approx(84.0, abs=0.5)


def test_wt_pct_must_sum_to_100():
    with pytest.raises(ValueError, match="100"):
        Mixture.from_wt_pct({"KNO3": 65, "SORBITOL": 30}, DB)


def test_unknown_ingredient():
    with pytest.raises(KeyError):
        Mixture.from_wt_pct({"KNO3": 50, "NOPE": 50}, DB)


def test_bad_memo_reaction_is_detected():
    # 메모에 적힌 잘못된 반응식: 6.25 KNO3 + C6H14O6 -> 3.125 K2CO3 + 2.875 CO2 + 7 H2O + 3.125 N2
    rep = check_balance(
        [(6.25, "KNO3"), (1, "C6H14O6")],
        [(3.125, "K2CO3"), (2.875, "CO2"), (7, "H2O"), (3.125, "N2")],
    )
    assert not rep.balanced
    assert set(rep.mismatches) == {"O"}
    assert rep.left["O"] == pytest.approx(24.75)
    assert rep.right["O"] == pytest.approx(22.125)


def test_correct_reaction_is_balanced():
    rep = check_balance(
        [(5.2, "KNO3"), (1, "C6H14O6")],
        [(2.6, "K2CO3"), (3.4, "CO2"), (7, "H2O"), (2.6, "N2")],
    )
    assert rep.balanced


def test_element_moles_and_enthalpy():
    m = Mixture.from_wt_pct({"KNO3": 65, "SORBITOL": 35}, DB)
    b = m.element_moles_per_kg()
    # 1 kg: KNO3 0.65/0.1011032 mol, 소르비톨 0.35/0.1821718 mol
    n_k = 650 / DB["KNO3"].M_g_per_mol
    n_s = 350 / DB["SORBITOL"].M_g_per_mol
    assert b["K"] == pytest.approx(n_k)
    assert b["C"] == pytest.approx(6 * n_s)
    assert b["O"] == pytest.approx(3 * n_k + 6 * n_s)
    h0 = m.h0_kJ_per_kg()
    expected = n_k * DB["KNO3"].dHf_kJ_per_mol + n_s * DB["SORBITOL"].dHf_kJ_per_mol
    assert h0 == pytest.approx(expected)


def test_oxygen_balance_zero_at_stoichiometric_and_signs():
    s = stoichiometry(DB["KNO3"], DB["SORBITOL"])
    st = Mixture.from_wt_pct({"KNO3": s.mass_pct_oxidizer, "SORBITOL": s.mass_pct_fuel}, DB)
    assert st.oxygen_balance_pct() == pytest.approx(0.0, abs=1e-6)
    rich = Mixture.from_wt_pct({"KNO3": 65, "SORBITOL": 35}, DB)
    assert rich.oxygen_balance_pct() < 0  # 연료 과잉
    lean = Mixture.from_wt_pct({"KNO3": 85, "SORBITOL": 15}, DB)
    assert lean.oxygen_balance_pct() > 0


def test_moisture_split():
    m = Mixture.from_wt_pct({"KNO3": 65, "CELLULOSE_COTTON": 35}, DB).with_moisture("CELLULOSE_COTTON", 6.0, DB)
    w = dict(m.wt_pct())
    assert w["CELLULOSE_COTTON"] == pytest.approx(35 * 0.94)
    assert w["WATER"] == pytest.approx(35 * 0.06)
    assert sum(w.values()) == pytest.approx(100.0)
