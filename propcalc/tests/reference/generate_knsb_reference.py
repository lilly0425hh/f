"""KNSB 65/35 회귀 기준값 생성기. 범위 검사 (tests/test_backend_knsb.py) 통과 후에만 실행할 것.

    python tests/reference/generate_knsb_reference.py
"""

import datetime
import json
from pathlib import Path

import cea

from propcalc.backends import get_backend
from propcalc.mixture import Mixture, load_ingredients

PC, PE = 6.895, 0.101325
db = load_ingredients()
mix = Mixture.from_wt_pct({"KNO3": 65, "SORBITOL": 35}, db)
r = get_backend("cea").run(mix, PC, Pe_MPa=PE)

out = {
    "_comment": "CEA 결과 자체를 회귀 기준값으로 저장한 것이다. 문헌 검증값이 아니다.",
    "generated": datetime.date.today().isoformat(),
    "cea_version": cea.__version__,
    "assumptions": r.assumptions(),
    "inputs": {
        "composition_wt_pct": dict(mix.wt_pct()),
        "Pc_MPa": PC, "Pe_MPa": PE, "mode": "equilibrium",
        "dHf_kJ_per_mol": {k: db[k].dHf_kJ_per_mol for k in ("KNO3", "SORBITOL")},
        "verification": {k: db[k].verification for k in ("KNO3", "SORBITOL")},
    },
    "outputs": {
        "Tc_K": r.Tc_K, "c_star_m_s": r.c_star_m_s, "isp_s": r.isp_s, "isp_vac_s": r.isp_vac_s,
        "Cf": r.Cf, "M_eff_g_mol": r.M_eff_g_mol, "MW_g_mol": r.MW_g_mol, "gamma_s": r.gamma_s,
        "area_ratio": r.area_ratio, "Te_K": r.Te_K, "Tt_K": r.Tt_K,
        "condensed_mass_fraction": r.condensed_mass_fraction,
        "exit_condensed_mass_fraction": r.exit_condensed_mass_fraction,
        "species_moles_per_kg_top": dict(r.top_species(10)),
    },
}
path = Path(__file__).with_name("knsb_65_35.json")
path.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
print("wrote", path)
