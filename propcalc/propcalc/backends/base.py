"""평형 계산 백엔드 공통 인터페이스."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from ..mixture import Mixture

G0 = 9.80665  # m/s^2, 표준 중력가속도 (정의값)
MODES = ("equilibrium", "frozen")


class ConvergenceError(RuntimeError):
    """백엔드가 수렴하지 못했거나 계산이 중단됨. 결과를 사용하면 안 된다."""


@dataclass
class Result:
    """연소실 평형 + 노즐 팽창 결과. 가정(assumptions)을 항상 함께 가진다."""

    # --- 계산 가정 ---------------------------------------------------------
    backend: str
    mode: str                     # "equilibrium" | "frozen"
    Pc_MPa: float
    Pe_MPa: float
    area_ratio: float             # Ae/At
    composition_wt_pct: dict[str, float]
    # --- 연소실 -----------------------------------------------------------
    Tc_K: float
    M_eff_g_mol: float            # 시스템 질량 / 기체 몰수 (응축상 질량 포함). CEA 의 M(1/n)
    MW_g_mol: float               # 모든 종(응축상 포함) 몰수 기준 평균 분자량. CEA 의 MW
    gamma_s: float                # 연소실 등엔트로피 지수
    species_moles_per_kg: dict[str, float]      # 연소실 종별 몰수 (mol/kg 추진제)
    condensed_species: frozenset[str]           # 응축상 종 이름
    condensed_mass_fractions: dict[str, float]  # 연소실 응축상 질량분율
    # --- 성능 -------------------------------------------------------------
    c_star_m_s: float
    Cf: float                     # Pe = 설정값, 이상팽창(Pa=Pe) 가정의 CEA 추력계수
    isp_s: float                  # Pa = Pe 일 때의 비추력 (= ve/g0)
    isp_vac_s: float              # 진공 비추력
    ve_m_s: float
    # --- 출구/목 -----------------------------------------------------------
    Te_K: float
    gamma_e: float
    M_e_g_mol: float
    Tt_K: float
    exit_condensed_mass_fraction: float
    warnings: list[str] = field(default_factory=list)

    @property
    def condensed_mass_fraction(self) -> float:
        return sum(self.condensed_mass_fractions.values())

    @property
    def gas_mole_fractions(self) -> dict[str, float]:
        """기체 종만으로 정규화한 연소실 몰분율."""
        gas = {k: v for k, v in self.species_moles_per_kg.items() if k not in self.condensed_species}
        total = sum(gas.values())
        return {k: v / total for k, v in gas.items()} if total > 0 else {}

    def assumptions(self) -> str:
        """결과 표기용 계산 가정 한 줄 요약."""
        comp = "/".join(f"{k} {v:g}" for k, v in self.composition_wt_pct.items())
        flow = {"equilibrium": "평형 유동", "frozen": "동결 유동"}[self.mode]
        return (f"{flow}, Pc={self.Pc_MPa:g} MPa, Pe={self.Pe_MPa:.4g} MPa (Ae/At={self.area_ratio:.3f}), "
                f"백엔드={self.backend}, 배합(wt%)={comp}")

    def top_species(self, n: int = 8, include_condensed: bool = True) -> list[tuple[str, float]]:
        items = [(k, v) for k, v in self.species_moles_per_kg.items()
                 if include_condensed or k not in self.condensed_species]
        return sorted(items, key=lambda kv: -kv[1])[:n]


class Backend(ABC):
    name: str = "base"

    @abstractmethod
    def run(
        self,
        mixture: Mixture,
        Pc_MPa: float,
        Pe_MPa: float | None = None,
        area_ratio: float | None = None,
        mode: str = "equilibrium",
    ) -> Result:
        """Pe_MPa 와 area_ratio 중 정확히 하나를 지정한다."""


def validate_run_args(Pc_MPa: float, Pe_MPa: float | None, area_ratio: float | None, mode: str) -> None:
    if mode not in MODES:
        raise ValueError(f"mode 는 {MODES} 중 하나여야 한다: {mode!r}")
    if (Pe_MPa is None) == (area_ratio is None):
        raise ValueError("Pe_MPa 와 area_ratio 중 정확히 하나만 지정해야 한다")
    if Pc_MPa <= 0:
        raise ValueError("Pc 는 양수여야 한다")
    if Pe_MPa is not None and not 0 < Pe_MPa < Pc_MPa:
        raise ValueError("0 < Pe < Pc 이어야 한다")
    if area_ratio is not None and area_ratio <= 1:
        raise ValueError("면적비는 1 보다 커야 한다")
