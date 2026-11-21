"""Particelle: classe base e muone."""

from __future__ import annotations

import math
from dataclasses import dataclass

MASSA_MUONE_GEV = 0.10566  # GeV/c^2
TAU_MUONE_S = 2.197e-6  # vita media a riposo, secondi


@dataclass
class Particella:
    """Particella generica: massa e energia totale in GeV, impulso in GeV/c."""

    massa: float
    energia: float
    px: float = 0.0
    py: float = 0.0
    pz: float = 0.0

    def __post_init__(self) -> None:
        if self.massa < 0:
            raise ValueError("la massa non puo' essere negativa")
        if self.energia < self.massa:
            raise ValueError("l'energia totale non puo' essere minore della massa")

    def tipo(self) -> str:
        return "particella"

    @property
    def impulso(self) -> float:
        """Modulo dell'impulso |p| in GeV/c."""
        return math.sqrt(self.px**2 + self.py**2 + self.pz**2)


class Muone(Particella):
    """Muone al suolo.

    Ha la posizione (x, y, z) in metri rispetto all'asse dello sciame. In uno
    sciame verticale l'impulso e' diretto verso il basso: px = py = 0, pz = -|p|.
    """

    tau = TAU_MUONE_S

    def __init__(self, energia: float, x: float = 0.0, y: float = 0.0, z: float = 0.0):
        if energia < MASSA_MUONE_GEV:
            raise ValueError("energia del muone minore della sua massa")
        p = math.sqrt(energia**2 - MASSA_MUONE_GEV**2)
        super().__init__(MASSA_MUONE_GEV, energia, 0.0, 0.0, -p)
        self.x = x
        self.y = y
        self.z = z

    def tipo(self) -> str:
        return "muone"

    @property
    def distanza_asse(self) -> float:
        """Distanza dall'asse dello sciame (piano x-y), in metri."""
        return math.hypot(self.x, self.y)

    def __repr__(self) -> str:
        return (f"Muone(energia={self.energia:.3f} GeV, "
                f"x={self.x:.2f} m, y={self.y:.2f} m, z={self.z:.2f} m)")
