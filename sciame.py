"""Un singolo sciame atmosferico."""

from __future__ import annotations

from typing import List, Optional

import numpy as np

from .nuclei import Composizione, Nucleo
from .particelle import Muone


class Sciame:
    """Sciame generato da un nucleo primario, con i muoni che arrivano al suolo."""

    def __init__(self, nucleo: Nucleo, energia_primaria: float,
                 composizione: Optional[Composizione] = None):
        self.nucleo = nucleo
        self.energia_primaria = energia_primaria  # GeV
        self.composizione = composizione
        self.muoni: List[Muone] = []

    @property
    def A(self) -> float:
        """Numero di massa del primario."""
        return self.nucleo.A

    @property
    def n_muoni(self) -> int:
        return len(self.muoni)

    def aggiungi_muone(self, muone: Muone) -> None:
        self.muoni.append(muone)

    def raggio(self, i: int) -> float:
        """Distanza dall'asse dello sciame del muone i-esimo, in metri."""
        return self.muoni[i].distanza_asse

    def raggi(self) -> np.ndarray:
        """Distanze dall'asse di tutti i muoni, in metri."""
        return np.array([m.distanza_asse for m in self.muoni])

    def __repr__(self) -> str:
        return (f"Sciame(primario={self.nucleo.nome}, "
                f"E={self.energia_primaria:.2f} GeV, muoni={self.n_muoni})")
