"""Nuclei primari, spettri a legge di potenza e composizione chimica."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np


def campiona_legge_di_potenza(rng: np.random.Generator, gamma: float,
                              emin: float, emax: float = math.inf) -> float:
    """Estrae E da una distribuzione ~ E^(-gamma) troncata in [emin, emax].

    Usa l'inversione della cumulativa. emax = inf e' ammesso solo se gamma > 1.
    """
    if not 0 < emin < emax:
        raise ValueError("servono 0 < emin < emax")
    r = rng.random()
    if math.isclose(gamma, 1.0):
        if math.isinf(emax):
            raise ValueError("con gamma = 1 serve un emax finito")
        return emin * (emax / emin) ** r
    a = 1.0 - gamma
    if math.isinf(emax) and a >= 0:
        raise ValueError("con emax infinito serve gamma > 1")
    return (emin**a + r * (emax**a - emin**a)) ** (1.0 / a)


@dataclass(frozen=True)
class Nucleo:
    """Nucleo primario con spettro differenziale  dN/dE = N0 * E^(-gamma)."""

    nome: str
    A: float
    Z: float
    N0: float
    gamma: float

    def __post_init__(self) -> None:
        if self.A <= 0:
            raise ValueError("A deve essere positivo")
        if self.N0 < 0:
            raise ValueError("N0 non puo' essere negativo")

    def flusso(self, e1: float, e2: Optional[float] = None) -> float:
        """Flusso integrato: tra e1 ed e2, oppure da e1 a infinito se e2 e' None.

        Le unita' sono quelle di N0 (per GeV^(gamma-1)).
        """
        if e1 <= 0:
            raise ValueError("e1 deve essere positivo")
        gp1 = self.gamma - 1.0
        if e2 is None:
            if gp1 <= 0:
                raise ValueError("il flusso integrato fino a infinito diverge per gamma <= 1")
            return self.N0 * e1 ** (-gp1) / gp1
        if e2 <= e1:
            raise ValueError("serve e2 > e1")
        if math.isclose(self.gamma, 1.0):
            return self.N0 * math.log(e2 / e1)
        return self.N0 * (e1 ** (-gp1) - e2 ** (-gp1)) / gp1

    def campiona_energia(self, rng: np.random.Generator, emin: float, emax: float) -> float:
        """Estrae un'energia primaria dallo spettro, troncato in [emin, emax]."""
        return campiona_legge_di_potenza(rng, self.gamma, emin, emax)


class Composizione:
    """Insieme di nuclei primari nell'intervallo di energia [emin, emax] (GeV)."""

    def __init__(self, nome: str, emin: float, emax: float,
                 nuclei: Optional[Iterable[Nucleo]] = None):
        if not 0 < emin < emax:
            raise ValueError("servono 0 < emin < emax")
        self.nome = nome
        self.emin = emin
        self.emax = emax
        self._nuclei = list(nuclei) if nuclei is not None else []

    @property
    def nuclei(self) -> tuple:
        return tuple(self._nuclei)

    def aggiungi(self, nucleo: Nucleo) -> None:
        self._nuclei.append(nucleo)

    def flusso(self, i: int, e1: float, e2: Optional[float] = None) -> float:
        """Flusso del nucleo di indice i (da 0). Vedi Nucleo.flusso."""
        return self._nuclei[i].flusso(e1, e2)

    def flusso_totale(self) -> float:
        """Somma dei flussi di tutti i nuclei in [emin, emax]."""
        return sum(n.flusso(self.emin, self.emax) for n in self._nuclei)

    def pesi(self) -> np.ndarray:
        """Probabilita' di estrarre ciascun nucleo (flussi in [emin, emax] normalizzati)."""
        if not self._nuclei:
            raise ValueError("la composizione non contiene nuclei")
        flussi = np.array([n.flusso(self.emin, self.emax) for n in self._nuclei])
        totale = flussi.sum()
        if totale <= 0:
            raise ValueError("il flusso totale e' nullo")
        return flussi / totale
