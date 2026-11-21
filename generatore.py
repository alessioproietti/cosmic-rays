"""Generatore Monte Carlo di sciami."""

from __future__ import annotations

import math
from typing import Iterator, Optional

import numpy as np

from .nuclei import Composizione, Nucleo, campiona_legge_di_potenza
from .particelle import MASSA_MUONE_GEV, Muone
from .sciame import Sciame


class Generatore:
    """Genera sciami da una composizione di primari.

    Il generatore fotografa la composizione al momento della creazione: se la
    modifichi, crea un nuovo Generatore. Con lo stesso `seed` i risultati sono
    riproducibili.

    Parametri
    ---------
    composizione : Composizione
    soglia_muoni : energia minima (GeV) dei muoni generati
    seed : seme del generatore di numeri casuali
    indice_muoni : indice spettrale (assunto) delle energie dei muoni
    """

    def __init__(self, composizione: Composizione, soglia_muoni: float = 1.0,
                 seed: Optional[int] = None, indice_muoni: float = 2.7):
        if soglia_muoni <= MASSA_MUONE_GEV:
            raise ValueError("soglia_muoni deve essere maggiore della massa del muone")
        self.composizione = composizione
        self.soglia_muoni = soglia_muoni
        self.indice_muoni = indice_muoni
        self.rng = np.random.default_rng(seed)
        self._cumulativa = np.cumsum(composizione.pesi())
        self._cumulativa[-1] = 1.0  # evita errori di arrotondamento

    # --- primario -----------------------------------------------------------
    def estrai_nucleo(self) -> Nucleo:
        """Sceglie il nucleo primario con probabilita' proporzionale al suo flusso."""
        i = int(np.searchsorted(self._cumulativa, self.rng.random(), side="right"))
        return self.composizione.nuclei[min(i, len(self._cumulativa) - 1)]

    # --- muoni --------------------------------------------------------------
    def numero_medio_muoni(self, A: float, energia: float) -> float:
        s = self.soglia_muoni
        if energia <= s:
            return 0.0
        re = energia / s
        g = re**0.775 * (1.0 - 1.0 / re) ** 5.96
        return 0.0142 * A * g / s

    def raggio_medio(self, energia: float) -> float:
        """Raggio medio (m) del profilo laterale esponenziale dei muoni."""
        s = self.soglia_muoni
        re = energia / s
        a = 4.2 * s**-0.56
        b = 20.5 * s**-0.94
        return a + b * re**-0.65

    def genera_muoni(self, sciame: Sciame) -> None:
        """Aggiunge allo sciame i muoni al suolo.

        - numero: Poisson attorno al valore medio;
        - posizione: raggio esponenziale, angolo azimutale uniforme;
        - energia: legge di potenza troncata tra soglia e l'energia ancora
          disponibile, cosi' la somma delle energie non supera quella del primario.
        """
        e_prim = sciame.energia_primaria
        s = self.soglia_muoni
        if e_prim <= s:
            return
        n = self.rng.poisson(self.numero_medio_muoni(sciame.A, e_prim))
        r_medio = self.raggio_medio(e_prim)
        residuo = e_prim
        for _ in range(n):
            if residuo <= s:
                break
            e = campiona_legge_di_potenza(self.rng, self.indice_muoni, s, residuo)
            residuo -= e
            r = self.rng.exponential(r_medio)
            phi = self.rng.uniform(0.0, 2.0 * math.pi)
            sciame.aggiungi_muone(Muone(e, r * math.cos(phi), r * math.sin(phi)))

    # --- sciame completo ----------------------------------------------------
    def genera_sciame(self, nucleo: Optional[Nucleo] = None,
                      energia: Optional[float] = None) -> Sciame:
        """Genera uno sciame. Si possono forzare nucleo ed energia (utile nei test)."""
        if nucleo is None:
            nucleo = self.estrai_nucleo()
        if energia is None:
            energia = nucleo.campiona_energia(self.rng, self.composizione.emin,
                                              self.composizione.emax)
        sciame = Sciame(nucleo, energia, self.composizione)
        self.genera_muoni(sciame)
        return sciame

    def genera(self, n: int) -> Iterator[Sciame]:
        """Genera n sciami."""
        for _ in range(n):
            yield self.genera_sciame()
