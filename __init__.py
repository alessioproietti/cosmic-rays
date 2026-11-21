"""raggicosmici: simulazione Monte Carlo semplificata di sciami atmosferici.

Dai nuclei primari (spettro a legge di potenza) ai muoni secondari al suolo.
Unita': energie in GeV, lunghezze in metri.
"""

from .generatore import Generatore
from .nuclei import Composizione, Nucleo, campiona_legge_di_potenza
from .particelle import MASSA_MUONE_GEV, TAU_MUONE_S, Muone, Particella
from .sciame import Sciame

__version__ = "0.1.0"

__all__ = [
    "Generatore",
    "Composizione",
    "Nucleo",
    "campiona_legge_di_potenza",
    "Particella",
    "Muone",
    "MASSA_MUONE_GEV",
    "TAU_MUONE_S",
    "Sciame",
]
