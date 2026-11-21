"""Esempio: genera N sciami e mostra alcune statistiche.

Uso (dalla cartella principale del repository):
    python examples/esempio.py --n 10000 --seed 42
    python examples/esempio.py --plot        # richiede matplotlib
"""

import argparse
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from raggicosmici import Composizione, Generatore, Nucleo  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n", type=int, default=10_000, help="numero di sciami")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--plot", action="store_true", help="mostra gli istogrammi")
    args = ap.parse_args()

    # N0 sono valori illustrativi, NON misure di flussi cosmici reali.
    comp = Composizione("Spettro primario", emin=5.0, emax=1000.0, nuclei=[
        Nucleo("Protoni", A=1, Z=1, N0=2.0, gamma=2.7),
        Nucleo("Elio", A=4, Z=2, N0=1.5, gamma=2.7),
        Nucleo("CNO", A=14, Z=7, N0=0.8, gamma=2.7),
        Nucleo("Magnesio", A=24, Z=12, N0=0.5, gamma=2.7),
        Nucleo("Ferro", A=56, Z=26, N0=0.2, gamma=2.7),
    ])
    print("Probabilita' di estrazione dei primari:")
    for nu, p in zip(comp.nuclei, comp.pesi()):
        print(f"  {nu.nome:<9} {p:6.3f}")

    gen = Generatore(comp, soglia_muoni=1.0, seed=args.seed)
    sciami = list(gen.genera(args.n))

    n_muoni = np.array([s.n_muoni for s in sciami])
    raggi = np.concatenate([s.raggi() for s in sciami if s.n_muoni])
    e_prim = np.array([s.energia_primaria for s in sciami])
    per_nucleo = Counter(s.nucleo.nome for s in sciami)

    print(f"\nSciami generati: {args.n} (seed {args.seed})")
    print("Primari estratti:", dict(per_nucleo))
    print(f"Muoni totali al suolo: {n_muoni.sum()}")
    print(f"Muoni per sciame: media {n_muoni.mean():.2f}, massimo {n_muoni.max()}")
    print(f"Sciami senza muoni: {np.mean(n_muoni == 0):.1%}")
    print(f"Distanza media dei muoni dall'asse: {raggi.mean():.1f} m")

    if args.plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 3, figsize=(14, 4))
        ax[0].hist(e_prim, bins=np.logspace(np.log10(5), 3, 40))
        ax[0].set(xscale="log", yscale="log", xlabel="E primario [GeV]", ylabel="sciami")
        ax[1].hist(n_muoni, bins=40)
        ax[1].set(yscale="log", xlabel="muoni per sciame", ylabel="sciami")
        ax[2].hist(raggi, bins=60)
        ax[2].set(yscale="log", xlabel="distanza dall'asse [m]", ylabel="muoni")
        fig.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
