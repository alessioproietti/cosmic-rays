# raggi-cosmici

Simulazione Monte Carlo **semplificata** di sciami atmosferici da raggi cosmici:
dai nuclei primari (protoni, elio, CNO, magnesio, ferro…) ai muoni secondari che
arrivano al suolo.

È un modello didattico: serve per capire come si costruisce un generatore
Monte Carlo (spettri a legge di potenza, estrazione per inversione della
cumulativa, fluttuazioni di Poisson, profilo laterale), non per fare stime
quantitative di flussi reali. Vedi [Limiti del modello](#limiti-del-modello).

## Cosa fa

1. Sceglie il nucleo primario con probabilità proporzionale al suo flusso
   integrato in `[emin, emax]`.
2. Estrae l'energia del primario da una legge di potenza troncata.
3. Genera i muoni al suolo: numero (Poisson), posizione (raggio esponenziale,
   azimut uniforme) ed energie (con la somma ≤ energia del primario).

Unità: **energie in GeV, lunghezze in metri**.

## Struttura

```
raggicosmici/
  particelle.py   Particella, Muone
  nuclei.py       Nucleo, Composizione, campionamento a legge di potenza
  sciame.py       Sciame (primario + muoni al suolo)
  generatore.py   Generatore (Monte Carlo)
examples/esempio.py      esempio da riga di comando
tests/test_raggicosmici.py
```

## Installazione

```bash
git clone https://github.com/<utente>/raggi-cosmici.git
cd raggi-cosmici
pip install -r requirements.txt
```

Serve Python ≥ 3.9 e NumPy. Il pacchetto si usa direttamente dalla cartella
del repository.

## Uso rapido

```python
from raggicosmici import Composizione, Generatore, Nucleo

# I valori di N0 sono illustrativi, non misure reali.
comp = Composizione("Primari", emin=5.0, emax=1000.0, nuclei=[
    Nucleo("Protoni", A=1,  Z=1,  N0=2.0, gamma=2.7),
    Nucleo("Ferro",   A=56, Z=26, N0=0.2, gamma=2.7),
])

gen = Generatore(comp, soglia_muoni=1.0, seed=42)   # seed -> riproducibile

sciame = gen.genera_sciame()
print(sciame)                    # Sciame(primario=..., E=... GeV, muoni=...)
print(sciame.raggi())            # distanze dei muoni dall'asse [m]

# Forzare primario ed energia (utile per test):
ferro = comp.nuclei[1]
s = gen.genera_sciame(nucleo=ferro, energia=200.0)
print(s.n_muoni, sum(m.energia for m in s.muoni))
```

Esempio completo (statistiche su molti sciami, grafici opzionali):

```bash
python examples/esempio.py --n 10000 --seed 42
python examples/esempio.py --plot          # richiede matplotlib
```

Output (con `--n 10000 --seed 42`):

```
Probabilita' di estrazione dei primari:
  Protoni    0.400
  Elio       0.300
  CNO        0.160
  Magnesio   0.100
  Ferro      0.040

Sciami generati: 10000 (seed 42)
Primari estratti: {'CNO': 1642, 'Protoni': 3982, 'Elio': 2954, 'Ferro': 404, 'Magnesio': 1018}
Muoni totali al suolo: 4265
Muoni per sciame: media 0.43, massimo 27
Sciami senza muoni: 78.0%
Distanza media dei muoni dall'asse: 7.9 m
```

## Modello

**Primari.** Spettro differenziale `dN/dE = N0 · E^(-γ)` per ogni nucleo.
Il flusso integrato tra `E1` ed `E2` è `N0 · (E1^(1-γ) − E2^(1-γ)) / (γ − 1)`
(per `γ = 1`: `N0 · ln(E2/E1)`). L'energia si estrae per inversione della
cumulativa, troncata in `[emin, emax]`.

**Muoni** (soglia `s` = `soglia_muoni`, `E0` = energia del primario, `RE = E0/s`):

```
G        = RE^0.775 · (1 − 1/RE)^5.96
<N_mu>   = 0.0142 · A · G / s            N_mu ~ Poisson(<N_mu>)
R_medio  = 4.2·s^-0.56 + 20.5·s^-0.94 · RE^-0.65      [m]
r ~ Esponenziale(R_medio),   φ ~ Uniforme(0, 2π)
E_mu ~ E^(-2.7) troncata tra s e l'energia ancora disponibile
```

> **TODO:** indicare qui la fonte (libro, articolo o dispense) da cui provengono
> le formule e i coefficienti di `<N_mu>` e `R_medio`. In questo repository
> non sono stati ricalibrati né validati.

## Limiti del modello

- Modello giocattolo: **non** è una simulazione di cascata (nessuna
  propagazione in atmosfera, perdite di energia, interazioni).
- Sciame verticale: nessuna dipendenza dall'angolo zenitale né dal campo
  geomagnetico.
- Il decadimento del muone non è simulato (`TAU_MUONE_S` è solo definita).
- L'indice spettrale delle energie dei muoni (`indice_muoni = 2.7`) è un
  **assunto**, modificabile nel costruttore di `Generatore`.
- I valori di `N0` negli esempi sono illustrativi, non flussi cosmici reali.
- A basse energie la maggior parte degli sciami non produce muoni sopra soglia
  (nell'esempio sopra, circa il 78%): è una conseguenza del modello.

## Test

```bash
pip install -r requirements-dev.txt
python -m pytest
```

## Origine del progetto

> **TODO:** completare, ad esempio: «Sviluppato a partire dalle specifiche
> del corso di *[nome corso]*, *[università]*, a.a. *[anno]*.»

## Licenza

MIT, vedi [LICENSE](LICENSE).
