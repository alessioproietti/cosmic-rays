import math

import numpy as np
import pytest

from raggicosmici import (Composizione, Generatore, Muone, Nucleo, Particella,
                          campiona_legge_di_potenza)


def composizione_esempio(n_nuclei=5):
    dati = [("Protoni", 1, 1, 2.0), ("Elio", 4, 2, 1.5), ("CNO", 14, 7, 0.8),
            ("Magnesio", 24, 12, 0.5), ("Ferro", 56, 26, 0.2)][:n_nuclei]
    return Composizione("test", 5.0, 1000.0,
                        [Nucleo(nome, a, z, n0, 2.7) for nome, a, z, n0 in dati])


# --- Nucleo -----------------------------------------------------------------
def test_flusso_gamma_uno_e_positivo():
    n = Nucleo("x", 1, 1, 1.0, 1.0)
    assert n.flusso(1.0, 100.0) == pytest.approx(math.log(100.0))


def test_flusso_additivo_negli_intervalli():
    n = Nucleo("x", 1, 1, 2.0, 2.7)
    assert n.flusso(1, 10) + n.flusso(10, 100) == pytest.approx(n.flusso(1, 100))


def test_flusso_integrato_a_infinito():
    n = Nucleo("x", 1, 1, 2.0, 2.7)
    assert n.flusso(1.0, 1e12) == pytest.approx(n.flusso(1.0), rel=1e-9)


def test_flusso_diverge_per_gamma_non_maggiore_di_uno():
    with pytest.raises(ValueError):
        Nucleo("x", 1, 1, 1.0, 1.0).flusso(1.0)


def test_flusso_intervallo_non_valido():
    with pytest.raises(ValueError):
        Nucleo("x", 1, 1, 1.0, 2.7).flusso(10.0, 1.0)


# --- campionamento ----------------------------------------------------------
def test_campionamento_dentro_i_limiti():
    rng = np.random.default_rng(0)
    e = [campiona_legge_di_potenza(rng, 2.7, 5.0, 1000.0) for _ in range(5000)]
    assert min(e) >= 5.0 and max(e) <= 1000.0


def test_campionamento_riproduce_la_media_teorica():
    # <E> = (g-1)/(g-2) * emin per gamma > 2 con emax -> infinito
    rng = np.random.default_rng(1)
    g, emin = 3.5, 2.0
    e = np.array([campiona_legge_di_potenza(rng, g, emin) for _ in range(200_000)])
    assert e.mean() == pytest.approx((g - 1) / (g - 2) * emin, rel=0.03)


# --- Particelle -------------------------------------------------------------
def test_tipo_muone_e_particella():
    assert Muone(1.0).tipo() == "muone"
    assert Particella(0.938, 10.0).tipo() == "particella"


def test_muone_impulso_coerente_con_energia():
    m = Muone(5.0)
    assert math.hypot(m.impulso, m.massa) == pytest.approx(5.0)


def test_particella_con_energia_minore_della_massa():
    with pytest.raises(ValueError):
        Particella(1.0, 0.5)


# --- Generatore -------------------------------------------------------------
@pytest.mark.parametrize("n_nuclei", [1, 3, 5])
def test_generatore_funziona_con_qualsiasi_numero_di_nuclei(n_nuclei):
    gen = Generatore(composizione_esempio(n_nuclei), seed=3)
    for s in gen.genera(50):
        assert s.nucleo in gen.composizione.nuclei


def test_stesso_seed_stessi_risultati():
    a = [s.n_muoni for s in Generatore(composizione_esempio(), seed=42).genera(200)]
    b = [s.n_muoni for s in Generatore(composizione_esempio(), seed=42).genera(200)]
    assert a == b


def test_energia_primaria_nei_limiti():
    comp = composizione_esempio()
    for s in Generatore(comp, seed=5).genera(500):
        assert comp.emin <= s.energia_primaria <= comp.emax


def test_frequenze_dei_nuclei_seguono_i_pesi():
    comp = composizione_esempio()
    gen = Generatore(comp, seed=7)
    n = 40_000
    conteggi = {nu.nome: 0 for nu in comp.nuclei}
    for _ in range(n):
        conteggi[gen.estrai_nucleo().nome] += 1
    for nu, p in zip(comp.nuclei, comp.pesi()):
        assert conteggi[nu.nome] / n == pytest.approx(p, abs=0.01)


def test_energia_dei_muoni_non_supera_quella_del_primario():
    comp = composizione_esempio()
    gen = Generatore(comp, seed=11)
    ferro = comp.nuclei[-1]
    for _ in range(100):
        s = gen.genera_sciame(nucleo=ferro, energia=200.0)
        assert sum(m.energia for m in s.muoni) <= 200.0 + 1e-9
        assert all(m.energia >= gen.soglia_muoni for m in s.muoni)


def test_muoni_distribuiti_in_azimut():
    comp = composizione_esempio()
    gen = Generatore(comp, seed=13)
    s = gen.genera_sciame(nucleo=comp.nuclei[-1], energia=500.0)
    assert s.n_muoni > 20
    x = np.array([m.x for m in s.muoni])
    y = np.array([m.y for m in s.muoni])
    assert y.std() > 0 and x.min() < 0 < x.max()


def test_nessun_muone_sotto_soglia():
    gen = Generatore(composizione_esempio(), soglia_muoni=1.0, seed=2)
    s = gen.genera_sciame(nucleo=gen.composizione.nuclei[0], energia=0.9)
    assert s.n_muoni == 0


def test_soglia_minore_della_massa_del_muone():
    with pytest.raises(ValueError):
        Generatore(composizione_esempio(), soglia_muoni=0.05)
