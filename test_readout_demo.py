import numpy as np
import pytest

import readout_demo as d


def test_shunt_impedance_matches_hand_result():
    assert d.s21_from_abcd(d.shunt(100.0)) == pytest.approx(2 * 100 / (2 * 100 + 50))


def test_series_impedance_matches_hand_result():
    assert d.s21_from_abcd(d.series(100.0)) == pytest.approx(2 * 50 / (2 * 50 + 100))


def test_exact_inductance_makes_the_impedance_real():
    r, cr, ccp = 20e3, 500e-15, 200e-15
    l = d.l_exact(r, cr, ccp)
    z = d.z_cap(d.F0, ccp) + d.z_tank(d.F0, r, l, cr)
    assert abs(z.imag) < 1e-6 * abs(z.real)


@pytest.mark.parametrize("r", [10e3, 30e3])
def test_search_and_formulas_agree(r):
    cr, ccp = 500e-15, 200e-15
    found = d.l_search(r, cr, ccp)
    assert found == pytest.approx(d.l_exact(r, cr, ccp), rel=0.01)
    assert found == pytest.approx(d.l_approx(cr, ccp), rel=0.01)


def depth(cr, r, r_l):
    return d.db(d.s21_single(d.F0, r, d.l_approx(cr, 200e-15), cr, 200e-15, r_l))


def test_lossless_inductor_makes_the_configurations_equal():
    assert depth(450e-15, 100e3, 0.0) == pytest.approx(depth(900e-15, 100e3, 0.0), abs=0.05)


def test_inductor_loss_separates_the_configurations():
    assert abs(depth(450e-15, 500e3, 0.3) - depth(900e-15, 500e3, 0.3)) > 1


def test_stronger_coupling_gives_lower_q_at_high_resistance():
    q = [d.q_factor(d.s21_single(d.FREQ, 1e6, d.l_approx(500e-15, c), 500e-15, c)) for c in (200e-15, 600e-15)]
    assert q[1] < q[0]


def test_parallel_loss_is_equal_to_the_series_loss():
    r_s, l, w = 0.3, 74e-9, 2 * np.pi * d.F0
    r_p = d.parallel_loss(r_s, l)
    l_p = (1 + (r_s / (w * l)) ** 2) * l
    assert 1 / r_p + 1 / (1j * w * l_p) == pytest.approx(1 / (r_s + 1j * w * l))


def test_formulas_agree_with_numerical_two_coupler_result():
    c1, c2, cr, r = 10.3e-15, 59.3e-15, 500e-15, 100e3
    l = 1 / ((2 * np.pi * d.F0) ** 2 * (c1 + c2 + cr))
    numeric = d.db(d.s21_double(d.FREQ, r, l, cr, c1, c2)).max()
    for value in (d.approx_abcd(r, c1, c2), d.exact_published(r, l, cr, c1, c2), d.approx_published(r, c1, c2)):
        assert d.db(value) == pytest.approx(numeric, rel=0.05)
