#!/usr/bin/env python3
"""Demonstration of the method from my bachelor's thesis (*): choose the circuit values that make the
transmitted signal S21 of a bolometer readout resonator respond most strongly to resistance.

The circuit values are round textbook-style examples. The two coupling capacitors are from
Karimi and Pekola, Phys. Rev. Applied 10, 054048 (2018).

(*) New Python code for this repo. The Matlab code that I used in the research group is not here.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths

Z0 = 50.0  # transmission line impedance (ohm)
F0 = 700e6  # target resonance (Hz)
FREQ = F0 + 80e6 * np.linspace(-0.5, 0.5, 1001)  # 660 to 740 MHz in 80 kHz steps


def series(z):
    z = np.asarray(z, dtype=complex)
    return np.stack([np.stack([np.ones_like(z), z], -1), np.stack([np.zeros_like(z), np.ones_like(z)], -1)], -2)


def shunt(z):
    z = np.asarray(z, dtype=complex)
    return np.stack([np.stack([np.ones_like(z), np.zeros_like(z)], -1), np.stack([1 / z, np.ones_like(z)], -1)], -2)


def s21_from_abcd(m):
    """Transmission of a two-port network between two lines of impedance Z0."""
    a, b, c, d = m[..., 0, 0], m[..., 0, 1], m[..., 1, 0], m[..., 1, 1]
    return 2 * Z0 / (a * Z0 + b + c * Z0 ** 2 + d * Z0)


def z_cap(f, c):
    return 1 / (2j * np.pi * f * c)


def z_tank(f, r, l, cr, r_l=0.0):
    """Parallel resonator: resistance r, capacitor cr and an inductor l with series loss r_l."""
    return 1 / (1 / (2j * np.pi * f * l + r_l) + 1 / r + 2j * np.pi * f * cr)


def s21_single(f, r, l, cr, ccp, r_l=0.0):
    """One coupling capacitor to the line: a dip in S21 at resonance."""
    return s21_from_abcd(shunt(z_cap(f, ccp) + z_tank(f, r, l, cr, r_l)))


def s21_double(f, r, l, cr, c1, c2, r_l=0.0):
    """Two coupling capacitors, in series with the line: a peak in S21 at resonance."""
    return s21_from_abcd(series(z_cap(f, c1)) @ shunt(z_tank(f, r, l, cr, r_l)) @ series(z_cap(f, c2)))


def db(s21):
    return 20 * np.log10(np.abs(s21))


def l_approx(cr, ccp, f0=F0):
    """Inductance for resonance at f0 when 1/R is small."""
    return 1 / ((2 * np.pi * f0) ** 2 * (ccp + cr))


def l_exact(r, cr, ccp, f0=F0):
    """Inductance that makes the imaginary part of the total impedance zero at f0 (no inductor loss)."""
    root = np.sqrt((f0 * ccp * np.pi * r) ** 2 - 1)
    return r * (f0 * (ccp + 2 * cr) * np.pi * r - root) / (2 * np.pi * f0 * (4 * np.pi ** 2 * f0 ** 2 * (ccp + cr) * cr * r ** 2 + 1))


def l_search(r, cr, ccp, r_l=0.0, f0=F0, tol=80e3, max_steps=200):
    """The search from the thesis: step L by 5% until f0 is bracketed, then halve the bracket."""
    l = 1 / ((2 * np.pi * f0) ** 2 * cr)
    low_f_l, high_f_l = np.nan, np.nan
    for _ in range(max_steps):
        found = FREQ[np.argmin(np.abs(s21_single(FREQ, r, l, cr, ccp, r_l)))]
        if abs(found - f0) < tol:
            return l
        if found < f0:
            low_f_l, l = l, l * 0.95
        else:
            high_f_l, l = l, l * 1.05
        if not np.isnan(low_f_l + high_f_l):
            l = (low_f_l + high_f_l) / 2
    raise RuntimeError("the inductance search did not converge")


def q_factor(s21, f=FREQ):
    """f_peak / width at half height of the power that the resonator takes from the line."""
    absorbed = 1 - np.abs(s21) ** 2
    peaks, _ = find_peaks(absorbed)
    peak = peaks[np.argmax(absorbed[peaks])]
    return f[peak] / (peak_widths(absorbed, [peak], rel_height=0.5)[0][0] * (f[1] - f[0]))


def parallel_loss(r_s, l, f=F0):
    """The parallel resistance that is equal to a series inductor loss r_s."""
    return (1 + (2 * np.pi * f * l / r_s) ** 2) * r_s


def approx_abcd(r, c1, c2, f0=F0):
    """The thesis approximation of the two-coupler S21 at resonance, from the ABCD matrix."""
    r1 = 1 / ((c1 ** 2 + c2 ** 2) * (2 * np.pi * f0) ** 2 * Z0)
    return -2 * c1 * c2 / (c1 ** 2 + c2 ** 2) / (1 + r1 / r)


def exact_published(r, l, cr, c1, c2, f0=F0):
    """The exact two-coupler S21 at resonance from Karimi and Pekola."""
    w = 2 * np.pi * f0
    top = 2j * w ** 3 * Z0 * c1 * c2 * l
    real = 1 - w ** 2 * (l * (cr + c1 + c2) + l * c2 * Z0 / r)
    imag = 1j * w * (l / r + Z0 * c2 - w ** 2 * l * (c1 + cr) * Z0 * c2)
    return top / (real + imag)


def approx_published(r, c1, c2, f0=F0):
    """The approximate two-coupler S21 at resonance from Karimi and Pekola."""
    r0 = 1 / ((2 * np.pi * f0) ** 2 * Z0 * c2 ** 2)
    return -2 * c1 / c2 / (1 + r0 / r)


def main():

    fig, ax = plt.subplots(figsize=(6, 3.6))
    cr, ccp = 550e-15, 200e-15
    for r in (10e3, 15e3, 20e3, 30e3):
        ax.plot(FREQ / 1e6, db(s21_single(FREQ, r, l_approx(cr, ccp), cr, ccp)), label=f"R = {r / 1e3:.0f} kΩ")
    ax.set(xlabel="f (MHz)", ylabel="S21 (dB)", title="The dip deepens with R, so S21 reads the bolometer")
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/spectroscopy.png", dpi=150)

    rs = np.arange(10e3, 1000e3, 5e3)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharex=True)
    for ax, r_l in zip(axes, (0.0, 0.3)):
        for cr in (450e-15, 600e-15, 750e-15, 900e-15):
            l = l_approx(cr, 200e-15)
            depth = [db(s21_single(FREQ, r, l, cr, 200e-15, r_l)).min() for r in rs]
            ax.plot(rs / 1e3, depth, label=f"Cr = {cr * 1e15:.0f} fF, L = {l * 1e9:.0f} nH")
        ax.set(xlabel="R (kΩ)", ylabel="S21 at the dip (dB)", title=f"Inductor loss = {r_l} Ω")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig("figures/lossy_inductor.png", dpi=150)

    fig, ax = plt.subplots(figsize=(6, 3.6))
    for ccp in (200e-15, 400e-15, 600e-15):
        l = l_approx(500e-15, ccp)
        q = [q_factor(s21_single(FREQ, r, l, 500e-15, ccp, 0.3)) for r in rs]
        ax.plot(rs / 1e3, q, label=f"Ccp = {ccp * 1e15:.0f} fF, L = {l * 1e9:.0f} nH")
    ax.set(xlabel="R (kΩ)", ylabel="Q factor", title="Stronger coupling limits Q")
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/coupling_q.png", dpi=150)

    c1, c2, cr = 10.3e-15, 59.3e-15, 500e-15
    l = 1 / ((2 * np.pi * F0) ** 2 * (c1 + c2 + cr))
    rs = np.arange(10e3, 1000e3, 2e3)
    numeric = np.array([db(s21_double(FREQ, r, l, cr, c1, c2)).max() for r in rs])
    curves = {"Thesis approximation from ABCD": db(approx_abcd(rs, c1, c2)),
              "Published exact": db(exact_published(rs, l, cr, c1, c2)),
              "Published approximation": db(approx_published(rs, c1, c2))}
    fig, ax = plt.subplots(figsize=(6, 3.6))
    for name, values in curves.items():
        error = np.abs(values - numeric) / np.abs(numeric) * 100
        ax.plot(rs / 1e3, error, label=name)
        print(f"{name}: largest difference from the numerical S21 in dB = {error.max():.2f}%")
    ax.set(xlabel="R (kΩ)", ylabel="Difference from numerical (%)", title="Three formulas against the numerical ABCD result")
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/formula_check.png", dpi=150)


if __name__ == "__main__":
    main()
