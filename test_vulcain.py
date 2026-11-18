import numpy as np
from vulcain_v5_complete import (
    build_l1_wheel,
    isolate_ntt_moduli,
    mul_poly_rns_exact,
    vulcain_solve,
)


def test_roue_miroir_L1():
    residues, gaps = build_l1_wheel()
    assert len(residues) == 480
    assert sum(gaps) == 2310


def test_multiplication_exacte_ntt_crt():
    rns = isolate_ntt_moduli(two_adicity=8, num_moduli=2)
    p = list((1, -2, 1))
    p4 = mul_poly_rns_exact(p, p, rns)
    p6 = mul_poly_rns_exact(p4, p, rns)
    assert p6 == list((1, -6, 15, -20, 15, -6, 1))


def test_factorisation_blocs_et_wilkinson():
    # Polynôme de degré 6 combinant racines réelles et paires complexes conjuguées
    poly_mixte = np.convolve(np.convolve(list((1.0, 0.0, 1.0)), list((1.0, -3.0, 2.0))), list((1.0, 2.0, 5.0)))
    res_mixte = vulcain_solve(poly_mixte)
    assert res_mixte["final_backward_error"] < 1e-13

    # Polynôme de Wilkinson W10
    w10 = np.poly(np.arange(1, 11, dtype=float))
    res_w10 = vulcain_solve(w10)
    assert res_w10["final_backward_error"] < 1e-12
