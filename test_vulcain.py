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


def test_factorisation_racines_multiples_et_wilkinson():
    res_mult = vulcain_solve(list((1, -6, 15, -20, 15, -6, 1)))
    assert res_mult["mbc_backward_error"] < 1e-12
    assert np.max(np.abs(res_mult["roots"] - 1.0)) < 1e-5

    w10 = np.poly(np.arange(1, 11, dtype=float))
    res_w10 = vulcain_solve(w10)
    assert res_w10["final_backward_error"] < 1e-13
