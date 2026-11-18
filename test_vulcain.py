def test_factorisation_racines_multiples_et_wilkinson():
    res_mult = vulcain_solve(list((1, -6, 15, -20, 15, -6, 1)))
    assert res_mult["final_backward_error"] < 1e-12
    assert np.max(np.abs(res_mult["roots"] - 1.0)) < 1e-2

    w10 = np.poly(np.arange(1, 11, dtype=float))
    res_w10 = vulcain_solve(w10)
    assert res_w10["final_backward_error"] < 1e-12
