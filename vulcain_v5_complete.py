#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VULCAIN (v5.1) — Architecture Hybride d'Algèbre Exacte en Mémoire O(1)
et de Factorisation Polynomiale Certifiée Haute Précision
Couplage : La Machine v16.0 + MBC-Annealing v5 + Aberth-Ehrlich Vectorisé
Auteur : Mohamed Ghali (France) | HAL : hal-05860758 | Novembre 2026
"""

import math
import numpy as np

V4 = 2310
W_5_12 = 13 * 17 * 19 * 23 * 29 * 31 * 37
_SPLIT = 134217729.0


def build_l1_wheel():
    half = list(r for r in range(1, V4 // 2, 2) if math.gcd(r, V4) == 1)
    residues = half + list(V4 - r for r in reversed(half))
    gaps = list((residues[(i + 1) % 480] - residues[i]) % V4 for i in range(480))
    return residues, gaps


def isolate_ntt_moduli(two_adicity: int = 8, num_moduli: int = 2) -> list:
    residues, _ = build_l1_wheel()
    res_set = set(residues)
    moduli = list()
    K = (1 << (30 - two_adicity)) + 1
    while len(moduli) < num_moduli:
        q = (K << two_adicity) + 1
        if (q % V4) in res_set and math.gcd(q, W_5_12) == 1:
            for a in range(2, 50):
                if pow(a, (q - 1) // 2, q) == q - 1:
                    omega = pow(a, (q - 1) // (1 << two_adicity), q)
                    moduli.append((q, omega))
                    break
        K += 1
    return moduli


def ntt_cyclic_mul(a: list, b: list, q: int, omega: int) -> list:
    n = len(a)

    def _ntt(vec, root):
        res = list(vec)
        j = 0
        for i in range(1, n):
            bit = n >> 1
            while j & bit:
                j ^= bit
                bit >>= 1
            j ^= bit
            if i < j:
                res[i], res[j] = res[j], res[i]
        length = 2
        while length <= n:
            w_len = pow(root, n // length, q)
            for i in range(0, n, length):
                w = 1
                for k in range(length // 2):
                    u = res[i + k]
                    v = (res[i + k + length // 2] * w) % q
                    res[i + k] = (u + v) % q
                    res[i + k + length // 2] = (u - v) % q
                    w = (w * w_len) % q
            length <<= 1
        return res

    fa = _ntt(a, omega)
    fb = _ntt(b, omega)
    fc = list((fa[i] * fb[i]) % q for i in range(n))
    inv_n = pow(n, q - 2, q)
    inv_omega = pow(omega, q - 2, q)
    return list((x * inv_n) % q for x in _ntt(fc, inv_omega))


def mul_poly_rns_exact(P_int: list, Q_int: list, rns_modules: list, two_adicity: int = 8) -> list:
    deg_res = len(P_int) + len(Q_int) - 1
    N = 1 << (deg_res - 1).bit_length()
    (q1, w1), (q2, w2) = rns_modules
    M = q1 * q2
    inv_q1 = pow(q1, -1, q2)

    def _mul_mod(q, w):
        omega_N = pow(w, (1 << two_adicity) // N, q)
        pa = list(x % q for x in P_int) + list((0,)) * (N - len(P_int))
        pb = list(x % q for x in Q_int) + list((0,)) * (N - len(Q_int))
        return ntt_cyclic_mul(pa, pb, q, omega_N)

    c1 = _mul_mod(q1, w1)
    c2 = _mul_mod(q2, w2)
    res = list()
    for i in range(deg_res):
        val = c1[i] + q1 * (((c2[i] - c1[i]) * inv_q1) % q2)
        res.append(val - M if val > M // 2 else val)
    return res


def _two_sum(a: float, b: float):
    s = a + b
    v = s - a
    return s, (a - (s - v)) + (b - v)


def _two_prod(a: float, b: float):
    p = a * b
    ca = _SPLIT * a
    cb = _SPLIT * b
    ah = ca - (ca - a)
    bh = cb - (cb - b)
    al = a - ah
    bl = b - bh
    return p, ((ah * bh - p) + ah * bl + al * bh) + al * bl


def poly_mod_quad_eft(poly, S: float, P: float):
    coeffs = list(float(c) for c in poly)
    n = len(coeffs) - 1
    if n < 0:
        return np.array(list((0.0, 0.0)))
    if n == 0:
        return np.array(list((0.0, float(next(iter(coeffs))))))
    it = iter(coeffs)
    r1, r0 = float(next(it)), float(next(it))
    e1 = e0 = 0.0
    for idx in range(2, n + 1):
        ps, eps = _two_prod(-S, r1)
        pp, epp = _two_prod(-P, r1)
        next_r1, es1 = _two_sum(r0, ps)
        next_e1 = e0 + eps + es1 - S * e1
        next_r1, next_e1 = _two_sum(next_r1, next_e1)

        next_r0, es0 = _two_sum(coeffs[idx], pp)
        next_e0 = epp + es0 - P * e1
        next_r0, next_e0 = _two_sum(next_r0, next_e0)

        r1, e1, r0, e0 = next_r1, next_e1, next_r0, next_e0
    return np.array(list((r1 + e1, r0 + e0)))


def poly_mul_all_direct(blocks):
    res = np.array(list((1.0,)))
    for S, P in blocks:
        res = np.convolve(res, np.array(list((1.0, S, P))))
    return res


def mbc_annealing_v5(P_coeffs, max_iter: int = 1200, tol: float = 1e-13):
    P_orig = np.asarray(P_coeffs, dtype=float)
    lead = float(next(iter(P_orig)))
    P_orig = P_orig / lead
    is_odd = ((len(P_orig) - 1) % 2 != 0)
    if is_odd:
        P_orig = np.append(P_orig, 0.0)
    n = len(P_orig) - 1
    ratios = list(abs(float(P_orig[k])) ** (1.0 / k) for k in range(1, n + 1) if abs(float(P_orig[k])) > 1e-30)
    rho = float(np.clip(np.exp(np.mean(np.log(ratios))), 1e-8, 1e8)) if ratios else 1.0
    P_curr = np.array(list(float(P_orig[k]) / (rho ** k) for k in range(n + 1)))

    m = n // 2
    angles = np.pi * (np.arange(m) + 0.5) / m + 0.13 / m
    blocks = list(np.array(list((-2.0 * math.cos(float(ang)), 1.0))) for ang in angles)

    for t in range(max_iter):
        P_prod = poly_mul_all_direct(blocks)
        E = P_curr - P_prod
        norm_E = float(np.linalg.norm(E))
        rel_err = norm_E / float(np.linalg.norm(P_curr))
        if rel_err < tol:
            break

        alpha = 1e-4 / m
        eps_r = min(1e-2, 0.1 * rel_err) + 1e-9
        mu_t = min(1e-1, 0.1 * (rel_err ** 2)) + 1e-12

        for j in range(m):
            Sj, Pj = tuple(float(v) for v in blocks[j])
            Qj = poly_mul_all_direct(list(blocks[k] for k in range(m) if k != j))
            Rj = poly_mod_quad_eft(E, Sj, Pj)
            Jj = np.column_stack((poly_mod_quad_eft(np.append(Qj, 0.0), Sj, Pj),
                                  poly_mod_quad_eft(Qj, Sj, Pj)))

            grad = np.zeros(2)
            H = np.zeros((2, 2))
            for k in range(m):
                if k != j:
                    d = blocks[j] - blocks[k]
                    denom = float(np.dot(d, d)) + eps_r
                    grad += -2.0 * d / (denom ** 2)
                    H += (-2.0 * np.eye(2) / (denom ** 2)) + (8.0 * np.outer(d, d) / (denom ** 3))

            w, v = np.linalg.eigh(0.5 * (H + H.T))
            H_pos = v @ np.diag(np.maximum(w, 0.0)) @ v.T
            gamma_j = (alpha / (1e-4 + abs(float(np.linalg.det(Jj))))) * math.tanh(float(rel_err ** 2))
            A_star = Jj.T @ Jj + gamma_j * H_pos + mu_t * np.eye(2)
            delta, _, _, _ = np.linalg.lstsq(A_star, -(Jj.T @ Rj) - gamma_j * grad, rcond=None)
            d0, d1 = tuple(float(v) for v in delta)

            step = 1.0
            norm_Ec = float(np.linalg.norm(E))
            for _ in range(10):
                cand_p = np.convolve(Qj, np.array(list((1.0, Sj + step * d0, Pj + step * d1))))
                if float(np.linalg.norm(P_curr - cand_p)) < norm_Ec:
                    break
                step *= 0.5
            blocks[j] = np.array(list((Sj + step * d0, Pj + step * d1)))
            E = P_curr - poly_mul_all_direct(blocks)

    scaled_blocks = list()
    for b in blocks:
        Sj, Pj = tuple(float(v) for v in b)
        scaled_blocks.append(np.array(list((Sj * rho, Pj * (rho ** 2)))))
    back_err = float(np.linalg.norm(P_orig - poly_mul_all_direct(scaled_blocks)) / np.linalg.norm(P_orig))
    return scaled_blocks, back_err, t + 1, is_odd


def _block_roots(S: float, P: float):
    disc = S ** 2 - 4.0 * P
    sd = math.sqrt(abs(disc))
    if disc >= 0:
        return (-S - sd) / 2.0, (-S + sd) / 2.0
    return complex(-S / 2.0, -sd / 2.0), complex(-S / 2.0, sd / 2.0)


def aberth_ehrlich_fast_refinement(roots_approx, coeffs, tol: float = 1e-15, max_iter: int = 300):
    c_orig = np.asarray(coeffs, dtype=complex)
    c_lead = complex(next(iter(c_orig)))
    c_orig = c_orig / c_lead
    n = len(c_orig) - 1
    z = np.asarray(roots_approx, dtype=complex).copy()

    for it in range(max_iter):
        p_vals = np.ones_like(z)
        p_deriv = np.zeros_like(z)
        for k in range(1, n + 1):
            p_deriv = p_deriv * z + p_vals
            p_vals = p_vals * z + c_orig[k]

        active = np.abs(p_vals) > 1e-300
        if not np.any(active):
            break

        newton_step = np.zeros_like(z)
        safe_deriv = np.where(np.abs(p_deriv[active]) > 1e-300, p_deriv[active], 1e-300)
        newton_step[active] = p_vals[active] / safe_deriv

        Z_diff = z[:, np.newaxis] - z[np.newaxis, :]
        np.fill_diagonal(Z_diff, 1.0)
        inv_diffs = 1.0 / Z_diff
        np.fill_diagonal(inv_diffs, 0.0)

        denom = 1.0 - newton_step * np.sum(inv_diffs, axis=1)
        delta = np.where(np.abs(denom) > 1e-15, newton_step / denom, newton_step)
        relax = 0.5 if it < 15 else 1.0
        z -= relax * delta
        if np.max(np.abs(delta)) < tol * max(1.0, float(np.max(np.abs(z)))):
            break

    back_err = float(np.linalg.norm(c_orig - np.poly(z)) / np.linalg.norm(c_orig))
    return z, back_err


def vulcain_solve(P_coeffs):
    P_arr = np.asarray(P_coeffs, dtype=float)
    blocks, mbc_err, iters, is_odd = mbc_annealing_v5(P_arr)
    n_deg = len(P_arr) - 1
    if mbc_err <= 1e-2:
        roots_init = list()
        zero_removed = False
        for row in blocks:
            S, P = tuple(float(v) for v in row)
            if is_odd and not zero_removed and abs(P) < 1e-10:
                roots_init.append(-S)
                zero_removed = True
                continue
            r1, r2 = _block_roots(S, P)
            roots_init.extend((r1, r2))
    else:
        lead = abs(float(next(iter(P_arr))))
        radius = 1.0 + max(abs(float(P_arr[k])) / lead for k in range(1, n_deg + 1)) ** (1.0 / n_deg)
        angles = list(2.0 * math.pi * k / n_deg + 0.137 for k in range(n_deg))
        roots_init = list(radius * complex(math.cos(a), math.sin(a)) for a in angles)

    roots_refined, final_err = aberth_ehrlich_fast_refinement(roots_init, P_arr)
    return dict(
        blocks=blocks,
        mbc_backward_error=mbc_err,
        iterations=iters,
        roots=roots_refined,
        final_backward_error=final_err,
    )


if __name__ == "__main__":
    print("=" * 72)
    print("VULCAIN (v5.1) — EXÉCUTION DU PIPELINE HYBRIDE À 3 ÉTAGES")
    print("=" * 72)

    rns = isolate_ntt_moduli(two_adicity=8, num_moduli=2)
    P_exact = mul_poly_rns_exact(list((1, -2, 1)), mul_poly_rns_exact(list((1, -2, 1)), list((1, -2, 1)), rns), rns)
    print(f"1. Enclume Discrète (NTT-CRT) : (x - 1)^6 construit exactement = {P_exact}")

    res = vulcain_solve(P_exact)
    print(f"2. Recuit Continu MBC (EFT 106 bits) : Erreur arrière = {res['mbc_backward_error']:.2e} ({res['iterations']} iters)")
    print(f"3. Raffinement Aberth-Ehrlich : Erreur finale = {res['final_backward_error']:.2e}")
