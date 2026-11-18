import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

namespace VulcainLean4

def V4 : Nat := 2310

theorem mirror_symmetry_V4 (r : Nat) (hr : r <= V4) :
    Nat.gcd r V4 = 1 → Nat.gcd (V4 - r) V4 = 1 := by
  intro h
  apply Nat.dvd_antisymm
  · have hd1 : Nat.gcd (V4 - r) V4 ∣ V4 - r := Nat.gcd_dvd_left (V4 - r) V4
    have hd2 : Nat.gcd (V4 - r) V4 ∣ V4 := Nat.gcd_dvd_right (V4 - r) V4
    have hd3 : Nat.gcd (V4 - r) V4 ∣ (V4 - r) + r := by
      rw [Nat.sub_add_cancel hr]
      exact hd2
    have hdr : Nat.gcd (V4 - r) V4 ∣ r := (Nat.dvd_add_right hd1).mp hd3
    rw [← h]
    exact Nat.dvd_gcd hdr hd2
  · exact Nat.one_dvd _

theorem goldbach_parametric_identity (S r m q : Nat) (hm : m <= q) :
    let N := 2 * (q * S + r) + 2
    let Q_m := 2 * (r + m * S) + 1
    let P_m := 2 * S * (q - m) + 1
    P_m + Q_m = N := by
  intro N Q_m P_m
  dsimp [N, Q_m, P_m]
  have h : S * (q - m) + m * S = q * S := by
    rw [Nat.mul_sub_left_distrib, Nat.mul_comm m S, Nat.sub_add_cancel (Nat.mul_le_mul_left S hm), Nat.mul_comm S q]
  linarith

theorem quadratic_block_expansion (S1 P1 S2 P2 x : Real) :
    (x ^ 2 + S1 * x + P1) * (x ^ 2 + S2 * x + P2) =
    x ^ 4 + (S1 + S2) * x ^ 3 + (P1 + P2 + S1 * S2) * x ^ 2 + (S1 * P2 + S2 * P1) * x + P1 * P2 := by
  ring

end VulcainLean4
