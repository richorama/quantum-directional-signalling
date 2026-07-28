import unittest
from fractions import Fraction

from quantum_coarse_graining.exact import is_unitary
from quantum_coarse_graining.ising import (
    diamond_autonomy_defect,
    flagged_state_witness_certificate,
    hidden_correlation_witness_holds,
    ising_decomposition_holds,
    ising_theorem_certificate,
    ising_unitary,
    product_state_lower_bound_certificate,
    witness_trace_norm_certificate,
)


class IsingAutonomyTheoremTests(unittest.TestCase):
    def test_pythagorean_ising_gates_are_exactly_unitary(self):
        for cosine, sine in (
            (Fraction(1), Fraction(0)),
            (Fraction(4, 5), Fraction(3, 5)),
            (Fraction(3, 5), Fraction(4, 5)),
            (Fraction(0), Fraction(1)),
        ):
            self.assertTrue(is_unitary(ising_unitary(cosine, sine)))

    def test_residual_decomposes_into_the_normalized_leakage_map(self):
        for cosine, sine in (
            (Fraction(1), Fraction(0)),
            (Fraction(4, 5), Fraction(3, 5)),
            (Fraction(3, 5), Fraction(4, 5)),
            (Fraction(0), Fraction(1)),
        ):
            self.assertTrue(ising_decomposition_holds(cosine, sine))

    def test_hidden_correlation_witness_saturates_the_lower_bound(self):
        self.assertTrue(
            hidden_correlation_witness_holds(Fraction(4, 5), Fraction(3, 5))
        )
        self.assertTrue(
            hidden_correlation_witness_holds(Fraction(3, 5), Fraction(4, 5))
        )
        self.assertTrue(
            witness_trace_norm_certificate(Fraction(4, 5), Fraction(3, 5))
        )

    def test_physical_product_states_saturate_the_lower_bound(self):
        self.assertTrue(
            product_state_lower_bound_certificate(
                Fraction(4, 5),
                Fraction(3, 5),
            )
        )

    def test_theorem_certificate_combines_upper_and_lower_bounds(self):
        self.assertTrue(
            ising_theorem_certificate(Fraction(4, 5), Fraction(3, 5))
        )

    def test_signed_witness_has_a_physical_flagged_state_realization(self):
        self.assertTrue(
            flagged_state_witness_certificate(Fraction(4, 5), Fraction(3, 5))
        )

    def test_sharp_diamond_defect(self):
        self.assertEqual(
            diamond_autonomy_defect(Fraction(4, 5), Fraction(3, 5)),
            Fraction(24, 25),
        )
        self.assertEqual(diamond_autonomy_defect(Fraction(1), Fraction(0)), 0)
        self.assertEqual(diamond_autonomy_defect(Fraction(0), Fraction(1)), 0)

    def test_defect_vanishes_at_product_points_not_only_zero_hamiltonian_angle(self):
        self.assertEqual(diamond_autonomy_defect(Fraction(1), Fraction(0)), 0)
        self.assertEqual(diamond_autonomy_defect(Fraction(0), Fraction(1)), 0)

    def test_rejects_non_unit_circle_parameters(self):
        with self.assertRaises(ValueError):
            ising_unitary(Fraction(1, 2), Fraction(1, 2))


if __name__ == "__main__":
    unittest.main()
