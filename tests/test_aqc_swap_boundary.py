import unittest
from fractions import Fraction

from quantum_coarse_graining.exact import (
    Gaussian,
    identity,
    is_unitary,
    matrix_unit,
    scale,
)
from quantum_coarse_graining.swap import (
    coarse_after_swap,
    completely_depolarizing,
    swap_ancilla_witness_certificate,
    swap_decomposition_holds,
    swap_diamond_autonomy_defect,
    swap_residual,
    swap_theorem_certificate,
    swap_unitary,
)


class SwapBoundaryTests(unittest.TestCase):
    def test_swap_is_an_involutive_unitary(self):
        swap = swap_unitary()
        self.assertTrue(is_unitary(swap))
        self.assertEqual(
            coarse_after_swap(matrix_unit(4, 1, 1)),
            matrix_unit(2, 1, 1),
        )

    def test_completely_depolarizing_channel(self):
        self.assertEqual(
            completely_depolarizing(matrix_unit(2, 0, 0)),
            scale(Gaussian(Fraction(1, 2)), identity(2)),
        )
        self.assertEqual(
            completely_depolarizing(matrix_unit(2, 0, 1)),
            scale(Gaussian(0), identity(2)),
        )

    def test_swap_residual_keeps_only_relabelled_b_information(self):
        correlated_basis_state = matrix_unit(4, 1, 1)
        self.assertEqual(
            swap_residual(correlated_basis_state),
            (
                (Gaussian(Fraction(-1, 2)), Gaussian(0)),
                (Gaussian(0), Gaussian(Fraction(1, 2))),
            ),
        )
        self.assertTrue(swap_decomposition_holds())

    def test_ancilla_witness_has_trace_norm_three_halves(self):
        self.assertTrue(swap_ancilla_witness_certificate())
        self.assertEqual(swap_diamond_autonomy_defect(), Fraction(3, 2))

    def test_full_swap_boundary_certificate(self):
        self.assertTrue(swap_theorem_certificate())


if __name__ == "__main__":
    unittest.main()
