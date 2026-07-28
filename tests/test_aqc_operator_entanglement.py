import unittest
from fractions import Fraction

from quantum_coarse_graining.autonomy import (
    closure_defect,
    linear_operator_entanglement,
    normalized_closure_defect,
    operator_entanglement_identity_holds,
)
from quantum_coarse_graining.exact import (
    I,
    ONE,
    ZERO,
    identity,
    is_unitary,
    kron,
    matrix,
    permutation_unitary,
)


def gates():
    z = matrix(((1, 0), (0, -1)))
    x = matrix(((0, 1), (1, 0)))
    cnot = permutation_unitary((0, 1, 3, 2))
    swap = permutation_unitary((0, 2, 1, 3))
    cz = matrix(
        (
            (1, 0, 0, 0),
            (0, 1, 0, 0),
            (0, 0, 1, 0),
            (0, 0, 0, -1),
        )
    )
    iswap = (
        (ONE, ZERO, ZERO, ZERO),
        (ZERO, ZERO, I, ZERO),
        (ZERO, I, ZERO, ZERO),
        (ZERO, ZERO, ZERO, ONE),
    )
    return (identity(4), kron(x, z), cnot, swap, cz, iswap)


class OperatorEntanglementIdentityTests(unittest.TestCase):
    def test_declared_gates_are_unitary(self):
        self.assertTrue(all(is_unitary(gate) for gate in gates()))

    def test_identity_holds_for_two_qubit_gates(self):
        for gate in gates():
            self.assertTrue(operator_entanglement_identity_holds(gate, 2, 2))

    def test_identity_holds_for_unequal_factor_dimensions(self):
        examples = (
            (2, 3, permutation_unitary((0, 2, 4, 1, 3, 5))),
            (3, 2, permutation_unitary((0, 3, 1, 4, 2, 5))),
        )
        for dim_a, dim_b, gate in examples:
            self.assertTrue(is_unitary(gate))
            self.assertTrue(
                operator_entanglement_identity_holds(gate, dim_a, dim_b)
            )

    def test_product_unitary_has_zero_defect_and_entanglement(self):
        product = gates()[1]
        self.assertEqual(closure_defect(product, 2, 2), 0)
        self.assertEqual(linear_operator_entanglement(product, 2, 2), 0)

    def test_known_two_qubit_values(self):
        _, _, cnot, swap, cz, iswap = gates()
        self.assertEqual(closure_defect(cnot, 2, 2), 4)
        self.assertEqual(closure_defect(cz, 2, 2), 4)
        self.assertEqual(closure_defect(swap, 2, 2), 6)
        self.assertEqual(closure_defect(iswap, 2, 2), 6)
        self.assertEqual(normalized_closure_defect(cnot, 2, 2), Fraction(1, 2))
        self.assertEqual(normalized_closure_defect(swap, 2, 2), Fraction(3, 4))


if __name__ == "__main__":
    unittest.main()
