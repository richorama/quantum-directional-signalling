"""Exact Pauli matrices over the Gaussian rationals."""

from .exact import I, ZERO, identity, matrix

PAULI_I = identity(2)
PAULI_X = matrix(((0, 1), (1, 0)))
PAULI_Y = ((ZERO, -I), (I, ZERO))
PAULI_Z = matrix(((1, 0), (0, -1)))
PAULIS = (PAULI_I, PAULI_X, PAULI_Y, PAULI_Z)
