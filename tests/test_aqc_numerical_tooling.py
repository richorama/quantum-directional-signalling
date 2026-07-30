"""Smoke tests for the optional independent numerical SDP checker."""

import ast
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPOSITORY_ROOT / "tools" / "check_numerical_sdp.py"
PYPROJECT = REPOSITORY_ROOT / "pyproject.toml"


class NumericalToolingSmokeTests(unittest.TestCase):
    def setUp(self):
        self.source = SCRIPT.read_text(encoding="utf-8")
        self.tree = ast.parse(self.source, filename=str(SCRIPT))

    def test_script_exists_and_parses(self):
        self.assertTrue(SCRIPT.is_file())
        self.assertTrue(ast.get_docstring(self.tree))

    def test_checker_solves_the_unreduced_problem(self):
        functions = {
            node.name
            for node in self.tree.body
            if isinstance(node, ast.FunctionDef)
        }
        self.assertIn("signalling_sdp", functions)
        self.assertIn("reduced_channel_choi", functions)
        self.assertIn("effective_after_trace_choi", functions)

    def test_checker_does_not_import_the_exact_package(self):
        imports = [
            node
            for node in ast.walk(self.tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
        ]
        self.assertFalse(
            any(
                isinstance(node, ast.ImportFrom)
                and node.module
                and node.module.startswith("quantum_coarse_graining")
                for node in imports
            )
        )

    def test_pyproject_pins_numerical_dependencies(self):
        text = PYPROJECT.read_text(encoding="utf-8")
        self.assertIn('cvxpy==1.5.2', text)
        self.assertIn('scipy==1.13.1', text)


if __name__ == "__main__":
    unittest.main()
