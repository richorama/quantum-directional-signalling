"""Smoke tests for the optional symbolic derivation script.

These tests never run the derivations themselves: they only parse the script
and check the packaging contract, so the standard suite stays dependency free.
"""

import ast
import contextlib
import importlib.util
import io
import unittest
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
PROOF_SCRIPT = REPOSITORY_ROOT / "tools" / "derive_xy_middle.py"
PYPROJECT = REPOSITORY_ROOT / "pyproject.toml"

REQUIRED_DERIVATIONS = (
    "check_m3_characteristic",
    "check_normal_derivative",
    "check_zero_multiplier_elimination",
    "check_t_parametrization",
    "check_sheet_selection",
    "check_quartic_discriminant",
    "check_root_window",
    "check_threshold_joining",
    "check_sturm_root_count",
)

# Root finding in the script must stay exact.
FORBIDDEN_CALLS = ("evalf", "nroots", "nsolve", "real_roots")


class ProofScriptSmokeTests(unittest.TestCase):
    def setUp(self):
        self.source = PROOF_SCRIPT.read_text(encoding="utf-8")
        self.tree = ast.parse(self.source, filename=str(PROOF_SCRIPT))

    def test_script_exists_and_parses(self):
        self.assertTrue(PROOF_SCRIPT.is_file())
        self.assertIsInstance(self.tree, ast.Module)
        self.assertTrue(ast.get_docstring(self.tree))

    def test_derivations_are_split_into_documented_functions(self):
        functions = {
            node.name: node
            for node in self.tree.body
            if isinstance(node, ast.FunctionDef)
        }
        for name in REQUIRED_DERIVATIONS:
            self.assertIn(name, functions)
            self.assertTrue(ast.get_docstring(functions[name]), name)
        self.assertIn("main", functions)

    def test_registry_lists_every_derivation(self):
        registry = [
            node
            for node in self.tree.body
            if isinstance(node, ast.AnnAssign)
            and isinstance(node.target, ast.Name)
            and node.target.id == "CHECKS"
        ]
        self.assertEqual(len(registry), 1)
        registered = {
            element.id
            for pair in registry[0].value.elts
            for element in pair.elts
            if isinstance(element, ast.Name)
        }
        self.assertEqual(registered, set(REQUIRED_DERIVATIONS))

    def test_root_counting_stays_exact(self):
        attributes = {
            node.attr for node in ast.walk(self.tree) if isinstance(node, ast.Attribute)
        }
        names = {node.id for node in ast.walk(self.tree) if isinstance(node, ast.Name)}
        for forbidden in FORBIDDEN_CALLS:
            self.assertNotIn(forbidden, attributes)
            self.assertNotIn(forbidden, names)
        self.assertIn("sturm", attributes)

    def test_pyproject_pins_the_optional_dependency(self):
        text = PYPROJECT.read_text(encoding="utf-8")
        self.assertIn("[project.optional-dependencies]", text)
        self.assertIn('proof = ["sympy==', text)

    def test_core_package_never_imports_sympy(self):
        for module in sorted((REPOSITORY_ROOT / "quantum_coarse_graining").glob("*.py")):
            self.assertNotIn("sympy", module.read_text(encoding="utf-8"), module.name)

    @unittest.skipIf(
        importlib.util.find_spec("sympy") is None, "sympy is an optional dependency"
    )
    def test_module_imports_without_running_derivations(self):
        specification = importlib.util.spec_from_file_location(
            "derive_xy_middle", PROOF_SCRIPT
        )
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        self.assertEqual(
            [name for name, _ in module.CHECKS],
            [
                "m3-characteristic",
                "normal-derivative",
                "zero-multiplier",
                "t-parametrization",
                "sheet-selection",
                "quartic-discriminant",
                "root-window",
                "threshold-joining",
                "sturm-count",
            ],
        )
        listing = io.StringIO()
        with contextlib.redirect_stdout(listing):
            self.assertEqual(module.main(["--list"]), 0)
        self.assertIn("sturm-count", listing.getvalue())


if __name__ == "__main__":
    unittest.main()
