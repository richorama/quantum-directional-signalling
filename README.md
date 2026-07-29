# Quantum Directional Signalling

[![tests](https://github.com/richorama/quantum-directional-signalling/actions/workflows/tests.yml/badge.svg)](https://github.com/richorama/quantum-directional-signalling/actions/workflows/tests.yml)

Exact diamond-norm evaluations of directional signalling for continuous
bipartite unitary families.

For a unitary \(U\) on \(A\otimes B\), the project studies

\[
\delta_A(U)=\inf_{\mathcal E\ \mathrm{CPTP}}
\left\|\operatorname{Tr}_B\operatorname{Ad}_U
-\mathcal E\operatorname{Tr}_B\right\|_\diamond .
\]

This asks how well the reduced output on \(A\) can be predicted from the
reduced input on \(A\) alone, even when the input is correlated with hidden
system \(B\) and an arbitrary reference.

## Plain-language picture

Imagine two quantum objects, A and B, interact. You then discard B and try to
describe A using only what was initially known about A. That description can
fail because B carried hidden information that affected A.

\(\delta_A(U)\) is the smallest unavoidable worst-case error, even after
choosing the best possible A-only rule. Zero means A has autonomous dynamics;
a larger value means B can signal more strongly into A. The diamond norm makes
this an adversarial test that also permits correlations with an external
reference.

The main qualitative result is that the best simplified rule can change as
the interaction gets stronger. For partial SWAP it eventually becomes better
to add deliberate noise than to pretend A evolved unchanged. For XY the exact
curve has three regimes:

![Exact XY directional-signalling curve](paper/xy_curve.svg)

## Results

- **Two-qubit Ising:** for \(U_\theta=e^{-i\theta Z\otimes Z}\),
  \(\delta_A(U_\theta)=|\sin 2\theta|\).
- **Partial SWAP in dimension \(d\):** the optimization reduces exactly to one
  depolarizing parameter. For
  \(\sin\phi\le (d^2-3)/(d^2-1)\),
  \(\delta_A(U_\phi)=2\sin\phi\); at SWAP,
  \(\delta_A=2(1-1/d^2)\).
- **Two-qubit XY/iSWAP:** the weak branch is
  \(\delta_A=\sin4\theta\), while near iSWAP it is
  \(\delta_A=\sin2\theta+\tfrac12\sin^2 2\theta\), giving
  \(\delta_A(\mathrm{iSWAP})=3/2\).
- **Intermediate XY regime:** an exact quartic and unique physical root are
  proved. Convex KKT analysis establishes the \(q_Z=0\) active facet, and an
  exact resultant plus root isolation selects the physical algebraic sheet.

The signalling measure itself and the exact isolated CNOT and SWAP values are
prior work by Barsse, Perinotti, Tosini, and Vaglini. The contribution here is
the continuous-family evaluation, thresholds, optimal effective channels, and
matching witnesses.

## Repository layout

- `quantum_coarse_graining/` — exact Gaussian-rational certificate library.
- `tests/` — focused theorem and regression certificates.
- `paper/aqc_signalling.tex` — manuscript source.
- `paper/aqc_signalling.pdf` — compiled manuscript.
- `tools/derive_xy_middle.py` — optional exact symbolic reconstruction.
- `AUDIT.md` — independent-review checklist.
- `RESEARCH_NOTES.md` — claims ledger, proof boundaries, and open questions.

## Reproduce

The exact certificate suite has no third-party Python dependencies:

```bash
python3 -m unittest discover -s tests -p 'test_aqc_*.py' -v
```

Reconstruct the load-bearing XY eliminations with the optional proof tools:

```bash
python3 -m pip install -e '.[proof]'
python3 tools/derive_xy_middle.py
```

Build the manuscript with:

```bash
cd paper
pdflatex -interaction=nonstopmode -halt-on-error aqc_signalling.tex
pdflatex -interaction=nonstopmode -halt-on-error aqc_signalling.tex
```

Regenerate the XY figure with:

```bash
python3 paper/generate_xy_plot.py
```

## Scope

This is mathematical quantum information, not a proposal for new fundamental
particles, constants, or cosmology. Code supplies exact finite certificates
and regression protection; analytical identities are proved in the manuscript.

## Provenance

The work was developed on the
`approximate-quantum-coarse-graining` branch of
[`richorama/cho`](https://github.com/richorama/cho), through source commit
`a81742b`, then extracted into this focused repository.

## License

Copyright 2026 Richard Astbury. Licensed under
[CC BY 4.0](LICENSE).
