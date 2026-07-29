# Independent proof audit

This checklist separates machine-verifiable algebra from mathematical and
literature-review judgments. It is intended for an independent reviewer.

## Reproduce

```bash
python3 -m unittest discover -s tests -p 'test_aqc_*.py' -v
python3 -m pip install -e '.[proof]'
python3 tools/derive_xy_middle.py
python3 paper/generate_xy_plot.py
cd paper
pdflatex -interaction=nonstopmode -halt-on-error aqc_signalling.tex
pdflatex -interaction=nonstopmode -halt-on-error aqc_signalling.tex
```

The core test suite has no third-party dependencies. SymPy is optional and is
used only to reconstruct the larger elimination identities.

## Load-bearing XY checks

- [ ] Confirm the axial symmetry reduction to channels with Bloch eigenvalues
      `(u,u,w)`.
- [ ] Confirm the Watrous-SDP reduction to the one-parameter invariant witness
      `A`.
- [ ] Confirm concavity in `A`, convexity in the channel, and the uniqueness
      arguments in Appendix A.
- [ ] Confirm the inertia of `M3` on `q_Z=0` and the identity
      `2(alpha + ||M3||_1) = 4 lambda_max`.
- [ ] Run the symbolic derivation and inspect every saturated factor in the
      zero-multiplier elimination.
- [ ] Confirm that the KKT multiplier points into the CPTP triangle.
- [ ] Confirm the `W_+ W_-` sheet factorization and the strict bound excluding
      `W_-`.
- [ ] Confirm the resultant factorization selecting `Q_+`.
- [ ] Confirm the Sturm root count and both threshold joins.

## Scope and novelty checks

- [ ] Verify that the signalling measure and exact CNOT/SWAP endpoints are
      attributed to Barsse, Perinotti, Tosini, and Vaglini.
- [ ] Search independently for prior evaluations of the continuous Ising,
      partial-SWAP, and XY/iSWAP families.
- [ ] Check that no numerical experiment is used as proof of a stated theorem.
- [ ] Record reviewer name, date, software versions, and any disagreements.

## Review record

Reviewer:

Date:

Commit:

Outcome:

Notes:
