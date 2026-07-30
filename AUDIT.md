# Independent proof audit

This checklist separates machine-verifiable algebra from mathematical and
literature-review judgments. It is intended for an independent reviewer.

## Reproduce

```bash
python3 -m unittest discover -s tests -p 'test_aqc_*.py' -v
python3 -m pip install -e '.[proof]'
python3 tools/derive_xy_middle.py
python3 tools/verify_block_reductions.py
python3 -m pip install -e '.[numerical]'
python3 tools/check_numerical_sdp.py --full
python3 paper/generate_xy_plot.py  # defect and optimizer phase diagrams
cd paper
pdflatex -interaction=nonstopmode -halt-on-error aqc_signalling.tex
pdflatex -interaction=nonstopmode -halt-on-error aqc_signalling.tex
```

The core test suite has no third-party dependencies. SymPy is optional and is
used only to reconstruct the larger elimination identities.

## Load-bearing XY checks

- [ ] Confirm the linear Watrous SDP in Eq. (watrous-linear) and the covariance
      action on its full feasible pair.
- [ ] Confirm the axial symmetry reduction to channels with Bloch eigenvalues
      `(u,u,w)`.
- [ ] Confirm the Watrous-SDP reduction to the one-parameter invariant witness
      `A`.
- [ ] Run the exact physical-Choi reconstruction for the partial-SWAP and XY
      block spectra.
- [ ] Confirm concavity in `A`, convexity in the channel, and the uniqueness
      arguments in Appendix A.
- [ ] Confirm the facet-selection lemma: the CPTP triangle description, the
      convexity of the fixed-channel norm, uniqueness and continuity of the
      convex-concave saddle, and normal-cone sufficiency for global
      optimality over the whole triangle rather than only the facet.
- [ ] Confirm the Danskin argument excluding the `A=0` witness boundary after
      the first threshold.
- [ ] Confirm the inertia of `M3` on `q_Z=0` and the identity
      `2(alpha + ||M3||_1) = 4 lambda_max`.
- [ ] Run the symbolic derivation and inspect every saturated factor in the
      zero-multiplier elimination.
- [ ] Confirm that the KKT multiplier points into the CPTP triangle.
- [ ] Confirm the `W_+ W_-` sheet factorization and the strict bound excluding
      `W_-`.
- [ ] Confirm the resultant factorization selecting `Q_+`.
- [ ] Confirm the Sturm root count and both threshold joins.
- [ ] Run the independent unreduced CVXPY optimization and compare every
      analytic regime, including a strong partial-SWAP interior point.

## Scope and novelty checks

- [ ] Verify that the signalling measure and exact CNOT/SWAP endpoints are
      attributed to Barsse, Perinotti, Tosini, and Vaglini.
- [ ] Search independently for prior evaluations of the continuous Ising,
      partial-SWAP, and XY/iSWAP families.
- [ ] Check the subsystem-label and normalization bridge to the signalling
      definition of Barsse et al. directly from their equations.
- [ ] Check that no numerical experiment is used as proof of a stated theorem.
- [ ] Record reviewer name, date, software versions, and any disagreements.

## Prior-art search record

Search repeated on 30 July 2026 using the primary texts and metadata for
arXiv:2309.07771 and arXiv:2505.14120, arXiv title/abstract and author
searches, and forward-citation queries.

Search families included:

- directional signalling, causal influence, semicausal and approximately
  semicausal channels;
- diamond distance to autonomous, reduced, or simulated subsystem dynamics;
- partial SWAP, collision models, homogenization, Ising, XY, and iSWAP;
- Cartan gates, operator entanglement, tensor-product-structure distance, and
  quantum model reduction.

No independent continuous-family evaluation of the manuscript's diamond-norm
quantity was found. This is a bounded search result, not proof of novelty:
title/abstract indexes may miss formulas appearing only in full text, and
alternative terminology may evade the listed queries.

Attribution verified from the primary papers:

- Barsse et al. (2024) define the measure and prove
  `S(CNOT) <= 1` while separating signalling from causal influence.
- Barsse et al. (2025) prove the matching CNOT lower bound, the exact
  arbitrary-dimensional SWAP value, the depolarizing SWAP optimizer, and the
  equal-dimension universal upper bound.

## Outstanding before publication

These are open items, not findings. Each must be closed before the work is
submitted or advertised as a stable citable artefact.

- [ ] **ORCID iD.** The author does not yet have one. Register at
      <https://orcid.org> and add it to `CITATION.cff` (`orcid:` under
      `authors`) and to the manuscript author block.
- [x] **Zenodo integration.** Sign in to Zenodo, connect the GitHub account,
      and enable `richorama/quantum-directional-signalling` before creating
      the release. Zenodo uses `CITATION.cff` as the archive metadata when no
      `.zenodo.json` file is present.
- [ ] **Archived DOI.** Deposit a release on Zenodo (or equivalent). After
      Zenodo mints the DOI, add `doi:` to `CITATION.cff` on `main` and add a
      badge in `README.md`; the Zenodo landing page remains the source of
      record for the first archived snapshot.
- [x] **Tagged release.** The archived commit records version `0.1.0` and its
      release date in `CITATION.cff`; the matching GitHub release is tagged
      `v0.1.0`.
- [ ] **arXiv posting.** Include the `paper/*.dat` files in the submission;
      the manuscript builds its figures from them through pgfplots.
- [ ] **Independent review.** The review record below is empty. An unexecuted
      checklist is not evidence of review.

## Computer-assisted proof status

Theorem 4 (the intermediate XY branch) is a computer-assisted result. The
convex analysis selecting the active facet is human-checkable and is stated in
the manuscript as a standalone lemma. The four eliminations that complete the
proof are verified in exact rational arithmetic by
`tools/derive_xy_middle.py`, which is part of the proof of record. A reviewer
must either re-run that script or re-derive the same resultants and Sturm
sequences in an independent computer algebra system. No other theorem in the
manuscript depends on computation.

## Review record

Reviewer:

Date:

Commit:

Outcome:

Notes:
