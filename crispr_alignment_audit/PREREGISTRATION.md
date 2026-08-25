# Preregistered feasibility gate: latent alignment ambiguity in CRISPR off-target prediction

Date (UTC): 2026-08-25  
Status: provisional BioML candidate; no method name or adoption claim  
Decision policy: every numerical and literature gate must pass without relaxation.

## Candidate problem

Current bulge-aware CRISPR/Cas9 off-target predictors consume one aligned sgRNA–DNA
pair. When several alignments have the same optimal edit cost, preprocessing selects one
representative and the neural model is therefore a function of an arbitrary tie-break.

The provisional method direction is **not** another sequence encoder. It is an
alignment-latent cleavage model whose likelihood is exactly marginalized over every
biologically admissible co-optimal alignment, with alignment entropy exposed as an
information-theoretic uncertainty variable.

## Known prior boundary

- CALITAS (2021) is a CRISPR-aware aligner and can report multiple alignments, but it
  does not train a cleavage-activity predictor by marginalizing them.
- SWOffinder/CRISPR-Bulge (NAR 2024) emits the alignment with fewer bulges when
  multiple optimal alignments exist. Its paper explicitly states that present predictors
  assume one alignment and proposes multiple-alignment or alignment-free analysis as
  future work.
- CHOPOFF (bioRxiv 2025) uses symbolic alignments for exhaustive off-target search and
  rapid guide-level ranking from off-target counts. It is not an alignment-marginalized
  site-level cleavage model.
- The manual novelty gate remains open until a broader 2024–2026 audit is completed.

Primary URLs:

- https://doi.org/10.1093/nar/gkae331
- https://doi.org/10.1089/crispr.2020.0036
- https://doi.org/10.1101/2025.01.06.603201

## Frozen public source

Repository: https://github.com/OrensteinLab/CRISPR-Bulge  
License: MIT  
LFS object: `files/datasets.zip`  
Expected bytes: **524,400,344**  
Expected SHA-256: `f892f70ba4ac3b05b03b2171b4ad38746630de08ad630e650f355dd61203eab0`

The audit will use the three core partitions referenced by the public loader code:

1. CHANGE-seq
2. Full GUIDE-seq
3. Refined TrueOT

It will report all discovered archive paths and freeze the exact resolved paths in the
machine-readable artifact. No dataset may be silently substituted after counting.

## Frozen ambiguity definition

For each row, remove existing gap symbols from `Align.sgRNA` and
`Align.off-target`. For raw sequences whose lengths differ by exactly one nucleotide,
enumerate every single-gap alignment whose gap is inside the protospacer (not the final
three PAM positions). Treat `N` in the sgRNA PAM as a wildcard. An alignment is
co-optimal when it attains the minimum mismatch count among this complete enumerated set.
A site is alignment-ambiguous only when at least two **distinct** alignments are co-optimal.

For CHANGE-seq and Full GUIDE-seq, an active site has `reads >= 100`, matching the
published default training threshold. For Refined TrueOT, use the published binary label
when present; otherwise use positive read/editing signal. All definitions and fallbacks
must be logged.

## Frozen gates

All must pass:

1. Source byte count and SHA-256 exactly match the pinned LFS object.
2. The three core partitions are present with required sequence and activity fields.
3. At least **50,000 active sites** and **500,000 inactive sites** are available overall.
4. At least **1,000 total sites** have two or more co-optimal alignments.
5. At least **100 active sites** have two or more co-optimal alignments.
6. At least **15 distinct sgRNAs** contribute active ambiguous sites.
7. Active ambiguous sites occur in at least **2 source partitions**, with at least
   **10 active ambiguous sites per contributing partition**.
8. At least **5% of active bulge sites** are alignment-ambiguous.
9. At least **25% of active ambiguous sites** have co-optimal gap positions spanning
   two or more nucleotide positions, ensuring a nontrivial positional representation change.
10. Required source files, archive manifest, counts, and hashes are fully recorded.
11. Manual literature search finds no prior site-level cleavage predictor trained by exact
    marginalization over all co-optimal CRISPR bulge alignments.

## Decision rule

Failure of any numerical gate rejects this candidate. Passing only authorizes a small
method proof with guide-disjoint splits; it does not establish novelty or superiority.
