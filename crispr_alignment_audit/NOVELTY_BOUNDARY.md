# Frozen novelty boundary: latent-alignment CRISPR off-target activity

Frozen: 2026-08-25 UTC  
Status: **manual literature gate passed provisionally; method proof not yet passed**

## Candidate problem

For one raw sgRNA/off-target sequence pair, a bulge can admit several co-optimal gapped alignments. Existing activity predictors normally receive one precomputed alignment, so the predicted cleavage risk can depend on an arbitrary tie-break even though the biological site is unchanged.

The candidate contribution is **not** a new genome searcher or a new aligner. It is a site-level probabilistic predictor that treats the alignment as a finite latent variable and exactly marginalizes the cleavage likelihood over the complete co-optimal alignment set. The model must also expose posterior alignment entropy and be invariant to enumeration order or to whichever alignment a caller happened to report.

## Empirical feasibility already measured

Pinned source:

- Repository: https://github.com/OrensteinLab/CRISPR-Bulge
- Archive URL: https://media.githubusercontent.com/media/OrensteinLab/CRISPR-Bulge/main/files/datasets.zip
- Bytes: 524,400,344
- SHA-256: `f892f70ba4ac3b05b03b2171b4ad38746630de08ad630e650f355dd61203eab0`
- Preregistered audit run: https://github.com/YellowJune/SVM-1VM/actions/runs/32815515855
- Artifact ID: `9551341189`
- Artifact ZIP SHA-256: `fa6dfc595bb84989c901ab2e99b861446a34b0a56888c6b097847c445cc6631c`

Measured across CHANGE-seq, FullGUIDE-seq, and Refined TrueOT:

- 72,056 active and 8,357,373 inactive sites
- 1,118,659 sites with at least two co-optimal one-bulge alignments
- 2,989 active ambiguous sites
- 6,472 active one-bulge sites, of which 46.1836% are ambiguous
- 98 distinct sgRNAs with active ambiguity
- all three partitions contain at least 10 active ambiguous sites
- 49.3142% of active ambiguous sites have a co-optimal gap-position span of at least two

All frozen numerical feasibility gates passed. These counts authorize only a method proof, not a paper claim.

## Direct-prior disqualification rule

Reject this candidate immediately if a publication, preprint, released method, thesis, or public implementation dated on or before 2026-08-25 does either of the following for CRISPR/Cas off-target **site-level cleavage/activity prediction**:

1. trains or evaluates a predictor by summing, integrating, or otherwise marginalizing a site label over multiple candidate sgRNA–DNA alignments; or
2. accepts a raw guide/target pair or an alignment set and explicitly guarantees invariance to co-optimal alignment choice while predicting cleavage/activity.

A searcher that merely finds loci or reports multiple alignments does not satisfy this direct-prior definition. A fixed-alignment uncertainty model also does not satisfy it.

## Closest primary work and the non-overlap boundary

| Work | What it does | Why it is adjacent rather than direct |
|---|---|---|
| Yaish and Orenstein, NAR 2024, https://doi.org/10.1093/nar/gkae331 | Builds large bulge datasets and GRU predictors from one aligned pair | States that current predictors assume one alignment, multiple optimal alignments may exist, and multiple/no-alignment input is future work |
| SWOffinder, iScience 2024, https://doi.org/10.1016/j.isci.2023.108557 | Exhaustive bulge-aware off-target search | Chooses a reported alignment; does not marginalize a cleavage likelihood |
| CALITAS, CRISPR Journal 2021, https://doi.org/10.1089/crispr.2020.0036 | CRISPR-aware gapped aligner/search | Returns a best alignment by default; no site-activity latent likelihood |
| CRISPRitz, Bioinformatics 2020, https://doi.org/10.1093/bioinformatics/btz867 | Variant-aware off-target search, including bulges | Search/reporting and locus merging, not learned activity marginalization |
| CHOPOFF, bioRxiv 2025, https://doi.org/10.1101/2025.01.06.603201 | Symbolic alignments and rapid guide ranking from off-target counts | Guide-level search/count ranking, not site-level cleavage prediction |
| Sassy, Bioinformatics 2026, https://doi.org/10.1093/bioinformatics/btag244 | Fast exhaustive approximate string search | Enumerates matches/alignments; does not learn cleavage likelihood |
| crispAI, NAR 2024, https://doi.org/10.1093/nar/gkae806 | ZINB predictive uncertainty for cleavage counts | Models assay/count uncertainty from one fixed pair encoding, not alignment uncertainty |
| CRISPR-MBTF, Briefings in Bioinformatics 2026, https://doi.org/10.1093/bib/bbag216 | Multi-branch Transformer with sequence/context modalities | Encodes one fixed 23-nt aligned pair |
| General pair-HMM/forward algorithms and latent biological alignment models | Marginalize sequence alignments in homology/evolution tasks | Establish general machinery, but not the CRISPR site-activity task or its benchmark |

## Search record

Primary-source searches included combinations of:

- `CRISPR off-target "multiple optimal alignments" prediction`
- `CRISPR off-target "alternative alignments" machine learning`
- `CRISPR off-target "alignment ensemble" activity`
- `CRISPR off-target "alignment uncertainty" bulge`
- `CRISPR off-target latent alignment marginalization`
- `CRISPR off-target co-optimal alignment cleavage`
- `CRISPR alignment-free off-target activity prediction`
- 2025–2026 model and preprint searches for CRISPR-MBTF, CRISMER, Guide-Guard, CHOPOFF, Sassy, and Sassy2

No direct prior satisfying the disqualification rule was found. This is a bounded search finding, not proof that no unpublished or unindexed work exists.

## Claim ceiling before experiments

Permitted provisional claim:

> To the best of a frozen primary-source search through 2026-08-25, this is the first CRISPR/Cas9 off-target activity formulation that trains a site-level cleavage likelihood by exact marginalization over all co-optimal one-bulge alignments and evaluates invariance to alignment tie-breaking.

Forbidden claims:

- first method to find multiple CRISPR alignments;
- first bulge-aware CRISPR model;
- first probabilistic or uncertainty-aware CRISPR predictor;
- proof of clinical safety;
- universal alignment-free sequence learning;
- superiority before the frozen proof gates and later full-seed experiments pass.

## Decision

**PROVISIONAL KEEP FOR METHOD PROOF.** If the published fixed-alignment ensemble is not materially unstable, or if learned exact marginalization does not beat both the canonical and simple invariant aggregators under the frozen gates, discard the entire candidate.
