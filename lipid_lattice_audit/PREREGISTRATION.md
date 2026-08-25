# Lipid structural-resolution lattice: feasibility preregistration

Date: 2026-08-25 (UTC)
Status: provisional candidate; not adopted

## Candidate problem

Tandem-MS lipid annotations carry different amounts of structural information. For example,
`PC 34:1` constrains total carbons and double bonds, while `PC 16:0_18:1` identifies an
unordered molecular species and slash-separated notation may additionally assert sn-position.
These are constraint sets over latent lipid components, not ordinary flat class labels.

The provisional method would retain every annotation at its experimentally supported resolution
and maximize the probability mass of all fine structures compatible with it. No hard promotion
of a sum-composition label to one molecular species is permitted.

## Why this audit is necessary

LipiDetective (Würf et al., Briefings in Bioinformatics, 2026) reports that one source was
originally annotated at sum-species level and then upgraded to a molecular species by selecting
the candidate with the greatest summed matching-fragment intensity. The authors explicitly note
residual ambiguity. This audit asks whether the underlying constraint classes are large and
diverse enough for a real learning problem.

## Fixed scope

Primary feasibility data: public GNPS PNNL-LIPIDS-POSITIVE and PNNL-LIPIDS-NEGATIVE
reference libraries. The main gate uses plain-acyl PC, PE, PG, PI, PS, PA, DG, and TG labels.
Ether/plasmalogen and sphingoid annotations are excluded from the main gate because their
additional chemistry must not be collapsed into the same sum-composition class.

## Gates fixed before download

1. At least 20,000 parseable fine-resolution spectra.
2. At least 750 distinct fine molecular species.
3. At least 250 distinct sum-composition equivalence classes with two or more fine candidates.
4. At least 30% of fine-resolution spectra lie in ambiguous sum-composition classes.
5. Median candidate count among ambiguous classes is at least 2.
6. Spectrum-weighted empirical conditional entropy H(fine species | class, sum composition)
   is at least 0.50 bits.
7. Source URLs, byte counts, and SHA-256 hashes are recorded.
8. Literature gate: no direct prior method may already train an MS model by exact
   constraint-marginalization of mixed lipid structural-resolution labels.

Failure of any numerical gate rejects this candidate without relaxing thresholds. Passing the
audit does not establish novelty or performance; it only authorizes implementation.
