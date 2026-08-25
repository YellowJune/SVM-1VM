# Native coarse-resolution source gate (preregistered)

Date (UTC): 2026-08-25  
Scope: provisional lipid structural-resolution supervision candidate  
Status: fixed before downloading or counting the six discovered lipid-specific GNPS libraries.

## Purpose

The earlier PNNL audit established a large ambiguous fine-label lattice, but it found zero
native coarse labels. Synthetic coarsening alone is not sufficient for adoption under the
project's no-weakness rule. This second audit asks whether public, independently contributed
MS/MS records exist at the native sum-composition resolution in enough quantity and diversity
to support a mixed-resolution learning paper.

## Frozen sources

The discovery run enumerated these public GNPS identifiers:

1. `PNNL-LIPIDS-POSITIVE`
2. `PNNL-LIPIDS-NEGATIVE`
3. `HCE-CELL-LYSATE-LIPIDS`
4. `GNPS-D2-AMINO-LIPID-LIBRARY`
5. `GNPS-N-ACYL-LIPIDS-MASSQL`
6. `GNPS-LIPID-MAPS-STANDARDS-SPECTRA-DB`

All six will be downloaded and hashed. No source may be silently removed after inspection.

## Frozen parsing scope

Target classes are PC, PE, PG, PI, PS, PA, DG, and TG. Ether/plasmalogen and
sphingoid-chain annotations are excluded because their composition algebra is different.
A native coarse record contains exactly one class-level carbon:double-bond pair for a class
whose molecular species requires two or three chains. A record is "mapped" only if its
(class, total carbon, total double bonds) key has at least two distinct fine molecular-species
candidates observed in the independently audited PNNL fine pool.

## Adoption gates

All gates must pass without relaxation:

1. At least **5,000 native coarse MS/MS spectra** in the frozen target classes.
2. At least **150 distinct mapped ambiguous sum-composition keys**.
3. At least **4 target lipid classes** represented among mapped native coarse spectra.
4. At least **80% of native coarse spectra** map to an ambiguous fine candidate set.
5. At least **2 independent GNPS library identifiers** contribute mapped native coarse spectra.
6. No single library contributes more than **90%** of mapped native coarse spectra.
7. Every fetched source has URL, byte count, and SHA-256 recorded.
8. Manual literature audit finds no direct prior that trains an MS2 molecular-species model by
   exactly marginalizing a native mixed-resolution lipid annotation lattice.

## Decision rule

Failure of any numerical gate rejects this candidate. Passing does not prove novelty or paper
quality; it only permits a method proof experiment. Gates will not be weakened after results
are observed.
