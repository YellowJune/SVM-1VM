# Decision: reject the lipid mixed-resolution candidate

Date (UTC): 2026-08-25  
Decision: **REJECTED — do not develop into a method paper**

## Evidence

The native-source gate was frozen before counting in commit
`adbaa6bf3a8969e4cfac1bdd7da457463dcf7c2b`.

Final audit:

- Workflow run: https://github.com/YellowJune/SVM-1VM/actions/runs/32814950944
- Audit code commit: `29e283d30e73a4a026d01f6a3c2421ef49d5e394`
- Artifact ID: `9551075096`
- Artifact SHA-256: `bd0e138eca6affce8ca63c147a13d548a3f8f8919a8b769a91ead372db7a2720`
- Records downloaded from six frozen source identifiers: 49,957
- Native coarse target-class spectra: **0**
- Mapped ambiguous coarse spectra: **0**
- Mapped ambiguous sum keys: **0**
- Mapped target classes: **0**
- Contributing native-coarse libraries: **0**
- New GNPS LIPID MAPS export: HTTP 404, recorded rather than silently removed

All seven preregistered numerical/reproducibility gates failed.

## Why the earlier positive audit is insufficient

The prior PNNL feasibility audit showed 34,754 fine spectra, 1,142 fine species,
271 ambiguous sum-composition classes, 70.19% of fine spectra in ambiguous classes,
and weighted conditional entropy of 0.7401 bits. Those measurements establish that
the resolution lattice is real. They do **not** establish that public native
mixed-resolution supervision exists. Synthetic label coarsening would therefore test
a constructed weak-supervision setting rather than the proposed native-data problem.

## Consequence

No method proof, baseline training, or paper writing will be performed for this topic.
The candidate is retired rather than weakened. The search moves to a different BioML
problem family.
