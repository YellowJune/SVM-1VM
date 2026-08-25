# Preregistration: latent-alignment CRISPR method proof

Frozen: 2026-08-25 UTC  
Novelty-boundary commit: `80c0f55091f1e0ddf4a68475f12c6d2ba5f7274c`  
Status: **frozen before any candidate-model training or alternative-alignment prediction audit**

No numerical gate below may be relaxed after results are observed. Technical fixes may correct downloads, parsing, dependency compatibility, or implementation bugs, but must be logged and must not change samples, splits, metrics, methods, or thresholds.

## 1. Source and task

Use only the pinned CRISPR-Bulge archive:

- URL: https://media.githubusercontent.com/media/OrensteinLab/CRISPR-Bulge/main/files/datasets.zip
- Bytes: 524,400,344
- SHA-256: `f892f70ba4ac3b05b03b2171b4ad38746630de08ad630e650f355dd61203eab0`

Core partitions:

- `CHANGEseq/include_on_targets/CHANGEseq_CR_Lazzarotto_2020_dataset.csv`
- `FullGUIDEseq/include_on_targets/FullGUIDEseq_CR_Lazzarotto_2020_dataset.csv`
- `Refined_TrueOT.csv`

Labels follow the prior feasibility audit: CHANGE-seq and FullGUIDE-seq active at reads >=100 and inactive at reads=0; intermediate reads are excluded. Refined TrueOT uses its binary label.

The unit is one raw sgRNA/off-target site. Existing gaps are removed before candidate alignment enumeration. For raw lengths differing by exactly one, enumerate every one-gap alignment whose gap is outside the final three PAM bases. Retain every alignment having the minimum nucleotide mismatch count, with N as a wildcard. Equal-length sites have a singleton alignment set for the proof model.

## 2. Phase A: published-model tie-break instability

Use the five released TL-GRU-Emb classification ensemble components from the CRISPR-Bulge repository:

`files/bulges/1_folds/5_revision_ensemble_{0..4}_exclude_RHAMPseq_continue_from_change_seq/read_ts_0/cleavage_models/aligned/FullGUIDEseq/classification/c_2/ln_x_plus_one_trans/model_fold_0`

For every active ambiguous site and a deterministic hash sample of 100,000 inactive ambiguous sites, score every co-optimal alignment with every ensemble component and their mean.

For each site, record:

- number and positions of co-optimal alignments;
- min, max, range, standard deviation, and rank of predictions;
- whether the released canonical alignment is among the co-optimal set;
- canonical prediction and ensemble-mean prediction;
- activity, sgRNA, partition, and source signal.

Phase-A problem-existence gates, all required:

1. exact archive bytes and SHA-256 match;
2. all five published model files load and give finite predictions;
3. at least 95% of audited canonical alignments are valid members of the enumerated co-optimal set, or any discrepancy is explained by a documented scoring-definition difference without changing enumeration;
4. at least 10% of active ambiguous sites have ensemble-mean prediction range >=0.05 across co-optimal alignments;
5. median ensemble-mean prediction range among active ambiguous sites is >=0.01;
6. at least 5% of active ambiguous sites change their top/bottom risk ordering relative to another active site under alternative tie-breaks, measured by a deterministic pair sample of 100,000 comparisons;
7. instability appears in at least two of the three partitions.

If any gate fails, reject the candidate without training the proposed method.

## 3. Phase B: five-seed method proof

### 3.1 Split and sampling

CHANGE-seq is the development corpus. Split by complete sgRNA identity, never by row. Guides are sorted by descending active-ambiguous count and greedily assigned to train/validation/test bins targeting 70%/15%/15% of active ambiguous sites; ties are broken by SHA-256 of `"split-v1|" + sgRNA`. The final guide lists and their hashes are written before training.

Use every active CHANGE-seq site. For training only, retain at most 20 inactive sites per active site within each sgRNA and edit class, selected by ascending SHA-256 of the complete raw pair under `"negative-v1|"`. Validation and test metrics use all eligible inactive sites so PR-AUC reflects native prevalence.

External tests use FullGUIDE-seq and Refined TrueOT after excluding any sgRNA present in the CHANGE-seq training bin. Report them separately and pooled; never tune on them.

Seeds: `11, 23, 47, 71, 101`.

### 3.2 Shared alignment encoder

All alignment-input neural methods use the same 24-position, 25-state paired nucleotide/gap encoding and the same parameter budget within +/-5%. The encoder is a small bidirectional GRU followed by a two-layer MLP. Optimizer, learning-rate search space, early-stopping rule, training examples, and class weighting are identical. Hyperparameters are selected once on the validation guides using seed 11, then frozen for all seeds.

### 3.3 Compared methods

1. **Canonical**: one released/canonical alignment per site.
2. **Random-tie**: one uniformly sampled co-optimal alignment per site per epoch; one random alignment at inference, averaged over 20 inference draws.
3. **Uniform-set**: arithmetic mean of alignment-level Bernoulli probabilities; exact permutation invariance.
4. **Max-set**: maximum alignment-level probability; exact permutation invariance.
5. **Raw-pair**: a parameter-matched raw unaligned guide/target cross-attention encoder.
6. **Learned latent mixture**: for alignment set A,
   `q_phi(a|g,t)=softmax(s_phi(a))` and
   `p(y=1|g,t)=sum_{a in A} q_phi(a|g,t) sigmoid(r_theta(a))`.
   Train by exact site-level Bernoulli likelihood. Report posterior entropy
   `H[q_phi(A|g,t)]` in bits.
7. **Entropy-calibrated latent mixture**: the learned latent mixture followed by a validation-fitted scalar calibration using its site logit and normalized posterior entropy. This is the proposed proof candidate.

Uniform-set, max-set, and both latent methods must be numerically invariant to alignment enumeration order.

### 3.4 Metrics

Primary:

- PR-AUC on the ambiguous-only native-prevalence test set;
- binary negative log likelihood on the same set.

Secondary:

- PR-AUC on all sites and on all one-bulge sites;
- Brier score and adaptive 15-bin expected calibration error;
- ROC-AUC only as a secondary diagnostic;
- prediction range under all alternative tie-breaks;
- Spearman correlation between alignment entropy and absolute calibration error;
- single-site and batch CPU latency, peak resident memory, and parameter count.

Report per seed, mean, standard deviation, paired bootstrap 95% confidence intervals, and raw predictions.

### 3.5 Phase-B keep gates

All required:

1. exact invariant methods have max order-induced prediction difference <=1e-7;
2. versus Canonical, the proposed entropy-calibrated latent mixture improves ambiguous-only PR-AUC by >=0.03 absolute **or** >=10% relative, and improves ambiguous-only NLL by >=5% relative;
3. versus the stronger of Uniform-set and Max-set, it improves ambiguous-only PR-AUC by >=0.01 absolute and NLL by >=2% relative;
4. the improvement in gate 2 has the same sign in at least 4/5 seeds and its paired bootstrap 95% interval excludes zero;
5. all-site PR-AUC is non-inferior to Canonical within -0.005 absolute;
6. at least one external partition shows positive PR-AUC and NLL improvement versus Canonical, with neither external partition losing more than 0.01 PR-AUC;
7. median CPU latency is <2 ms/site and <=5x Canonical for the observed one-bulge candidate-set sizes;
8. no train/validation/test sgRNA overlap, no duplicated raw site across splits, and all prediction rows have provenance hashes.

If any gate fails, reject the method and record the negative result. Do not proceed to 10–20 seeds or paper construction.

## 4. Full-study requirements after proof only

A successful proof authorizes, but does not itself satisfy, the final study. The full study must add:

- 10–20 seeds;
- GRU, LSTM, causal Transformer, XGBoost/random forest, and a resource-bounded Gaussian-process or neural-GP baseline;
- guide-disjoint, assay-transfer, and temporally independent evaluations where possible;
- learned-mixture ablations: uniform prior, no entropy calibration, no sequence context in q, all near-optimal alignments, and alignment-count-only;
- actual CPU and available GPU latency with warm-up, batch-size curves, memory, and FLOPs/parameter counts;
- sensitivity to alignment scoring parameters and to one versus multiple bulges;
- calibration, false-negative risk at fixed recall, and subgroup analyses by sgRNA, gap type, gap position, mismatch count, assay, and PAM;
- complete raw CSV/JSON logs, environment lock, checksums, figures, research notes, and an independent rerun.

## 5. Decision rule

Only a full Phase-A pass followed by a full Phase-B pass changes the status to **KEEP FOR FULL STUDY**. Any failure yields **REJECT**, regardless of qualitative appeal.
