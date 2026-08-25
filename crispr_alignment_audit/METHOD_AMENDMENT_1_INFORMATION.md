# Method Amendment 1: alignment-information decomposition and certificate

Frozen: 2026-08-25 UTC  
Parent preregistration: `b7fb312d253605509c6cf4ba0cb00a55e5175a33`  
Phase-A run: https://github.com/YellowJune/SVM-1VM/actions/runs/32816721644  
Freeze state: the Phase-A job was still `in_progress` in step 7; no Phase-A metric, prediction, report, artifact, or gate result had been observed.

This amendment strengthens the information-theoretic method definition before any candidate-model training and before observing the published-model instability result. It does **not** change the source data, candidate enumeration, samples, splits, seeds, encoders, comparison methods, metrics, numerical thresholds, or rejection rule in the parent preregistration.

## 1. Site-level latent alignment variables

For a raw guide/target site (x=(g,t)), let (mathcal A(x)) be the complete finite set of co-optimal alignments specified in the parent preregistration. The learned latent model produces

[
q_a=q_phi(amid x),qquad
p_a=p_	heta(Y=1mid a,x),qquad
ar p=sum_{ainmathcal A(x)}q_a p_a.
]

The site is trained once, by the exact Bernoulli likelihood of (ar p). No alignment receives a copied site label as an independent training observation.

For (|mathcal A|=1), all alignment-uncertainty quantities below are defined as zero.

## 2. Two distinct information quantities

### 2.1 Normalized alignment-posterior entropy

[
H_A(x)=
rac{-sum_a q_alog q_a}{log |mathcal A(x)|}.
]

Thus (H_Ain[0,1]). It measures uncertainty about which co-optimal alignment the model uses, but does not by itself show whether the alternatives matter to cleavage prediction.

### 2.2 Alignment-induced predictive information

Let (h(p)=-plog p-(1-p)log(1-p)) be binary entropy in nats. Define

[
I_A(x)
=h(ar p)-sum_a q_a h(p_a)
=sum_a q_a,D_{mathrm{KL}}!left[
operatorname{Ber}(p_a),Vert,operatorname{Ber}(ar p)
ight].
]

By concavity of entropy, (I_A(x)ge 0). It separates two cases that posterior entropy alone conflates:

- high (H_A), low (I_A): several alignments remain plausible but predict essentially the same cleavage risk;
- high (I_A): alignment choice materially changes the predicted biological outcome.

Probabilities are clipped only for logarithms to ([10^{-7},1-10^{-7}]). A computed value in ([-10^{-12},0)) is rounded to zero as floating-point error; any value below (-10^{-12}) is an implementation failure.

## 3. Tie-break sensitivity certificate

Pinsker's inequality gives, for every candidate alignment (a),

[
|p_a-ar p|
le
sqrt{rac{1}{2}
D_{mathrm{KL}}!left[
operatorname{Ber}(p_a),Vert,operatorname{Ber}(ar p)
ight]}
le
sqrt{rac{I_A(x)}{2q_{min}(x)}},
]

where (q_{min}(x)=min_a q_a). Report the clipped certificate

[
C_A(x)=minleft{1,sqrt{rac{I_A(x)}
{2max(q_{min}(x),10^{-12})}}ight}.
]

The certificate may be vacuous when the learned posterior assigns nearly zero mass to an alignment; this limitation must be reported rather than hidden. For every test site, verify numerically that
(max_a|p_a-ar p|le C_A+10^{-7}). Report certificate tightness and the fraction of vacuous certificates (C_A=1).

This is an average-to-worst-case information certificate for the model's finite alignment set, not a guarantee about unenumerated genomic sites or clinical safety.

## 4. Frozen proposed calibrator

Item 7 of Section 3.3 in the parent preregistration is refined as follows.

The proposed proof candidate is the **alignment-information calibrated latent mixture**. It fits one validation-only logistic calibration layer to the features

[
[operatorname{logit}(ar p),,H_A,,I_A].
]

Implementation is frozen to:

- clip (ar p) to ([10^{-6},1-10^{-6}]) before the logit;
- standardize all three features by validation-set mean and population standard deviation, replacing a zero standard deviation by one;
- fit unweighted logistic regression with an intercept, L2 regularization (C=1.0), LBFGS, and at most 1000 iterations;
- fit separately for each training seed using only that seed's guide-disjoint validation predictions;
- apply the frozen validation transformation and coefficients to test and external partitions.

The original entropy-only calibrator ([operatorname{logit}(ar p),H_A]), an information-only calibrator ([operatorname{logit}(ar p),I_A]), and the uncalibrated learned mixture are mandatory ablations. The full three-feature calibrator is the sole proposed method used for the existing keep gates; no post-result selection among calibrators is allowed.

The calibration layer adds four scalar parameters. It is excluded from the shared encoder's +/-5% parameter-budget comparison but included in total parameter and latency reports.

## 5. Required invariance and information checks

In addition to the unchanged parent metrics and gates:

1. permute each evaluated alignment set with ten deterministic SHA-256-derived permutations;
2. require maximum absolute differences in (ar p), (H_A), (I_A), (C_A), and the calibrated prediction to remain <= (10^{-7});
3. verify the identity between the entropy-gap and expected-KL forms of (I_A) to absolute error <= (10^{-7});
4. report (H_A), (I_A), (q_{min}), (C_A), exact prediction range, certificate tightness, and error by partition and activity state;
5. report Spearman correlations of (H_A) and (I_A) with absolute error and with false negatives at the frozen operating points.

These checks add diagnostics; they do not replace or relax any parent keep gate.

## 6. Claim boundary

A successful experiment may support a claim that exact latent-alignment marginalization produces a tie-break-invariant CRISPR off-target activity predictor and exposes a finite-set alignment-information certificate.

It may not support claims that the general mutual-information identity, Pinsker's inequality, latent-variable marginalization, or attention pooling are new. The candidate's novelty remains the task formulation, exact co-optimal CRISPR alignment set likelihood, empirical benchmark, and site-level use of the information decomposition and certificate.

## 7. Decision rule

All Phase-A and Phase-B gates in the parent preregistration remain mandatory and unchanged. Failure of any gate still rejects the candidate. This amendment cannot be used to reinterpret a failing result as a pass.
