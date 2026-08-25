# Phase-A technical fix log

Date: 2026-08-25 UTC

## Failed run

- Workflow run: https://github.com/YellowJune/SVM-1VM/actions/runs/32816721644
- Job: 97706404175
- Source commit: `3b126a0cdc4b976c54e34fe0a7b876df6d193405`
- Failure location: after all five published-model prediction passes and after CSV creation, while serializing `metrics.json`
- Exception: `TypeError: Object of type bool_ is not JSON serializable`

## Root cause

In `ordering_flip_fraction`, comparisons of NumPy interval scalars produced `numpy.bool_`. Adding that value to the Python integer accumulator promoted the accumulator and returned a NumPy scalar. The frozen numerical calculation was completed, but the resulting gate comparison remained a `numpy.bool_`, which Python's standard JSON encoder rejects.

## Permitted technical correction

Cast the unchanged strict interval-overlap predicate to the built-in Python `bool` before accumulation:

```python
can_reverse = bool((i_min < j_max) and (i_max > j_min))
```

This changes only scalar type, not truth value, selected sites, sampled pairs, model predictions, metrics, gates, or thresholds.

## Integrity rule

The failed run's partial CSV artifact is retained for audit. The corrected run must start from the same pinned archive, source repository commit, five released models, deterministic inactive sample, pair RNG seed, and preregistered thresholds. No Phase-A numerical result from the failed run was printed or observed before this correction was frozen.
