# E5-C7 age-information leverage audit

## Purpose

E5-C7 asks a value-of-information question that is different from the main present-state audit:

> If an independent age measurement restricted the model-conditional matching age to an interval of total width Δt, how wide could the cumulative-forcing identified set still be in the worst case?

The calculation uses only the frozen derived histories already selected by the paper's present-state constraints.

## Definition

For a model slice with retained histories `h`, matching ages `t(h)` and a positive historical functional `Q(h)`, E5-C7 evaluates the largest multiplicative max/min width that can occur inside any age interval of total width Δt.

Operationally, because the support is finite, this is the maximum `Q_i / Q_j` ratio over pairs of histories whose matching ages differ by no more than Δt.

## Interpretation

The result is:

- a deterministic finite-support sensitivity diagnostic;
- conditional on the adopted stellar-model slice;
- worst-case over the unknown location of the age interval.

It is **not**:

- a posterior credible interval;
- a claim that a specified M-dwarf age precision is currently attainable;
- a population statement about M dwarfs in general.

## Key result

At a 100-Myr total age-window width:

- TOI-700 0.40 Msun: cumulative-EUV width = 1.119;
- TOI-700 0.45 Msun: cumulative-EUV width = 1.115;
- LHS 1140 0.20 Msun: cumulative-EUV width = 1.064;
- pooled TOI-700 0.40/0.45 Msun support: cumulative-EUV width = 1.731.

Thus precise orthogonal age information can sharply contract historical support **within** a model slice, while age information alone need not remove ambiguity between model slices. This is a positive control on the identifiability audit: the method is not constructed to return historical degeneracy for every constraint set.
