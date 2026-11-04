# Changelog

## 0.6.0

Changes since 0.5.2.

### Breaking

**Products are sorted as they are formed.**
- Define `mul_effect` for each pair of types that can become neighbours, or their product throws a `MethodError`.
- To leave products unsorted, wrap the code in `with_autosort(false) do ... end`.
- `enable_autosort!` and `disable_autosort!` are no longer public. Delete calls to the first, and replace calls to the second with `with_autosort(false)`.
- `a * b` is now usually an `NCAdd`, so `(a * b).factors` and `convert(NCMul, a * b)` throw. Use `factors` and `coefficient`, which accept a single-term sum.

**`inner`, `dot` and `norm` sort their arguments**, unless autosort is on and factors can have keys.
- Sorting calls `mul_effect` on neighbouring factors, so kets multiplied inside `with_autosort(false)` can now throw.
- With autosort on, sort a sum formed inside `with_autosort(false)` before calling them if its factors can have keys. Otherwise a zero vector can get a norm around `1e-8`, and `eigsolve` can return wrong eigenvalues.

**`contraction_key` has stricter rules.**
- `mul_effect(k', l)` must be defined for kets whose keys differ, and must be a `Number` when the keys are equal, or `inner` throws.
- A `Number` from it asserts that `k` and `l` are in one subsystem. Across subsystems it must not be a `Number`, not even 0, or `inner` can silently miss terms.
- A pair of terms that doesn't reduce to a scalar throws an `ArgumentError`, unless kets of one subsystem already make it zero.
- Any order of kets gives the right result, but mixing orders can be tens of times slower.

**The macros reject more input with an `ArgumentError`.**
- `@nc`: `Number` subtypes, `Any`, and unions that include a number.
- `@commutative`: repeated or overlapping types, types not yet registered with `@nc`, and the same types in the opposite order to an earlier call. Reordering them needs a Julia restart.

**`ncmap` and `ncmapreduce` treat an atom as `1 * atom`.**
- `ncmap(f, atom)` is `1 * f(atom)`, so `ncmap(string, atom)` now throws.
- `ncmapreduce(f, (add, mul), atom)` is `mul(scalarmap(1), f(atom))` instead of `f(atom)`.

### New

- `Swap` and `AddTerms` are public.
- `add!!` accepts atoms and numbers in any position, as in `add!!(k1, k2)` and `add!!(1, a * c, 2, 3)`.

Documented but not yet `public`:
- `atom_group`, `foreach_segment`, `split_by_group` and `segment_ranges` split products into runs of factors in the same group. See [Multiple types](docs/src/multiple_types.md).
- `truncate_terms` and `opnorm_bound` drop small terms and return a bound on the error. See [Linear algebra](docs/src/linear_algebra.md).

### Fixed

- `eigsolve` could return wrong eigenvalues with autosort off: -18 for a Hubbard dimer with ground energy -0.83.
- `inner` gave 0 for keyed kets of one state listed in different orders.
- `Inf * x` and `NaN * x` no longer turn a zero scalar part of a sum into `NaN`.

### Performance

- `inner` compares `contraction_key` along the leading kets, not just the first, so it skips many more pairs of terms.
- `ncmap` and `VectorInterface.add!` with `β = 1` no longer take time proportional to the whole sum for each term added.
