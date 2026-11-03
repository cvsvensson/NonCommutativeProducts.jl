# Changelog

## 0.6.0

Changes since 0.5.2.

### Breaking

**Products are sorted as they are formed.**
- Define `mul_effect` for each pair of types that can become neighbours, or the first product throws a `MethodError`.
- To leave products unsorted, wrap the code in `with_autosort(false) do ... end`.
- `enable_autosort!` and `disable_autosort!` are no longer public. Delete calls to the first, and replace calls to the second with `with_autosort(false)`.
- `a * b` is now usually an `NCAdd`, so `(a * b).factors` and `convert(NCMul, a * b)` throw. Use `factors` and `coefficient`, which accept a single-term sum.

**`inner`, `dot` and `norm` sort their arguments**, unless autosort is on and factors can have keys.
- Sorting needs `mul_effect` for neighbouring factors, so kets multiplied inside `with_autosort(false)` can now throw a `MethodError`.
- With autosort on, sort a sum formed inside `with_autosort(false)` before calling them if its factors can have keys. Otherwise a zero vector can get a norm around `1e-8`, and `eigsolve` can return wrong eigenvalues.

**`contraction_key` is compared along the leading kets, not just the first.**
- `mul_effect(k', l)` must now be defined for kets whose keys differ.
- A `Number` from it asserts that `k` and `l` are in one subsystem. Across subsystems it must not be a `Number`, not even 0, or `inner` can silently miss terms.
- Any order of kets gives the right result, but one order is fastest. Mixing orders can be tens of times slower.

**The macros throw an `ArgumentError` for more input.**
- `@nc`: `Number` subtypes, `Any`, unions that include a number, and tuple syntax.
- `@commutative`: repeated or overlapping types, types not yet registered with `@nc`, tuple syntax, and a pair in the opposite order to an earlier call. Reordering a pair needs a Julia restart.

**`mul_effect` return values.**
- An invalid return value now throws an `ArgumentError`. This includes `Swap(missing)`, invalid values inside `AddTerms`, and atoms not registered with `@nc`.
- `mul_effect` methods defined by `@commutative` return `Swap(1)` instead of `y * x`.
- `H^0` and `H^1` keep the type of `H` when `mul_effect` returns a wider coefficient or atom type.

**`ncmapreduce` and `ncmap` treat an atom as `1 * atom`**, as they do a product.
- `ncmapreduce(f, (add, mul), atom)` is `mul(scalarmap(1), f(atom))` instead of `f(atom)`.
- `ncmap(f, atom)` is `1 * f(atom)`, so `ncmap(string, atom)` now throws.

### New

- `Swap` and `AddTerms` are public.
- `atom_group`, `foreach_segment`, `split_by_group` and `segment_ranges` split products into runs of factors in the same group. See [Multiple types](docs/src/multiple_types.md).
- `truncate_terms` and `opnorm_bound` drop small terms and return a bound on the error. See [Linear algebra](docs/src/linear_algebra.md).
- `add!!` accepts atoms and numbers in any position, as in `add!!(k1, k2)` and `add!!(1, a * c, 2, 3)`.

The functions in the second and third items are documented but not yet `public`.

### Fixed

- `eigsolve` could return wrong eigenvalues with autosort off: -18.0 for a Hubbard dimer with ground energy -0.83.
- `inner` gave 0 for keyed kets of one state listed in different orders.
- `termpairs` and `map_coefficients` accept atoms.
- `Inf * x` and `NaN * x` keep a zero scalar part zero.
- A zero product is `isequal` only to `zero(prefactor)`, matching `hash`, so `isequal(-(0.0 * a), -0.0)` is false.
- Sums with a `typemin` coefficient print correctly.
- `H^k` for a sum, `ncmap`, and `VectorInterface.add!` with `β = 1` are much faster on large sums.
