# Changelog

## 0.6.0

Changes since 0.5.2.

### Breaking

**Products are sorted as they are formed.**
- Define `mul_effect` for each pair of types that can become neighbours, or the first product throws a `MethodError`.
- To leave products unsorted, wrap the code in `with_autosort(false) do ... end`.
- `enable_autosort!` and `disable_autosort!` are no longer public. Delete calls to the first, and replace calls to the second with `with_autosort(false)`.
- `a * b` is now usually an `NCAdd`, so `(a * b).factors` and `convert(NCMul, a * b)` throw. Use `factors` and `coefficient`, which accept a single-term sum.

**`inner`, `dot` and `norm` sort their arguments**, also inside `with_autosort(false)`.
- Kets multiplied inside `with_autosort(false)` now need a ket-ket `mul_effect`, or these functions throw a `MethodError`.
- If `mul_effect` throws for a pair of neighbours, so do these functions. In FermionicHilbertSpaces, a product of two kets on the same space now throws an `ArgumentError` instead of giving 0.
- Input that is already sorted pays for the check: `inner` took up to a third longer in our measurements, depending on the type and size of the kets.
- `contraction_key` is now compared at every position up to the first factor without a key. Keys must differ between subsystems, and `mul_effect` must sort kets of different subsystems into one order, since other orders count as orthogonal.

**The macros throw an `ArgumentError` for more input.**
- `@nc`: `Number` subtypes, `Any`, unions that include a number, and tuple syntax.
- `@commutative`: repeated or overlapping types, types not yet registered with `@nc`, tuple syntax, and a pair in the opposite order to an earlier call. Reordering a pair needs a Julia restart.

**`mul_effect` return values.**
- An invalid return value now throws an `ArgumentError`. This includes `Swap(missing)`, invalid values inside `AddTerms`, and atoms not registered with `@nc`.
- `mul_effect` methods defined by `@commutative` return `Swap(1)` instead of `y * x`.

### New

- `Swap` and `AddTerms` are public.
- `atom_group`, `foreach_segment`, `split_by_group` and `segment_ranges` split products into runs of factors in the same group. See [Multiple types](docs/src/multiple_types.md).
- `truncate_terms` and `opnorm_bound` drop small terms and return a bound on the error. See [Linear algebra](docs/src/linear_algebra.md).
- `add!!` accepts atoms and numbers in any position, as in `add!!(k1, k2)` and `add!!(1, a * c, 2, 3)`.

The functions in the second and third items are documented but not yet `public`.

### Fixed

- `eigsolve` with an expression as the operator could return wrong eigenvalues with autosort off: -18.0 for a Hubbard dimer with ground energy -0.83.
- `Inf * x` and `NaN * x` keep a zero scalar part zero.
- A zero product is `isequal` only to `zero(prefactor)`, matching `hash`, so `isequal(-(0.0 * a), -0.0)` is false.
- Sums with a `typemin` coefficient print correctly.
- `H^k` for a sum, `ncmap`, and `VectorInterface.add!` with `β = 1` are much faster on large sums.
