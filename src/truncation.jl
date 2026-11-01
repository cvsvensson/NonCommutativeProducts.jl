# Truncation with a certified error bound. For x = Σᵢ cᵢPᵢ, dropping a set D of terms changes x by an operator of norm
# at most Σ_{i∈D} |cᵢ| b(Pᵢ), by the triangle inequality, where b(P) bounds the operator norm of the product P. The
# bounds are computed in floating point, rounded up, so that they are never smaller than the exact values.

"""
    opnorm_bound(x)

An upper bound on the operator norm of `x`, or `nothing` if no bound is known.

Define methods for your atom types. A method must return a nonnegative real number that bounds the operator norm of the
atom in every representation of the atoms by bounded operators that you will use, or `nothing`, the default, when there
is no such bound. Terms with an atom without a bound are never dropped by [`truncate_terms`](@ref). For example:
- a fermion `c` or `c'` with `c' * c + c * c' == 1`: `1`
- a Majorana `γ` with `γ^2 == 1`, a Pauli matrix: `1`
- a component of a spin `s`, such as `Sz`: `s`
- a ket-bra `|a⟩⟨b|` of normalized states: `1`
- a boson `a`: no method, since `a` is unbounded. Return a bound only for an explicit cutoff, e.g. `sqrt(nmax)`.

For products and sums, `opnorm_bound` is computed from the bounds of the atoms: a product `c * a₁ * … * aₙ` has the
bound [`coefficient_bound`](@ref)`(c) *` [`opnorm_bound_factors`](@ref)`([a₁, …, aₙ])`, and a sum has the sum of the
bounds of its terms, with `1` for the factors of its scalar part. The result is `nothing` if the bound of a term is
unknown, and `Inf` if a coefficient is `NaN` or infinite. These are floating-point numbers, rounded up.

# Example
```julia
NonCommutativeProducts.opnorm_bound(::Fermion) = 1
opnorm_bound(0.5 * c[1]' * c[2] - 2 * c[1] + 1) == 3.5
opnorm_bound(c[1] + 0.5 * b) === nothing # if b has no bound
```
"""
opnorm_bound(x) = nothing
opnorm_bound(x::Number) = _bound_value(coefficient_bound(x))
opnorm_bound(x::NCMul) = _sum_weights(_term_weight(fs, c) for (fs, c) in termpairs(x))
opnorm_bound(x::NCAdd) = _sum_weights(_term_weight(fs, c) for (fs, c) in termpairs(x))

"""
    opnorm_bound_factors(factors)

An upper bound on the operator norm of the product of `factors`, which is an `AbstractVector` of atoms, or `nothing`
if no bound is known. The default is the product of [`opnorm_bound`](@ref) of the factors, by
submultiplicativity of the operator norm, and `nothing` if one of them is `nothing`. The product of no factors has the
bound `1`.

Define methods to give tighter bounds for whole products, valid in the same representations as `opnorm_bound`. For
example, a product of fermion operators is zero or a partial isometry, so its norm is at most `1` however long it is.
A method can dispatch on the element type of `factors`, which may be abstract, e.g. `Any`, for products of atoms of
several types.

# Example
```julia
NonCommutativeProducts.opnorm_bound_factors(fs::AbstractVector{<:Fermion}) = 1
```
"""
function opnorm_bound_factors(factors)
    b = 1.0
    for f in factors
        bf = _bound_value(opnorm_bound(f))
        isnothing(bf) && return nothing
        b = _mul_up(b, bf)
    end
    return b
end

"""
    coefficient_bound(c)

An upper bound on the absolute value of the coefficient `c`, or `nothing` if it is unknown. The default is `abs(c)` for
numbers, rounded up for complex numbers, and `nothing` for other values. A coefficient is also treated as unknown if
its bound is not a real number that can be compared with `0`, as for symbolic coefficients. With SymbolicUtils loaded,
symbolic constants have the bound of their value.

Terms with coefficients without a bound are never dropped by [`truncate_terms`](@ref). Define methods for your own
coefficient types if `abs` doesn't give a bound.
"""
coefficient_bound(c) = nothing
coefficient_bound(c::Number) = _coefficient_bound(c, abs(c))
# widened, since abs(typemin(Int)) < 0
coefficient_bound(c::Integer) = abs(widen(c))
coefficient_bound(c::Bool) = c
coefficient_bound(::VectorInterface.Zero) = false
coefficient_bound(::One) = true
function coefficient_bound(c::Complex)
    a = _bound_value(coefficient_bound(real(c)))
    b = _bound_value(coefficient_bound(imag(c)))
    (isnothing(a) || isnothing(b)) && return nothing
    iszero(b) && return a
    iszero(a) && return b
    # hypot is accurate to about one ulp
    h = hypot(a, b)
    return isfinite(h) ? nextfloat(h, 2) : h
end
# A bound that doesn't compare with 0 is symbolic. The SymbolicUtils extension adds a method for _unwrap_symbolic, so
# that wrapped symbolic constants such as Symbolics.Num get the bound of their value.
function _coefficient_bound(c, b)
    _is_comparable(b) && return b
    applicable(_unwrap_symbolic, c) || return nothing
    u = _unwrap_symbolic(c)
    return u === c ? nothing : coefficient_bound(u)
end
function _unwrap_symbolic end
_is_comparable(b) = b isa Real && (b >= 0) isa Bool

# `b` as a floating-point number rounded up, Inf for NaN, or nothing if it is unknown
_bound_value(::Nothing) = nothing
function _bound_value(b)
    _is_comparable(b) || return nothing
    f = _float_up(b)
    isnan(f) && return oftype(f, Inf)
    f >= 0 || throw(ArgumentError("bounds on norms must be nonnegative, got $b"))
    return f
end
_float_up(x::AbstractFloat) = x
function _float_up(x::Real)
    f = float(x)
    # the comparison of a float with an integer, rational or irrational number is exact
    return f < x ? nextfloat(f) : f
end

# a * b and a + b for nonnegative floats, rounded up. Error-free transformations detect inexact results.
function _mul_up(a::AbstractFloat, b::AbstractFloat)
    a, b = promote(a, b)
    p = a * b
    isnan(p) && return oftype(p, Inf) # 0 * Inf
    return isfinite(p) && fma(a, b, -p) > 0 ? nextfloat(p) : p
end
function _add_up(a::AbstractFloat, b::AbstractFloat)
    a, b = promote(a, b)
    s = a + b
    isfinite(s) || return s
    bb = s - a
    err = (a - (s - bb)) + (b - bb)
    return err > 0 ? nextfloat(s) : s
end

# The bound |c| b(factors) of a term, or nothing if it is unknown
function _term_weight(factors, c)
    cb = _bound_value(coefficient_bound(c))
    isnothing(cb) && return nothing
    fb = _bound_value(opnorm_bound_factors(factors_vector(factors)))
    isnothing(fb) && return nothing
    return _mul_up(cb, fb)
end
function _sum_weights(weights)
    total = 0.0
    for w in weights
        isnothing(w) && return nothing
        total = _add_up(total, w)
    end
    return total
end

"""
    truncate_terms(x; atol=0, rtol=0)

Drop terms of `x` whose total contribution has a certified bound, and return the named tuple `(; expr, error)` of the
truncated expression `expr` and a bound `error` with
```
opnorm(x - expr) ≤ error ≤ max(atol, rtol * B)
```
in every representation of the atoms by bounded operators where [`opnorm_bound`](@ref) holds. `x` can be a sum, a
product, a number or an atom.

For `x = Σᵢ cᵢ Pᵢ`, the weight of the term `cᵢ Pᵢ` is `|cᵢ| b(Pᵢ)`, where `|cᵢ|` is [`coefficient_bound`](@ref)`(cᵢ)`
and `b(Pᵢ)` is [`opnorm_bound_factors`](@ref) of its factors, or `1` for the scalar part. By the triangle inequality,
dropping a set of terms changes `x` by at most the sum of their weights, which is `error`. The terms are dropped in
order of increasing weight while the sum stays within the budget `max(atol, rtol * B)`, where `B` is the sum of the
weights of the terms that can be dropped, which is `opnorm_bound(x)` if all terms have finite weights. Dropping the
smallest weights first drops as many terms as possible within the budget, except that terms with the same weight are
dropped together or not at all, so that the result doesn't depend on the order of the terms.

Terms with an unknown weight are never dropped: those with an atom without a bound, such as a boson, and those with a
coefficient without a bound, such as a symbolic one. Terms with a `NaN` or infinite coefficient are never dropped
either. With `atol = rtol = 0`, only terms of weight exactly zero are dropped.

The weights and `error` are computed in floating point and rounded up, so `error` is never smaller than the exact sum
of the dropped weights. The bound holds whether or not the terms of `x` are sorted and combined, but combining like
terms first, e.g. with autosort enabled, makes the best use of the budget, since cancelling terms then don't use it.
The bound is often pessimistic, by about the square root of the number of dropped terms.

# Example
```julia
NonCommutativeProducts.opnorm_bound(::Fermion) = 1
expr, error = truncate_terms(c[1]' * c[2] + 0.01 * c[1] + 0.001 * c[2]; atol=0.005)
# expr == c[1]' * c[2] + 0.01 * c[1] and error == 0.001
```
"""
function truncate_terms(x::NCAdd{C,K}; atol::Real=0, rtol::Real=0) where {C,K}
    _check_tolerance(atol, :atol)
    _check_tolerance(rtol, :rtol)
    coeff = additive_coeff(x)
    ps = collect(x.dict)
    # index 0 is the scalar part, with weight nothing if it is zero, which keeps it
    weights = [i == 0 ? (iszero(coeff) ? nothing : _term_weight(_scalar_factors(x), coeff)) :
               _term_weight(ps[i][1].factors, ps[i][2]) for i in 0:length(ps)]
    dropped, err = _greedy_drop(weights, atol, rtol)
    newcoeff = dropped[1] ? zero(coeff) : coeff
    dict = Dict{K,C}(ps[i] for i in eachindex(ps) if !dropped[i+1])
    result = NCAdd{C,K}(newcoeff, dict; filter_zeros=false, filter_scalars=false)
    return (; expr=_with_zero_grades(result, () -> _grades(x)), error=err)
end
function truncate_terms(x::NCMul; atol::Real=0, rtol::Real=0)
    _check_tolerance(atol, :atol)
    _check_tolerance(rtol, :rtol)
    iszero(x) && return (; expr=x, error=0.0)
    dropped, err = _greedy_drop([_term_weight(x.factors, prefactor(x))], atol, rtol)
    # a zero product keeps its factors, for its grade
    return (; expr=only(dropped) ? NCMul(zero(prefactor(x)), x.factors) : x, error=err)
end
function truncate_terms(x::Number; atol::Real=0, rtol::Real=0)
    _check_tolerance(atol, :atol)
    _check_tolerance(rtol, :rtol)
    iszero(x) && return (; expr=x, error=0.0)
    dropped, err = _greedy_drop([_bound_value(coefficient_bound(x))], atol, rtol)
    return (; expr=only(dropped) ? zero(x) : x, error=err)
end
function truncate_terms(x; kwargs...)
    _is_atom(typeof(x)) || throw(ArgumentError("truncate_terms needs a sum, a product, a number or an atom, got $x"))
    return truncate_terms(NCMul(x); kwargs...)
end

function _check_tolerance(tol, name)
    (tol >= 0) === true || throw(ArgumentError("`$name` must be a nonnegative number, got $tol"))
    return nothing
end

# Which of the terms with `weights` to drop, and the bound on the dropped part. Finite weights are dropped in
# increasing order, a group of equal weights at a time, while their sum is within the budget.
function _greedy_drop(weights, atol, rtol)
    droppable = [i for i in eachindex(weights) if !isnothing(weights[i]) && isfinite(weights[i])]
    budget = atol
    if rtol > 0
        scale = foldl((s, i) -> _add_up(s, weights[i]), droppable; init=0.0)
        budget = max(budget, rtol * scale)
    end
    sort!(droppable; by=i -> weights[i])
    dropped = falses(length(weights))
    err = 0.0
    i = 1
    while i <= length(droppable)
        w = weights[droppable[i]]
        j = i
        total = err
        while j <= length(droppable) && weights[droppable[j]] == w
            total = _add_up(total, w)
            j += 1
        end
        total <= budget || break
        dropped[droppable[i:j-1]] .= true
        err = total
        i = j
    end
    return dropped, err
end
