# Truncation with a certified error bound. For x = Σᵢ cᵢPᵢ, dropping a set D of terms changes x by an operator of norm
# at most Σ_{i∈D} |cᵢ| b(Pᵢ), by the triangle inequality, where b(P) bounds the operator norm of the product P. The
# bounds are computed in floating point, rounded up, so that they are never smaller than the exact values.

"""
    opnorm_bound(x)

An upper bound on the operator norm of `x`, or `nothing` if no bound is known. Define it for your atom types, with a
bound that holds in every representation you use, so that [`truncate_terms`](@ref) can drop terms containing them.
Products and sums get their bound from their atoms and coefficients.

# Example
```julia
NonCommutativeProducts.opnorm_bound(::Fermion) = 1
```
"""
opnorm_bound(x) = nothing
opnorm_bound(x::Number) = _bound_value(coefficient_bound(x))
opnorm_bound(x::NCMul) = _sum_weights(_term_weight(fs, c) for (fs, c) in termpairs(x))
function opnorm_bound(x::NCAdd)
    # the scalar part separately, since iterating over termpairs(x) isn't inferred
    coeff = additive_coeff(x)
    scalar_weight = iszero(coeff) ? 0.0 : _term_weight(_scalar_factors(x), coeff)
    weight = _sum_weights(_term_weight(k.factors, c) for (k, c) in x.dict)
    (isnothing(scalar_weight) || isnothing(weight)) && return nothing
    return _add_up(scalar_weight, weight)
end

"""
    opnorm_bound_factors(factors)

An upper bound on the operator norm of the product of the atoms in the vector `factors`, or `nothing` if no bound is
known. The default multiplies their [`opnorm_bound`](@ref)s. Define a method to give a tighter bound.

# Example
A product of fermion operators has norm at most 1, however long it is:
```julia
NonCommutativeProducts.opnorm_bound_factors(::AbstractVector{<:Fermion}) = 1
```
"""
function opnorm_bound_factors(factors)
    isconcretetype(eltype(factors)) || return _opnorm_bound_runs(factors)
    b = 1.0
    for f in factors
        bf = _bound_value(opnorm_bound(f))
        isnothing(bf) && return nothing
        b = _mul_up(b, bf)
    end
    return b
end
# The product of the bounds of maximal runs of factors of one type, a bound by submultiplicativity
function _opnorm_bound_runs(factors)
    fs = collect(factors)
    b = 1.0
    i = firstindex(fs)
    while i <= lastindex(fs)
        S = typeof(fs[i])
        j = i
        while j < lastindex(fs) && typeof(fs[j+1]) === S
            j += 1
        end
        br = _bound_value(opnorm_bound_factors(convert(Vector{S}, fs[i:j])))
        isnothing(br) && return nothing
        b = _mul_up(b, br)
        i = j + 1
    end
    return b
end

"""
    coefficient_bound(c)

An upper bound on the absolute value of the coefficient `c`, or `nothing` if it is unknown, as for a symbolic
coefficient. Define it for your own coefficient types if `abs` doesn't give a bound.
"""
coefficient_bound(c) = nothing
coefficient_bound(c::Number) = _coefficient_bound(c, abs(c))
# widened, since abs(typemin(Int)) < 0
coefficient_bound(c::Integer) = abs(widen(c))
coefficient_bound(c::Rational) = abs(widen(c))
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
    isnothing(f) && return nothing
    isnan(f) && return oftype(f, Inf)
    f >= 0 || throw(ArgumentError("bounds on norms must be nonnegative, got $b"))
    return f
end
_float_up(x::AbstractFloat) = x
function _float_up(x::Real)
    f = float(x)
    f isa AbstractFloat || return nothing # e.g. a dual number
    # the comparison of a float with an integer, rational or irrational number is exact
    return f < x ? nextfloat(f) : f
end
# x as a T, rounded down
function _float_down(::Type{T}, x::Real) where T
    f = convert(T, x)
    return f > x ? prevfloat(f) : f
end

# a * b and a + b for nonnegative floats, rounded up, and a * b rounded down. fma and TwoSum detect inexact results,
# but the fma residual of a product below floatmin may underflow, so such products are rounded outward.
function _mul_up(a::AbstractFloat, b::AbstractFloat)
    a, b = promote(a, b)
    p = a * b
    isnan(p) && return oftype(p, Inf) # 0 * Inf
    (iszero(a) || iszero(b) || !isfinite(p)) && return p
    p < floatmin(p) && return nextfloat(p)
    return fma(a, b, -p) > 0 ? nextfloat(p) : p
end
function _mul_down(a::AbstractFloat, b::AbstractFloat)
    a, b = promote(a, b)
    (iszero(a) || iszero(b)) && return zero(a) # also for 0 * Inf
    p = a * b
    isfinite(p) || return p
    p < floatmin(p) && return max(prevfloat(p), zero(p))
    return fma(a, b, -p) < 0 ? prevfloat(p) : p
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
    total = 0.0 # promoted by _add_up if a weight is wider
    for w in weights
        isnothing(w) && return nothing
        total = _add_up(total, w)
    end
    return total
end

"""
    truncate_terms(x; atol=0, rtol=0)

Drop small terms of `x`. Returns `(; expr, err)`, the truncated expression and a bound on the error:
```
opnorm(x - expr) ≤ err ≤ max(atol, rtol * B)
```
in every representation where [`opnorm_bound`](@ref) holds. A term `c * P` has the weight `|c|` times the bound on `P`,
and `B` is the total weight of the terms that can be dropped. Terms are dropped from the smallest weight up, while
their total fits the budget. Terms with an unknown or infinite weight are always kept, such as those containing an atom
without a bound.

# Example
```julia
NonCommutativeProducts.opnorm_bound(::Fermion) = 1
expr, err = truncate_terms(c[1]' * c[2] + 0.01 * c[1] + 0.001 * c[2]; atol=0.005)
# expr == c[1]' * c[2] + 0.01 * c[1] and err == 0.001
```
"""
function truncate_terms(x::NCAdd{C,K}; atol::Real=0, rtol::Real=0) where {C,K}
    coeff = additive_coeff(x)
    ps = collect(x.dict)
    # a zero scalar part has weight nothing, which keeps it
    scalar_weight = iszero(coeff) ? nothing : _term_weight(_scalar_factors(x), coeff)
    weights = [scalar_weight; [_term_weight(k.factors, c) for (k, c) in ps]]
    dropped, err = _greedy_drop(weights, atol, rtol)
    newcoeff = dropped[1] ? zero(coeff) : coeff
    dict = Dict{K,C}(ps[i] for i in eachindex(ps) if !dropped[i+1])
    result = NCAdd{C,K}(newcoeff, dict; filter_zeros=false, filter_scalars=false)
    return (; expr=_with_zero_grades(result, () -> _grades(x)), err)
end
function truncate_terms(x::NCMul; atol::Real=0, rtol::Real=0)
    weight = iszero(x) ? nothing : _term_weight(x.factors, prefactor(x))
    dropped, err = _greedy_drop([weight], atol, rtol)
    # a zero product keeps its factors, for its grade
    return (; expr=only(dropped) ? NCMul(_zero_coefficient(prefactor(x)), x.factors) : x, err)
end
# zero(One()) is a Zero, which isn't a Number
_zero_coefficient(c) = zero(c)
_zero_coefficient(::One) = false
function truncate_terms(x::Number; atol::Real=0, rtol::Real=0)
    weight = iszero(x) ? nothing : _bound_value(coefficient_bound(x))
    dropped, err = _greedy_drop([weight], atol, rtol)
    return (; expr=only(dropped) ? zero(x) : x, err)
end
function truncate_terms(x; kwargs...)
    _is_atom(typeof(x)) || throw(ArgumentError("truncate_terms needs a sum, a product, a number or an atom, got $x"))
    return truncate_terms(NCMul(x); kwargs...)
end

# max(atol, rtol * scale), rounded down. rtol * 0 is 0 also for rtol = Inf.
function _budget(atol, rtol, scale::T) where T
    for (tol, name) in ((atol, :atol), (rtol, :rtol))
        (tol >= 0) === true || throw(ArgumentError("`$name` must be a nonnegative number, got $tol"))
    end
    return max(_float_down(T, atol), _mul_down(_float_down(T, rtol), scale))
end

# The float type of the bounds: Float64, or wider if a weight is wider
function _weight_type(weights::AbstractVector{W}) where W
    T = typeintersect(W, AbstractFloat)
    isconcretetype(T) && return promote_type(Float64, T)
    return mapreduce(typeof, promote_type, (w for w in weights if w isa AbstractFloat); init=Float64)
end

# Which of the terms with `weights` to drop, and the bound on the dropped part. Finite weights are dropped in
# increasing order, a group of equal weights at a time, while their sum is within the budget.
function _greedy_drop(weights, atol, rtol)
    T = _weight_type(weights)
    droppable = [i for i in eachindex(weights) if !isnothing(weights[i]) && isfinite(weights[i])]
    scale = foldl((s, i) -> _add_up(s, T(weights[i])), droppable; init=zero(T))
    budget = _budget(atol, rtol, scale)
    sort!(droppable; by=i -> weights[i])
    dropped = falses(length(weights))
    err = zero(T)
    i = 1
    while i <= length(droppable)
        w = T(weights[droppable[i]])
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
