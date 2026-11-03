# Truncation with an error bound. For x = Σᵢ cᵢPᵢ, dropping a set D of terms changes x by an operator of norm at most
# Σ_{i∈D} |cᵢ| b(Pᵢ), by the triangle inequality, where b(P) bounds the operator norm of the product P.

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
opnorm_bound(x::Number) = _term_weight((), x)
opnorm_bound(x::NCMul) = _sum_weights(_term_weight(fs, c) for (fs, c) in termpairs(x))
function opnorm_bound(x::NCAdd)
    # the scalar part separately, since iterating over termpairs(x) isn't inferred
    coeff = additive_coeff(x)
    scalar_weight = iszero(coeff) ? 0.0 : _term_weight(_scalar_factors(x), coeff)
    weight = _sum_weights(_term_weight(k.factors, c) for (k, c) in x.dict)
    (isnothing(scalar_weight) || isnothing(weight)) && return nothing
    return scalar_weight + weight
end

# |c| as a float, or nothing if it is unknown, as for a symbolic coefficient. The SymbolicUtils extension adds a method
# for _unwrap_symbolic, so that wrapped symbolic constants such as Symbolics.Num get the bound of their value.
_coefficient_bound(c) = nothing
_coefficient_bound(::VectorInterface.Zero) = 0.0
_coefficient_bound(::One) = 1.0
# float first, since abs(typemin(Int)) < 0
_coefficient_bound(c::Union{Integer,Rational}) = abs(float(c))
function _coefficient_bound(c::Number)
    b = abs(c)
    b isa Real && (b >= 0) isa Bool && return float(b)
    applicable(_unwrap_symbolic, c) || return nothing
    u = _unwrap_symbolic(c)
    return u === c ? nothing : _coefficient_bound(u)
end
function _coefficient_bound(c::Complex)
    a = _coefficient_bound(real(c))
    b = _coefficient_bound(imag(c))
    return isnothing(a) || isnothing(b) ? nothing : hypot(a, b)
end
function _unwrap_symbolic end

# The bound |c| b(factors) of a term, Inf for NaN, or nothing if it is unknown
function _term_weight(factors, c)
    w = _coefficient_bound(c)
    isnothing(w) && return nothing
    bounds = [opnorm_bound(f) for f in factors]
    any(isnothing, bounds) && return nothing
    any(<(0), bounds) && throw(ArgumentError("bounds on norms must be nonnegative, got $bounds"))
    # sorted, so that a product and its adjoint get exactly the same weight
    for b in sort!(bounds; rev=true)
        w *= b
    end
    # the weight of a nonzero term must stay positive when it underflows, since dropping the term is an error
    iszero(w) && !iszero(c) && !any(iszero, bounds) && return nextfloat(zero(w))
    return isnan(w) ? oftype(w, Inf) : w
end
function _sum_weights(weights)
    total = 0.0
    for w in weights
        isnothing(w) && return nothing
        total += w
    end
    return total
end

"""
    truncate_terms(x; atol=0, rtol=0)

Drop small terms of `x`. Returns `(; expr, err)`, the truncated expression and a bound on the error. Up to rounding,
```
opnorm(x - expr) ≤ err ≤ max(atol, rtol * B)
```
in every representation where [`opnorm_bound`](@ref) holds. A term `c * P` has the weight `|c|` times the bound on `P`,
and `B` is the total weight of the terms that can be dropped. Terms are dropped from the smallest weight up, while
their total fits the budget, and terms of equal weight are dropped together or not at all. Terms with an unknown or
infinite weight are always kept, such as those containing an atom without a bound.

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
    weight = iszero(x) ? nothing : _term_weight((), x)
    dropped, err = _greedy_drop([weight], atol, rtol)
    return (; expr=only(dropped) ? zero(x) : x, err)
end
function truncate_terms(x; kwargs...)
    _is_atom(typeof(x)) || throw(ArgumentError("truncate_terms needs a sum, a product, a number or an atom, got $x"))
    return truncate_terms(NCMul(x); kwargs...)
end

# Which of the terms with `weights` to drop, and the bound on the dropped part. Finite weights are dropped in
# increasing order, a group of equal weights at a time, while their sum is within the budget.
function _greedy_drop(weights, atol, rtol)
    for (tol, name) in ((atol, :atol), (rtol, :rtol))
        (tol >= 0) === true || throw(ArgumentError("`$name` must be a nonnegative number, got $tol"))
    end
    droppable = [i for i in eachindex(weights) if !isnothing(weights[i]) && isfinite(weights[i])]
    sort!(droppable; by=i -> weights[i])
    # summed in the order of the loop below, so that rtol = 1 drops everything, and scaled term by term, since the total
    # can overflow
    budget = max(atol, foldl((s, i) -> s + rtol * weights[i], droppable; init=0.0))
    dropped = falses(length(weights))
    err = 0.0
    i = 1
    while i <= length(droppable)
        w = weights[droppable[i]]
        j = i
        total = err
        while j <= length(droppable) && weights[droppable[j]] == w
            total += w
            j += 1
        end
        total <= budget || break
        dropped[droppable[i:j-1]] .= true
        err = total
        i = j
    end
    return dropped, err
end
