"""
    AddTerms(terms)

Return value for [`mul_effect`](@ref) that replaces `a*b` by the sum of `terms`. `terms` is an iterable such as a
vector or a tuple, and its elements can be any other allowed return value of `mul_effect` except `nothing`, including
sums and nested `AddTerms`.
"""
struct AddTerms{T}
    terms::T
end
"""
    Swap(λ)

Return value for [`mul_effect`](@ref) that replaces `a*b` by `λ*b*a`.
"""
struct Swap{T}
    λ::T
end


"""
    sort(x::Union{NCMul,NCAdd})

Sort `x` by repeatedly applying [`mul_effect`](@ref) to neighbouring factors. The result is always an `NCAdd`, also
when `x` is a product that sorts to a single term.
"""
Base.sort(x::MulAdd) = bubble_sort(x)
function bubble_sort(a::NCMul)
    sorted = _without_autosort(() -> sort_into!!(zero(_sum_type(a)), prefactor(a), copy_factors(a.factors)))
    return _with_zero_grades(sorted, () -> _grades(a))
end
function bubble_sort(x::NCAdd)
    buffer = _factor_eltype(typeof(x))[]
    sorted = _without_autosort() do
        fold_terms!!(zero(_sum_type(x)), x) do acc, t
            sort_into!!(acc, prefactor(t), append!(empty!(buffer), t.factors))
        end
    end
    return _with_zero_grades(sorted, () -> _grades(x))
end
# Sorts the product `a` into a sum, using its factors as scratch space, so they must not be shared. A product with at
# most one factor is returned as is.
function bubble_sort!(a::NCMul)
    length(a.factors) <= 1 && return a
    return _without_autosort(() -> sort_into!!(zero(_sum_type(a)), prefactor(a), factors_vector(a.factors)))
end
_without_autosort(f) = with(f, _autosort => false)

# foldl(f, itr; init=acc), where the loop continues in a new specialization whenever `f` returns an accumulator of
# another type
fold_retyping!!(f, acc, itr) = _fold_retyping!!(f, acc, itr, ())
# `state` is () to start the iteration, or (s,) to continue it from the state s
function _fold_retyping!!(f, acc::A, itr, state) where {A}
    next = iterate(itr, state...)
    while !isnothing(next)
        x, s = next
        newacc = f(acc, x)
        newacc isa A || return _fold_retyping!!(f, newacc, itr, (s,))
        acc = newacc
        next = iterate(itr, s)
    end
    return acc
end
# fold_retyping!! over terms(x)
fold_terms!!(f, acc, x::NCMul) = iszero(x) ? acc : f(acc, x)
function fold_terms!!(f, acc, x::NCAdd)
    iszero(additive_coeff(x)) || (acc = f(acc, NCMul(additive_coeff(x), _scalar_factors(x))))
    return fold_retyping!!(f, acc, NCterms(x))
end

# Adds the term `coeff * factors` to the sum `acc`, reusing its storage when the term fits, and returns the result.
# The factors are copied, so callers can pass a scratch vector. A term without factors is added to the scalar part.
function add_term!!(acc::NCAdd, coeff, factors)
    iszero(coeff) && return acc
    isempty(factors) && return set_coeff!!(acc, additive_coeff(acc) + coeff)
    K = _term_key_type(acc, factors)
    K === keytype(acc.dict) || return _add_factors!!(convert(NCAdd{valtype(acc.dict),K}, acc), coeff, factors)
    return _add_factors!!(acc, coeff, factors)
end
function _add_factors!!(acc::NCAdd{C,K}, coeff, factors) where {C,K}
    dict = add_key!!(acc.dict, convert(K, NCMul(1, collect(_factor_eltype(K), factors))), coeff)
    return dict === acc.dict ? acc : NCAdd(additive_coeff(acc), dict)
end
# The key type of `acc` with a term with `factors` added. The first term of an empty sum sets the key type, with the
# factors narrowed to a concrete type they all share.
function _term_key_type(acc::NCAdd{C,K}, factors) where {C,K}
    T = _narrow_eltype(eltype(factors), factors)
    S = _factor_eltype(K)
    (T === S || (T <: S && !isempty(acc.dict))) && return K
    termkey = NCMul{Int,T,Vector{T}}
    return isempty(acc.dict) ? termkey : promote_type(K, termkey)
end
# The type of all elements of `factors` if their element type T is abstract and they share a concrete type. Any is
# kept, since sums from `anyadd` must keep their type.
function _narrow_eltype(::Type{T}, factors) where {T}
    (isconcretetype(T) || T === Any) && return T
    U = typeof(first(factors))
    return all(f -> f isa U, factors) ? U : T
end

# Add `coeff` to the term with key `key` in `dict`, and delete the term if it becomes zero. `key` is stored as given.
# Mutates `dict` if the result fits, and otherwise returns a widened copy (e.g. true + true isa Int). A result that
# fits is converted before the zero check, since the conversion can round it to zero (e.g. Rational to Float16).
function add_key!!(dict::AbstractDict, key, coeff)
    newdict, _ = modify!!(dict, key) do val
        newval = _fit_value(valtype(dict), isnothing(val) ? coeff : something(val) + coeff)
        iszero(newval) ? nothing : newval
    end
    return newdict
end
_fit_value(::Type{V}, x) where {V} = promote_type(V, typeof(x)) == V ? convert(V, x) : x

# Sorts `coeff * factors` and adds the result to `acc`, using `factors` as scratch space. Callers disable autosort, so
# that the products that mul_effect forms are not sorted recursively.
sort_into!!(acc, coeff, factors) = iszero(coeff) ? acc : bubble!!(acc, coeff, factors, nothing, 1)

# Applies mul_effect to neighbouring factors from position `start` until none applies, and adds the sorted product
# and then the products in `pending` to `acc`. Effects that split the product (AddTerms, NCAdd) add the other products
# to `pending`. Continues in a new specialization when an effect changes the type of the coefficient or the factors.
function bubble!!(acc, coeff::C, factors::F, pending, start) where {C,F}
    i = start
    while i < length(factors)
        effect = _mul_effect(factors, factors[i], factors[i+1])
        isnothing(effect) && (i += 1; continue)
        effect isa Union{AddTerms,NCAdd} && return _split!!(acc, coeff, factors, i, effect, pending)
        newcoeff, newfactors = apply_effect!!(coeff, factors, i, effect)
        iszero(newcoeff) && return _sort_pending!!(acc, pending)
        i = max(i - 1, 1)
        if newcoeff isa C && newfactors isa F
            coeff, factors = newcoeff, newfactors
        else
            return bubble!!(acc, newcoeff, newfactors, pending, i)
        end
    end
    return _sort_pending!!(add_term!!(acc, coeff, factors), pending)
end
_sort_pending!!(acc, ::Nothing) = acc
_sort_pending!!(acc, pending) = fold_retyping!!((acc, p) -> sort_into!!(acc, prefactor(p), p.factors), acc, pending)
# Mixed-type products dispatch mul_effect at runtime anyway. The inference barrier keeps the compiled sorting loop
# independent of the mul_effect methods, so downstream methods don't invalidate the precompiled code.
_mul_effect(::Vector{Any}, a, b) = Base.inferencebarrier(mul_effect)(a, b)
_mul_effect(factors, a, b) = mul_effect(a, b)

# Continues bubble!! after a sum of effects on factors[i:i+1]: every part but the first is applied to a copy of the
# factors and added to `pending`, and the first is applied in place.
function _split!!(acc, coeff, factors, i, effect, pending)
    parts = effect.terms
    next = iterate(parts)
    isnothing(next) && return _sort_pending!!(acc, pending)
    first_part, state = next
    pending = _push_parts!!(pending, coeff, factors, i, Base.rest(parts, state))
    return _apply_part!!(acc, coeff, factors, i, first_part, pending)
end
# A sum splits into its scalar part, applied in place when nonzero, and its products
function _split!!(acc, coeff, factors, i, effect::NCAdd, pending)
    λ = additive_coeff(effect)
    iszero(λ) && return _split!!(acc, coeff, factors, i, AddTerms(NCterms(effect)), pending)
    pending = _push_parts!!(pending, coeff, factors, i, NCterms(effect))
    return _apply_part!!(acc, coeff, factors, i, λ, pending)
end
# Continues bubble!! after applying `part` to factors[i:i+1] in place
function _apply_part!!(acc, coeff, factors, i, part, pending)
    part isa Union{AddTerms,NCAdd} && return _split!!(acc, coeff, factors, i, part, pending)
    newcoeff, newfactors = apply_effect!!(coeff, factors, i, part)
    iszero(newcoeff) && return _sort_pending!!(acc, pending)
    return bubble!!(acc, newcoeff, newfactors, pending, max(i - 1, 1))
end
# Tuples are split one element at a time, so that each element is applied in a call specialized on its type
_push_parts!!(pending, coeff, factors, i, ::Tuple{}) = pending
function _push_parts!!(pending, coeff, factors, i, parts::Tuple)
    pending = _push_part!!(pending, coeff, factors, i, first(parts))
    return _push_parts!!(pending, coeff, factors, i, Base.rest(parts, 2))
end
function _push_parts!!(pending, coeff, factors, i, parts)
    for part in parts
        pending = _push_part!!(pending, coeff, factors, i, part)
    end
    return pending
end
# A sum adds its scalar part, when nonzero, and its products
function _push_parts!!(pending, coeff, factors, i, parts::NCAdd)
    λ = additive_coeff(parts)
    iszero(λ) || (pending = _push_part!!(pending, coeff, factors, i, λ))
    return _push_parts!!(pending, coeff, factors, i, NCterms(parts))
end
_push_parts!!(pending, coeff, factors, i, parts::AddTerms) = _push_parts!!(pending, coeff, factors, i, parts.terms)
# `pending` with the products from applying `part` to a copy of `factors` added
function _push_part!!(pending, coeff, factors, i, part)
    part isa Union{AddTerms,NCAdd} && return _push_parts!!(pending, coeff, factors, i, part)
    newcoeff, newfactors = apply_effect!!(coeff, copy(factors), i, part)
    iszero(newcoeff) && return pending
    product = NCMul(newcoeff, newfactors)
    return isnothing(pending) ? [product] : push!!(pending, product)
end

# The coefficient and factors after replacing factors[i:i+1] according to `effect`
function apply_effect!!(coeff, factors, i, effect)
    λ, replacement = _effect_parts(effect, factors[i], factors[i+1])
    return coeff * λ, replace_pair!!(factors, i, replacement)
end

# The scale factor and the replacement factors of an effect on the pair a*b
_effect_parts(λ::Number, a, b) = (λ, ())
_effect_parts(swap::Swap, a, b) = (swap.λ, (b, a))
_effect_parts(term::NCMul, a, b) = (prefactor(term), term.factors)
# NCMul(atom) is defined by @nc, so values of other types throw a MethodError
function _effect_parts(atom, a, b)
    _check_atom_effect(atom)
    return _is_atom(typeof(atom)) ? (1, (atom,)) : _effect_parts(NCMul(atom)::NCMul, a, b)
end
_check_atom_effect(term) = nothing
# `nothing` inside AddTerms would keep a*b as a term, which gives the same effect again and never terminates
_check_atom_effect(::Nothing) = throw(ArgumentError("unsupported mul_effect value: `nothing` can only be returned as the whole result of mul_effect, not inside AddTerms."))

# `factors` with factors[i:i+1] replaced by the elements of `replacement`. Mutates `factors` if the elements fit, and
# otherwise returns a new vector with a wider element type. Replacements of at most two factors are written in place.
function replace_pair!!(factors, i, replacement)
    T = eltype(factors)
    S = promote_type(eltype(replacement), T)
    S <: T || return append!(S[], view(factors, 1:i-1), replacement, view(factors, i+2:lastindex(factors)))
    m = length(replacement)
    m > 2 && (splice!(factors, i:i+1, replacement); return factors)
    for (k, x) in enumerate(replacement)
        factors[i+k-1] = x
    end
    m < 2 && deleteat!(factors, i+m:i+1)
    return factors
end

"""
    mul_effect(a, b)

When sorting a non-commutative expression, `mul_effect(a, b)` gives the replacement for `a*b`. 
It must have a branch that returns `nothing` in order for the sorting to terminate. 

This package defines no methods for `mul_effect`, so define them for your types.

The return value can be `nothing`, a number, an object, a sum or product of these, or one of the helper types `Swap(λ)`
(replace `a*b` by `λ*b*a`) and `AddTerms(terms)` (replace `a*b` by a sum of `terms`). `Swap` and `AddTerms` are
available as `NonCommutativeProducts.Swap` and `NonCommutativeProducts.AddTerms`.

# Example
```julia
# Anticommuting Grassmann numbers
struct θ id::Int end
NonCommutativeProducts.@nc θ
NonCommutativeProducts.mul_effect(a::θ, b::θ) = a.id < b.id ? nothing : -b * a * (a.id==b.id ? 0 : 1)
NonCommutativeProducts.enable_autosort!()
θ(2)*θ(1) # = -θ(1)*θ(2)
θ(1)*θ(1) # = 0
```
"""
function mul_effect end
