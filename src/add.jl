function filter_ncadd_dict!(d::AbstractDict{K,V}; filter_zeros=true, filter_scalars=true) where {K<:NCMul,V}
    coeff = zero(V)
    !filter_zeros && !filter_scalars && return d, coeff
    for (k, v) in d
        if filter_zeros && iszero(v)
            delete!(d, k)
            continue
        end
        if filter_scalars && isscalar(k)
            delete!(d, k)
            coeff += prefactor(k) * v
        end
    end
    return d, coeff
end

filter_zeros!(d::AbstractDict) = filter!(kv -> !iszero(last(kv)), d)

# The terms are a Dict, since equality and key lookup assume keys are compared by value.
# `zero_grades` holds the grades of a zero sum, or `nothing` (see grading.jl). Equality and hashing ignore it.
"""
    NCAdd

A sum of products (`NCMul`) with a scalar term. Create sums by adding atoms and products, e.g. `a*b + 2*c + 1`. Use
`terms(x)` to iterate over the terms of a sum.

The constructors of `NCAdd` are internal. They may take ownership of the dictionary passed to them and mutate it.
"""
mutable struct NCAdd{C,K}
    coeff::C
    const dict::Dict{K,C}
    zero_grades::Union{Nothing,Set{Any}}
    # throws an InexactError if the coefficients don't fit in C
    function NCAdd{C,K}(coeff, dict::AbstractDict; zero_grades=nothing, kwargs...) where {C,K}
        newdict = dict isa Dict{K,C} ? dict : Dict{K,C}(dict)
        _, addcoeff = filter_ncadd_dict!(newdict; kwargs...)
        new{C,K}(convert(C, coeff + addcoeff), newdict, _as_grades(zero_grades))
    end
end
function NCAdd(coeff, dict::AbstractDict; zero_grades=nothing, kwargs...)
    _, addcoeff = filter_ncadd_dict!(dict; kwargs...)
    newcoeff = coeff + addcoeff
    T = promote_type(typeof(newcoeff), valtype(dict))
    return NCAdd{T,keytype(dict)}(newcoeff, dict; zero_grades, filter_zeros=false, filter_scalars=false)
end
NCAdd{C,K}(ncadd::NCAdd{C,K}) where {C,K} = ncadd
NCAdd(ncmul::NCMul{C}) where C = _with_zero_grades(NCAdd(zero(C), to_add_dict(ncmul)), () -> _grades(ncmul))

additive_coeff(x::NCAdd) = x.coeff
function set_coeff!(a::NCAdd, x::Number)
    a.coeff = x
    return a
end
function set_coeff!!(a::NCAdd{C}, x::Number) where {C}
    # use set_coeff! only where it can't fail
    promote_type(typeof(x), C) <: C ? set_coeff!(a, x) : NCAdd(x, a.dict)
end
# Every key type promotes to keys with Vector{Any} factors, so solvers such as KrylovKit can use this type throughout.
anyadd(x::Union{NCMul{C},NCAdd{C}}) where {C} = convert(NCAdd{C,NCMul{Int,Any,Vector{Any}}}, x)

const MulAdd = Union{NCMul,NCAdd}
Base.iszero(x::NCAdd) = iszero(additive_coeff(x)) && all(iszero, values(x.dict))
Base.:(==)(a::NCAdd, b::NCAdd) = additive_coeff(a) == additive_coeff(b) && a.dict == b.dict
Base.:(==)(a::NCAdd, b::Number) = additive_coeff(a) == b && isempty(a.dict)
Base.:(==)(a::Number, b::NCAdd) = a == additive_coeff(b) && isempty(b.dict)
function Base.hash(a::NCAdd, h::UInt)
    # if it's only a number, hash should equal hash of number
    if isempty(a.dict)
        return hash(additive_coeff(a), h)
    end
    # if coeff is zero and there is only one term, it should hash equals to the corresponding NCMul
    if iszero(additive_coeff(a)) && length(a.dict) == 1
        ncmul, coeff = only(a.dict)
        return hash(NCMul(coeff, ncmul.factors), h)
    end
    return hash(additive_coeff(a), hash(a.dict, h))
end

isscalar(x::NCAdd) = length(x.dict) == 0 || all(isscalar, keys(x.dict)) || all(iszero, values(x.dict))
scalar(x::NCAdd) = isscalar(x) ? additive_coeff(x) : throw(ArgumentError("NCAdd is not a scalar"))
Base.copy(x::NCAdd) = NCAdd(copy(additive_coeff(x)), copy(x.dict); zero_grades=x.zero_grades)

function print_coeff(io, coeff)
    if isreal(coeff)
        print(io, real(coeff), "I")
    else
        print(io, "(", coeff, ")", "I")
    end
end
# Terms are sorted by their signed printed form, so the output doesn't depend on dictionary order.
function _term_parts(io::IO, x::NCAdd)
    parts = Tuple{Bool,String}[]
    for (k, v) in pairs(x.dict)
        r = isreal(v) ? real(v) : v
        isneg = isreal(v) ? r < 0 : nothing
        if isneg isa Bool
            neg = isneg
            text = sprint(print, neg ? -r * k : r * k; context=io)
        else
            neg = false
            text = sprint(print, "(", v, ")*", k; context=io)
        end
        push!(parts, (neg, text))
    end
    sort!(parts; by=((neg, text),) -> neg ? "-" * text : text)
    return parts
end

function _show_terms(io::IO, x::NCAdd, parts, compact::Bool, ellipsis::Bool)
    print_one = !iszero(additive_coeff(x)) || length(x.dict) == 0
    print_sign(s) = compact ? print(io, s) : print(io, " ", s, " ")
    print_one && print_coeff(io, additive_coeff(x))
    for (n, (neg, text)) in enumerate(parts)
        if neg
            n > 1 || print_one ? print_sign("-") : print(io, "-")
        elseif n > 1 || print_one
            print_sign("+")
        end
        print(io, text)
    end
    ellipsis && print(io, " + ...")
    return nothing
end

Base.show(io::IO, x::NCAdd) = _show_terms(io, x, _term_parts(io, x), false, false)

function Base.show(io::IO, ::MIME"text/plain", x::NCAdd; max_terms=3)
    parts = _term_parts(io, x)
    if get(io, :limit, false) && length(parts) > max_terms
        println(io, "Sum with ", length(parts) + !iszero(additive_coeff(x)), " terms: ")
        _show_terms(io, x, parts[1:max_terms], true, true)
    else
        _show_terms(io, x, parts, false, false)
    end
end

Base.:+(a::Number, b::NCAdd) = _with_zero_grades(NCAdd(a + additive_coeff(b), copy(b.dict)), () -> _union_grades(_grades(a), _grades(b)))
Base.:+(a::UniformScaling, b::NCAdd) = a.λ + b
Base.:+(a::NCAdd, b::B) where B<:Union{Number,UniformScaling} = b + a
Base.:+(a::NCAdd, b::B) where B<:NCMul = b + a
Base.:/(a::NCMul, b::Number) = NCMul(prefactor(a) / b, a.factors)
Base.://(a::NCMul, b::Number) = NCMul(prefactor(a) // b, a.factors)
Base.:/(a::NCAdd, b::Number) = _divide_coefficients(/, a, b)
Base.://(a::NCAdd, b::Number) = _divide_coefficients(//, a, b)
function _divide_coefficients(op, a::NCAdd{C,K}, b::Number) where {C,K}
    coeff = op(additive_coeff(a), b)
    dict = copy_dict(a, K, typeof(coeff))
    map!(v -> op(v, b), values(dict))
    return _with_zero_grades(NCAdd(coeff, dict), () -> _grades(a))
end
Base.:-(a::Union{Number,UniformScaling}, b::MulAdd) = a + (-b)
Base.:-(a::MulAdd, b::Union{Number,MulAdd,UniformScaling}) = a + (-b)
Base.:-(a::NCAdd) = (-1) * a
function Base.:+(a::NCAdd{C1,K1}, b::NCAdd{C2,K2}) where {C1,K1,C2,K2}
    C = promote_type(C1, C2)
    dict = copy_dict(a, promote_type(K1, K2), C)
    for (k, v) in b.dict
        dict = add_key!!(dict, k, v)
    end
    # The terms have factors and add_key!! drops zeros, so only converting the coefficients of a can leave zero terms
    # (e.g. Rational to Float16).
    sum = NCAdd(additive_coeff(a) + additive_coeff(b), dict; filter_zeros=C1 !== C, filter_scalars=false)
    return _with_zero_grades(sum, () -> _union_grades(_grades(a), _grades(b)))
end
"""
    add!!(a, b)
    add!!(a, b, α)
    add!!(a, b, α, β)

Compute `a + b`, or `β * a + α * b`, reusing the storage of `a` when possible. Use the returned value, since `a` may
be left unchanged or modified. This function belongs to this package and is not `VectorInterface.add!!`. It differs
in having the 2-argument form, and in the weights `α` and `β` defaulting to `One()`, so that `add!!(a, b)` is `a + b`.
Methods for `VectorInterface.add!!` call it.
"""
add!!(a::NCMul, b::MulAdd, α::Number=One(), β::Number=One()) = add!!(a + 0, b, α, β)


# Computes β * a + α * b. A zero result has the grades of b, plus the stored grades of a zero a unless β isa Zero.
# A nonzero β * a cancels against b, so its grades are among those of b.
function add!!(_a::NCAdd, b::Union{MulAdd,Number,UniformScaling}, α::Number=One(), β::Number=One())
    # scaling _a in place would also scale b if they alias, e.g. in scale!!(x, x, α)
    b = _a === b ? copy(b) : b
    a = scale!!(_a, β)
    agrades = β isa VectorInterface.Zero ? nothing : _zero_grades(a)
    return _with_zero_grades(_add_scaled!!(a, b, α), () -> _union_grades(agrades, _grades(b)))
end
# a + α * b without maintaining the grades of a zero result
function _add_scaled!!(a::NCAdd, b::NCMul, α)
    # a scalar product belongs in the coefficient, as in the NCAdd constructor
    isscalar(b) && return _add_scaled!!(a, prefactor(b), α)
    newdict = add_key!!(a.dict, term_key(b), α * prefactor(b))
    newcoeff = additive_coeff(a)
    if newdict === a.dict
        return set_coeff!!(a, newcoeff)
    end
    return NCAdd(newcoeff, newdict)
end
function _add_scaled!!(a::NCAdd, b::NCAdd, α)
    newdict = a.dict
    for (k, v) in b.dict
        newdict = add_key!!(newdict, k, α * v)
    end
    newcoeff = additive_coeff(a) + additive_coeff(b) * α
    if newdict === a.dict
        return set_coeff!!(a, newcoeff)
    end
    return NCAdd(newcoeff, newdict)
end
_add_scaled!!(a::NCAdd, b::Number, α) = set_coeff!!(a, additive_coeff(a) + α * b)
_add_scaled!!(a::NCAdd, b::UniformScaling, α) = set_coeff!!(a, additive_coeff(a) + α * b.λ)

# Sets `y` to `coeff` plus `terms`, key => value pairs that may read from `y.dict`. Everything is converted before `y`
# is mutated, so a failed conversion leaves `y` unchanged.
function _set_ncadd!(y::NCAdd{C,K}, coeff, terms) where {C,K}
    newcoeff = convert(C, coeff)
    newterms = [convert(K, k) => convert(C, v) for (k, v) in terms]
    y.coeff = newcoeff
    empty!(y.dict)
    for (k, v) in newterms
        iszero(v) || (y.dict[k] = v)
    end
    return y
end
function _set_ncadd!(y::NCAdd, x::NCAdd)
    _set_ncadd!(y, additive_coeff(x), x.dict)
    y.zero_grades = x.zero_grades
    return y
end
# The terms of `x` already have the types of `y`, so only the coefficient can fail to convert
function _set_ncadd!(y::NCAdd{C,K}, x::NCAdd{C,K}) where {C,K}
    newcoeff = convert(C, additive_coeff(x))
    if x.dict !== y.dict
        empty!(y.dict)
        for (k, v) in x.dict
            iszero(v) || (y.dict[k] = v)
        end
    end
    y.coeff = newcoeff
    y.zero_grades = x.zero_grades
    return y
end
_set_ncadd!(y::NCAdd, x::NCMul) = _set_ncadd!(y, NCAdd(x))

# If the result doesn't fit in the coefficient type, add! and scale! throw and leave `a` unchanged. An error from the
# arithmetic itself (e.g. an OverflowError with Rational{Int}) may leave `a` partly updated.
function add!(a::NCAdd{C}, b::Number, α::Number=One(), β::Number=One()) where {C}
    # β * a + α * b, with grades as in add!!, except that a scaled by zero contributes all its grades
    agrades = β isa VectorInterface.Zero ? nothing : iszero(β) ? _grades(a) : _zero_grades(a)
    newgrades = () -> _union_grades(agrades, _grades(b))
    newcoeff = additive_coeff(a) * β + α * b
    if !(promote_type(typeof(β), C) <: C)
        return _with_zero_grades(_set_ncadd!(a, newcoeff, (k => v * β for (k, v) in a.dict)), newgrades)
    end
    # the scaled terms fit in C, so only the coefficient can fail to convert
    set_coeff!(a, convert(C, newcoeff))
    if !(β isa One)
        map!(v -> v * β, values(a.dict))
        # a nonzero β can also give zero terms, through underflow
        filter_zeros!(a.dict)
    end
    return _with_zero_grades(a, newgrades)
end
scale!(x::NCAdd, α::Number) = add!(x, false, One(), α)

function scale!!(x::NCAdd{CS}, α::C) where {CS,C<:Number}
    # use scale! only where it can't fail
    if promote_type(C, CS) <: CS
        scale!(x, α)
    else
        scale(x, α)
    end
end
function scale!!(y::NCAdd, x::MulAdd, α::C) where C<:Number
    return add!!(y, x, α, VectorInterface.Zero())
end
scale(x::NCAdd, α::Number) = α * x

NCterms(a::NCAdd) = (NCMul(v, k.factors) for (k, v) in pairs(a.dict))

function Base.:*(x::Number, a::NCAdd{C,K}) where {C,K}
    dict = copy_dict(a, K, promote_type(typeof(x), C))
    map!(v -> x * v, values(dict))
    return _with_zero_grades(NCAdd(x * additive_coeff(a), dict), () -> _grades(a))
end
Base.:*(a::NCAdd, x::Number) = x * a

# a * b accumulated in a fresh zero sum, so a zero result has the grades of a and b only
Base.:*(a::NCAdd, b::MulAdd) = _shrink!(mul!!(zero(_sum_type(a, b)), a, b))
Base.:*(a::NCMul, b::NCAdd) = _shrink!(mul!!(zero(_sum_type(a, b)), a, b))
_shrink!(x::NCAdd) = (sizehint!(x.dict, length(x.dict)); x)
# Computes c + a * b. A zero result has the stored grades of c and the products of the grades of a and b.
function mul!!(c::NCAdd, a::MulAdd, b::MulAdd)
    cgrades = _zero_grades(c)
    return _with_zero_grades(_mul!!(c, a, b), () -> _union_grades(cgrades, _mul_grades(_grades(a), _grades(b))))
end
# c + a * b. With autosort, c is sorted and each term x of a is multiplied into it by mul_terms!!. Otherwise the
# catenated products are added unsorted. c is first converted to the promoted type of c, a and b, so the coefficient
# type of the result doesn't depend on which products vanish.
function _mul!!(c::NCAdd, a::MulAdd, b::MulAdd)
    autosort() && return _without_autosort(() -> _mul_sorted!!(c, a, b))
    c = convert(promote_type(typeof(c), _sum_type(a, b)), c)
    buffer = _factor_eltype(typeof(c))[]
    return fold_terms!!(c, a) do acc, x
        fold_terms!!((acc, y) -> add_term!!(acc, prefactor(x) * prefactor(y), _catenate!(buffer, x, y)), acc, b)
    end
end
function _mul_sorted!!(c::NCAdd, a::MulAdd, b::MulAdd)
    c = convert(promote_type(typeof(c), _sum_type(a, b)), _sorted(c))
    return fold_terms!!((acc, x) -> mul_terms!!(acc, x, b), c, a)
end
# `x` with its terms sorted, which is `x` itself if they already are
function _sorted(x::NCAdd)
    all(k -> _is_sorted(k.factors), keys(x.dict)) && return x
    return bubble_sort(x)
end
_is_sorted(factors) = all(i -> isnothing(_mul_effect(factors, factors[i], factors[i+1])), 1:length(factors)-1)

"""
    mul_terms!!(c, x::NCMul, b)

Add the sorted product of the term `x` and the product or sum `b` to the sum `c`, and return the result, which may be
`c` mutated or a new sum. `*` and `mul!!` call it for each term of the left factor when autosort is enabled.

Define methods for your own factor types to compute these products directly, e.g. to apply an operator string to every
basis state in a sum of kets with one precomputation. Use [`termpairs`](@ref) to read the terms of `b` and
[`add!!`](@ref) to add results to `c`. The result must equal `c + x * b` sorted, with a type that does not depend on
which products vanish.

Autosort is disabled while a `mul_terms!!` method runs, as it is for [`mul_effect`](@ref), so `*` inside the method
does not sort. A method must therefore produce terms that are already sorted and add them with `add!!`.
"""
function mul_terms!!(c, x::NCMul, b::MulAdd)
    buffer = _factor_eltype(typeof(c), typeof(x), typeof(b))[]
    return fold_terms!!((c, y) -> _mul_term_pair!!(c, x, y, buffer), c, b)
end
# c + x * y sorted, with the product formed in `buffer`
function _mul_term_pair!!(c, x, y, buffer)
    coeff = prefactor(x) * prefactor(y)
    iszero(coeff) && return c
    return sort_into!!(c, coeff, _catenate!(buffer, x, y))
end
_catenate!(buffer, x::NCMul, y::NCMul) = append!(empty!(buffer), x.factors, y.factors)

# The type of a sum of the terms of `xs`, with coefficients of the type of a sum of coefficients (Bool + Bool isa Int)
_sum_type(x) = _sum_type(typeof(x))
_sum_type(x, y) = _sum_type(promote_type(typeof(x), typeof(y)))
function _sum_type(::Type{T}) where {T<:MulAdd}
    C = _coeff_type(T)
    return NCAdd{Base.promote_op(+, C, C),to_add_dict_type(T)}
end
_coeff_type(::Type{<:NCMul{C}}) where {C} = C
_coeff_type(::Type{<:NCAdd{C}}) where {C} = C
# The element type of factors that can hold the factors of products and sums of the types `Ts`
_factor_eltype(Ts::Type...) = promote_type(map(_factor_eltype, Ts)...)
_factor_eltype(::Type{<:NCMul{<:Any,S}}) where {S} = S
_factor_eltype(::Type{<:NCAdd{<:Any,K}}) where {K} = _factor_eltype(K)
_factor_eltype(::Type) = Any
_factor_eltype(::Type{Union{}}) = Union{}

add!!(c::NCMul, term) = c + term
add!!(c::Number, term) = c + term

function Base.adjoint(x::NCAdd)
    adjoint_sum = autosort() ? _without_autosort(() -> _adjoint_sum(sort_into!!, x)) : _adjoint_sum(add_term!!, x)
    return _with_zero_grades(adjoint_sum, () -> _adjoint_grades(_grades(x)))
end
# The sum of the adjoints of the terms of x, each added to the sum with `add`
function _adjoint_sum(add::F, x::NCAdd) where {F}
    return fold_terms!!(zero(_sum_type(x)), x) do acc, t
        add(acc, adjoint(prefactor(t)), adjoint_factors(t.factors))
    end
end

Base.zero(::Type{NCAdd{C,K}}) where {C,K} = NCAdd{C,K}(zero(C), Dict{K,C}())
function Base.one(::Type{NCAdd{C,K}}) where {C,K}
    O = typeof(one(C))
    return NCAdd{O,K}(one(C), Dict{K,O}())
end
Base.oneunit(::T) where {T<:NCAdd} = oneunit(T)
Base.oneunit(::Type{NCAdd{C,K}}) where {C,K} = NCAdd{C,K}(oneunit(C), Dict{K,C}())
