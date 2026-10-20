# The rules rely on two invariants:
#   NCMul{C,S,F}: S == eltype(F)
#   NCAdd{C,K}: the keys K are NCMul{Int,S,F}, and the terms are a Dict{K,C}
# Coefficient and factor types are promoted independently. Each rule is defined for one argument order only, since
# promote_type tries both.

# promote_type(Vector{A}, Vector{B}) is the abstract `Vector` when A and B promote to Any, so build the type explicitly.
# Other mixes, such as a Tuple with a Vector, also give Vector{S}.
promote_factors_type(::Type{S}, ::Type{<:Vector}, ::Type{<:Vector}) where {S} = Vector{S}
function promote_factors_type(::Type{S}, ::Type{F1}, ::Type{F2}) where {S,F1,F2}
    F = promote_type(F1, F2)
    return isconcretetype(F) && eltype(F) == S ? F : Vector{S}
end

# convert has no method between Tuples and Vectors
convert_factors(::Type{F}, factors) where {F} = convert(F, factors)
convert_factors(::Type{F}, factors::Tuple) where {F<:AbstractVector} = convert(F, collect(factors))
convert_factors(::Type{F}, factors::AbstractVector) where {F<:Tuple} = convert(F, Tuple(factors))

function Base.promote_rule(::Type{NCMul{C1,S1,F1}}, ::Type{NCMul{C2,S2,F2}}) where {C1,S1,F1,C2,S2,F2}
    S = promote_type(S1, S2)
    return NCMul{promote_type(C1, C2),S,promote_factors_type(S, F1, F2)}
end
function Base.promote_rule(::Type{NCAdd{C1,K1}}, ::Type{NCAdd{C2,K2}}) where {C1,K1,C2,K2}
    return NCAdd{promote_type(C1, C2),promote_type(K1, K2)}
end
function Base.promote_rule(::Type{NCMul{C1,S1,F1}}, ::Type{NCAdd{C2,K2}}) where {C1,S1,F1,C2,K2}
    return NCAdd{promote_type(C1, C2),promote_type(to_add_dict_type(NCMul{C1,S1,F1}), K2)}
end

Base.convert(::Type{NCMul{C,S,F}}, x::NCMul{C,S,F}) where {C,S,F} = x
Base.convert(::Type{NCMul{C,S,F}}, x::NCMul) where {C,S,F} = NCMul{C,S,F}(convert(C, prefactor(x)), convert_factors(F, x.factors))

Base.convert(::Type{NCAdd{C,K}}, x::NCAdd{C,K}) where {C,K} = x
function Base.convert(::Type{NCAdd{C,K}}, x::NCAdd) where {C,K}
    dict = Dict{K,C}(convert(K, k) => convert(C, v) for (k, v) in pairs(x.dict))
    NCAdd{C,K}(additive_coeff(x), dict; zero_grades=x.zero_grades)
end
function Base.convert(::Type{NCAdd{C,K}}, x::NCMul) where {C,K}
    ncadd = NCAdd{C,K}(zero(C), Dict{K,C}(term_key(K, x) => convert(C, prefactor(x))))
    return _with_zero_grades(ncadd, () -> _grades(x))
end
Base.convert(::Type{NCAdd{C,K}}, x::Number) where {C,K} = NCAdd{C,K}(x, Dict{K,C}())

# A copy of the terms of `x` that can hold keys of type K and coefficients of type C
function copy_dict(x::NCAdd, ::Type{K}, ::Type{C}) where {K,C}
    keytype(x.dict) === K && valtype(x.dict) === C && return copy(x.dict)
    return Dict{K,C}(x.dict)
end
