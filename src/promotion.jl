# Promotion and conversion between NCMul, NCAdd, and the atoms registered with @nc.
#
# The rules rely on two invariants:
#   NCMul{C,S,F}: S == eltype(F)
#   NCAdd{C,K,D}: the keys K are NCMul{Int,S,F}, and valtype(D) == C (enforced by the NCAdd constructor)
#
# Coefficient types and factor types are promoted independently. A product promotes with a sum as the sum
# with the product's factors as a key, and an atom of type W promotes as the product NCMul{Int,W,Vector{W}}
# (see @nc_common). Each rule is defined for one argument order only, since promote_type tries both.

# promote_type(Vector{A}, Vector{B}) typejoins to the abstract `Vector` when A and B promote to Any, so build the container type explicitly.
promote_factors_type(::Type{S}, ::Type{<:Vector}, ::Type{<:Vector}) where {S} = Vector{S}
promote_factors_type(::Type{S}, ::Type{F1}, ::Type{F2}) where {S,F1,F2} = promote_type(F1, F2)

ncadd_type(::Type{C}, ::Type{K}) where {C,K} = NCAdd{C,K,Dict{K,C}}

function Base.promote_rule(::Type{NCMul{C1,S1,F1}}, ::Type{NCMul{C2,S2,F2}}) where {C1,S1,F1,C2,S2,F2}
    S = promote_type(S1, S2)
    return NCMul{promote_type(C1, C2),S,promote_factors_type(S, F1, F2)}
end
function Base.promote_rule(::Type{NCAdd{C1,K1,D1}}, ::Type{NCAdd{C2,K2,D2}}) where {C1,K1,D1,C2,K2,D2}
    return ncadd_type(promote_type(C1, C2), promote_type(K1, K2))
end
function Base.promote_rule(::Type{NCMul{C1,S1,F1}}, ::Type{NCAdd{C2,K2,D2}}) where {C1,S1,F1,C2,K2,D2}
    return ncadd_type(promote_type(C1, C2), promote_type(NCMul{Int,S1,F1}, K2))
end

Base.convert(::Type{NCMul{C,S,F}}, x::NCMul) where {C,S,F} = NCMul{C,S,F}(convert(C, prefactor(x)), convert(F, x.factors))

Base.convert(::Type{NCAdd{C,K,D}}, x::NCAdd{C,K,D}) where {C,K,D} = x
function Base.convert(::Type{NCAdd{C,K,D}}, x::NCAdd) where {C,K,D}
    dict = D(convert(K, k) => convert(valtype(D), v) for (k, v) in pairs(x.dict))
    NCAdd(convert(C, additive_coeff(x)), dict)
end
function Base.convert(::Type{NCAdd{C,K,D}}, x::NCMul) where {C,K,D}
    NCAdd(zero(C), D(convert(K, NCMul(1, x.factors)) => convert(valtype(D), prefactor(x))))
end
Base.convert(::Type{NCAdd{C,K,D}}, x::Number) where {C,K,D} = NCAdd(convert(C, x), D())

# A copy of the terms of `x` that can hold keys of type K and coefficients of type C
function copy_dict(x::NCAdd, ::Type{K}, ::Type{C}) where {K,C}
    keytype(x.dict) === K && valtype(x.dict) === C && return copy(x.dict)
    return Dict{K,C}(x.dict)
end
