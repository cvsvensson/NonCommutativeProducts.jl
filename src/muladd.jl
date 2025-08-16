NCMul(f::NCAdd) = (length(f.dict) == 0 && iszero(f.coeff) && return prod(only(f.dict))) || throw(ArgumentError("Cannot convert NCAdd to NCMul: $f"))




Base.:+(a::NCMul, b) = a + NCMul(b)
Base.:+(a, b::NCMul) = NCMul(a) + b
function Base.:+(a::NCMul{CA,KA}, b::NCMul{CB,KB}) where {CA,CB,KA,KB}
    C = promote_type(CA, CB)
    K = Union{NCMul{C,KA},NCMul{C,KB}}
    D = Dict{K,C}
    if a.factors == b.factors
        coeff = a.coeff + b.coeff
        return canonicalize!(NCAdd(0, D(NCMul(1, a.factors) => coeff)))
    end
    at, bt = to_add_tuple(a), to_add_tuple(b)
    return canonicalize!(NCAdd(0, D(at..., bt...)))
end

to_add(a::NCMul, coeff=1) = Dict(NCMul(1, a.factors) => a.coeff * coeff)
to_add(a, coeff=1) = Dict(NCMul(a) => coeff)
to_add_tuple(a::NCMul, coeff=1) = (NCMul(1, a.factors) => a.coeff * coeff,)
to_add_tuple(a, coeff=1) = (NCMul(a) => coeff,)

Base.:^(a::Union{NCMul,NCAdd}, b) = Base.power_by_squaring(a, b)




##
TermInterface.head(a::Union{NCMul,NCAdd}) = operation(a)
TermInterface.iscall(::Union{NCMul,NCAdd}) = true
TermInterface.isexpr(::Union{NCMul,NCAdd}) = true

TermInterface.operation(::NCMul) = (*)
TermInterface.operation(::NCAdd) = (+)
TermInterface.arguments(a::NCMul) = [a.coeff, a.factors...]
TermInterface.arguments(a::NCAdd) = iszero(a.coeff) ? NCterms(a) : allterms(a)
TermInterface.sorted_arguments(a::NCAdd) = iszero(a.coeff) ? sort(NCterms(a), by=x -> x.factors) : [a.coeff, sort(NCterms(a); by=x -> x.factors)...]
TermInterface.children(a::Union{NCMul,NCAdd}) = arguments(a)
TermInterface.sorted_children(a::Union{NCMul,NCAdd}) = sorted_arguments(a)

TermInterface.maketerm(::Type{<:NCMul}, ::typeof(*), args, metadata) = *(args...)
TermInterface.maketerm(::Type{<:NCAdd}, ::typeof(+), args, metadata) = +(args...)

TermInterface.maketerm(::Type{Q}, head::Type{T}, args, metadata) where {Q<:Union{<:NCMul,<:NCAdd},T} = T(args...)

