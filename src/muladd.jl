NCMul(f::NCAdd) = (length(f.dict) == 1 && iszero(f.coeff) && return prod(only(f.dict))) || throw(ArgumentError("Cannot convert NCAdd to NCMul: $f"))


function Base.:(==)(a::NCAdd, b::NCMul)
    iszero(a.coeff) || return false
    length(a.dict) == 1 || return false
    ncmul, coeff = only(a.dict)
    NCMul(coeff, ncmul.factors) == b
end
Base.:(==)(a::NCMul, b::NCAdd) = b == a

Base.:+(a::NCMul, b) = a + NCMul(b)
Base.:+(a, b::NCMul) = NCMul(a) + b
function Base.:+(a::NCMul{CA,KA}, b::NCMul{CB,KB}) where {CA,CB,KA,KB}
    C = promote_type(CA, CB)
    K = Union{NCMul{C,KA},NCMul{C,KB}}
    D = Dict{K,C}
    if a.factors == b.factors
        coeff = a.coeff + b.coeff
        return (NCAdd(0, D(NCMul(1, a.factors) => coeff)))
    end
    at, bt = to_add_tuple(a), to_add_tuple(b)
    return (NCAdd(0, D(at..., bt...)))
end

to_add(a::NCMul, coeff=1) = Dict(NCMul(1, a.factors) => a.coeff * coeff)
to_add(a, coeff=1) = Dict(NCMul(a) => coeff)
to_add_tuple(a::NCMul, coeff=1) = (NCMul(1, a.factors) => a.coeff * coeff,)
to_add_tuple(a, coeff=1) = (NCMul(a) => coeff,)

Base.:^(a::Union{NCMul,NCAdd}, b) = Base.power_by_squaring(a, b)

macro nc(T)
    quote
        Base.:+(x::$(esc(T)), y::$(esc(T))) = NCMul(1, [x]) + NCMul(1, [y])
        Base.:+(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = NCMul(1, [x]) + y
        Base.:+(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = x + NCMul(1, [y])

        Base.:-(x::$(esc(T)), y::$(esc(T))) = NCMul(1, [x]) - NCMul(1, [y])
        Base.:-(x::$(esc(T))) = NCMul(-1, [x])
        Base.:-(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = x - NCMul(1, [y])
        Base.:-(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = NCMul(1, [x]) - y


        Base.:*(x::$(esc(T)), y::$(esc(T))) = NCMul(1, [x, y])
        Base.:*(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = x * NCMul(1, [y])
        Base.:*(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = NCMul(1, [x]) * y
        Base.:*(x::$(esc(T)), y::NCMul) = NCMul(y.coeff, pushfirst!!(copy(y.factors), x))
        Base.:*(x::NCMul, y::$(esc(T))) = NCMul(x.coeff, push!!(copy(x.factors), y))
        Base.:^(a::$(esc(T)), b) = Base.power_by_squaring(a, b)
        Base.convert(::Type{NCMul{C,S}}, x::$(esc(T))) where {C,S<:$(esc(T))} = NCMul(one(C), [x])
    end
end
# Base.convert(::Type{NCMul{C,F}}, x::Fermion) where {C,F<:Fermion} = NCMul(one(C), [x])


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

