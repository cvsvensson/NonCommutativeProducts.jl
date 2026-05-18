
struct NCMul{C,S,F}
    coeff::C
    factors::F
end
function NCMul(coeff::C, factors::F) where {C,F}
    NCMul{C,eltype(factors),F}(coeff, factors)
end
Base.convert(::Type{NCMul{C,S,F}}, x::NCMul{<:Any,S,F}) where {C,S,F} = NCMul{C,S,F}(convert(C, prefactor(x)), x.factors)
Base.convert(::Type{NCMul{C,S,F}}, x::NCMul) where {C,S,F} = NCMul{C,S,F}(convert(C, prefactor(x)), F(x.factors))
Base.promote_rule(::Type{NCMul{C,S,F}}, x::Type{NCMul{C2,S2,F2}}) where {C,S,F,C2,S2,F2} = NCMul{promote_type(C, C2),promote_type(S, S2),promote_type(F, F2)}

Base.zero(::Type{NCMul{C,S,F}}) where {C,S,F} = NCAdd(zero(C), Dict{NCMul{Int,S,F},C}())
Base.copy(x::NCMul) = NCMul(copy(prefactor(x)), copy(x.factors))
isscalar(x::NCMul) = length(x.factors) == 0
additive_coeff(::NCMul{C}) where C = zero(C)
additive_coeff(::NCMul{Any}) = 0
prefactor(x::NCMul) = getfield(x, :coeff)

function Base.show(io::IO, x::NCMul)
    print_coeff = !isone(prefactor(x)) || isscalar(x)
    if print_coeff
        v = prefactor(x)
        if isreal(v)
            neg = real(v) < 0
            if neg isa Bool
                print(io, real(v))
            else
                print(io, "(", v, ")")
            end
        else
            print(io, "(", v, ")")
        end
    end
    for (n, x) in enumerate(x.factors)
        if print_coeff || n > 1
            print(io, "*")
        end
        print(io, x)
    end
end
Base.iszero(x::NCMul) = iszero(prefactor(x))

Base.:(==)(a::NCMul, b::Number) = (isscalar(a) && prefactor(a) == b) || iszero(a) && iszero(b)
Base.:(==)(a::NCMul, b::NCMul) = prefactor(a) == prefactor(b) && a.factors == b.factors
Base.hash(a::NCMul, h::UInt) = isone(prefactor(a)) && length(a.factors) == 1 ? hash(only(a.factors), h) : hash(prefactor(a), hash(a.factors, h))
NCMul(f::NCMul) = f

NCterms(a::NCMul) = (a,)
Base.:-(a::NCMul) = NCMul(-prefactor(a), a.factors)

Base.:*(x::Number, a::NCMul) = NCMul(x * prefactor(a), a.factors)
Base.:*(m::NCMul, x::Number) = x * m
function Base.:*(a::NCMul, b::NCMul)
    ncmul = catenate(a, b)
    if autosort()
        return bubble_sort!(ncmul)
    end
    return ncmul
end
catenate(x::NCMul, others...) = NCMul(prefactor(x) * prod(y -> prefactor(y), others), vcat(x.factors, map(y -> y.factors, others)...))

function Base.adjoint(x::NCMul)
    length(x.factors) == 0 && return NCMul(adjoint(prefactor(x)), x.factors)
    ncmul = NCMul(adjoint(prefactor(x)), collect(Iterators.reverse(Iterators.map(adjoint, x.factors))))
    if autosort()
        return bubble_sort!(ncmul)
    end
    return ncmul
end