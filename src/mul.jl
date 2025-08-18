abstract type AbstractOrdering end
struct NaiveOrdering <: AbstractOrdering end

mutable struct NCMul{C,S,F}
    coeff::C
    factors::F
    function NCMul(coeff::C, factors::F) where {C,F}
        new{C,eltype(factors),F}(coeff, factors)
    end
end
Base.convert(::Type{NCMul{C,S,F}}, x::NCMul{<:Any,S,F}) where {C,S,F} = NCMul(convert(C, x.coeff), x.factors)
Base.copy(x::NCMul) = NCMul(copy(x.coeff), copy(x.factors))
isscalar(x::NCMul) = length(x.factors) == 0


function Base.show(io::IO, x::NCMul)
    #isscalar(x) && print(io, x.coeff)
    print_coeff = !isone(x.coeff)
    if print_coeff
        v = x.coeff
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
Base.iszero(x::NCMul) = iszero(x.coeff)

Base.:(==)(a::NCMul, b::Number) = isscalar(a) && a.coeff == b
Base.:(==)(a::NCMul, b::NCMul) = a.coeff == b.coeff && a.factors == b.factors
Base.hash(a::NCMul, h::UInt) = hash(a.coeff, hash(a.factors, h))
NCMul(f::NCMul) = f
NCMul(f) = NCMul(1, [f])

NCterms(a::NCMul) = (a,)
Base.:-(a::NCMul) = NCMul(-a.coeff, a.factors)

Base.:*(x::Number, a::NCMul) = NCMul(x * a.coeff, a.factors)
Base.:*(m::NCMul, x::Number) = x * m

ordered_product(x::Number, a::NCMul, ordering) = NCMul(x * a.coeff, a.factors)
ordered_product(a::NCMul, b::NCMul, ::NaiveOrdering) = NCMul(a.coeff * b.coeff, vcat(a.factors, b.factors))

ordered_product(a::NCMul, bs::NCMul, ordering) = ordered_product(a, bs, ordering)

Base.adjoint(x::NCMul) = length(x.factors) == 0 ? NCMul(adjoint(x.coeff), x.factors) : adjoint(x.coeff) * foldr(*, Iterators.reverse(Iterators.map(adjoint, x.factors)))

