abstract type AbstractOrdering end
struct NaiveOrdering <: AbstractOrdering end

struct NCMul{C,S}
    coeff::C
    factors::Vector{S}
    function NCMul(coeff::C, factors) where {C}
        new{C,eltype(factors)}(coeff, factors)
    end
end
ordered_product(a, b, ::NaiveOrdering) = NCMul(1, [a, b])
Base.convert(::Type{NCMul{C,S}}, x::NCMul{<:Any,S}) where {C,S} = NCMul(convert(C, x.coeff), x.factors)
function canonicalize!(a::NCMul)
    if iszero(a.coeff)
        return a.coeff
    elseif length(a.factors) == 1 && isone(a.coeff)
        return only(a.factors)
    end
    return a
end

function Base.show(io::IO, x::NCMul)
    isscalar(x) && print(io, x.coeff)
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


ordered_product(x::Number, a::NCMul, ordering) = NCMul(x * a.coeff, a.factors)
ordered_product(a::NCMul, b::NCMul, ::NaiveOrdering) = NCMul(a.coeff * b.coeff, vcat(a.factors, b.factors))

ordered_product(bs::NCMul, ordering) = ordered_product(NCMul(a), bs, ordering)
ordered_product(as::NCMul, b, ordering) = ordered_product(as, NCMul(b), ordering)

Base.adjoint(x::NCMul) = length(x.factors) == 0 ? NCMul(adjoint(x.coeff), x.factors) : adjoint(x.coeff) * foldr(*, Iterators.reverse(Iterators.map(adjoint, x.factors)))

# isscalar(x::NCMul) = iszero(x.coeff) || (length(x.factors) == 0)
