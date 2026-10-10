
# The factors of an NCMul may be shared with other products and with the keys of sums, e.g. by `2 * x` or `x + y`.
# They must therefore never be mutated, except by bubble_sort! on a product whose factors were just created.
struct NCMul{C,S,F}
    coeff::C
    factors::F
end
NCMul{C,S,F}(ncmul::NCMul{C,S,F}) where {C,S,F} = ncmul

function NCMul(coeff::C, factors::F) where {C,F}
    NCMul{C,eltype(factors),F}(coeff, factors)
end

# Tuple factors are meant for small products that are cheap to construct. Only cheap operations such as
# equality, hashing and adjoint accept them as they are. Operations that multiply, sort or mutate collect
# the factors into a Vector first, so their results, as well as zero and one, are Vector-backed.
factors_vector(factors::AbstractVector) = factors
factors_vector(factors::Tuple) = collect(eltype(factors), factors)
copy_factors(factors) = copy(factors)
copy_factors(factors::Tuple) = factors_vector(factors)

Base.zero(::Type{NCMul{C,S,F}}) where {C,S,F} = NCMul(zero(C), S[])
Base.one(::Type{NCMul{C,S,F}}) where {C,S,F} = NCMul(one(C), S[])
Base.oneunit(::T) where {T<:NCMul} = oneunit(T)
Base.oneunit(::Type{NCMul{C,S,F}}) where {C,S,F} = NCMul(one(C), S[])
Base.copy(x::NCMul) = NCMul(copy(prefactor(x)), copy_factors(x.factors))
isscalar(x::NCMul) = length(x.factors) == 0
scalar(x::NCMul) = isscalar(x) ? prefactor(x) : throw(ArgumentError("NCMul is not a scalar"))
scalar(x::Number) = x
additive_coeff(::NCMul{C}) where C = zero(C)
additive_coeff(::NCMul{Any}) = 0
prefactor(x::NCMul) = x.coeff

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
Base.:(==)(a::Number, b::NCMul) = b == a
Base.:(==)(a::NCMul, b::NCMul) = prefactor(a) == prefactor(b) && factors_equal(a.factors, b.factors)
# Equality and hashing ignore the factors container, so that Tuple-backed and Vector-backed products
# can be used interchangeably as keys of the same Dict
factors_equal(a::F, b::F) where {F} = a == b
factors_equal(a, b) = length(a) == length(b) && all(splat(==), zip(a, b))
hash_factors(factors, h::UInt) = foldl((h, f) -> hash(f, h), factors; init=hash(length(factors), h))
function Base.hash(a::NCMul, h::UInt)
    # scalar and zero products are == to numbers, so they must hash like them
    iszero(a) && return hash(zero(prefactor(a)), h)::UInt
    isscalar(a) && return hash(prefactor(a), h)::UInt
    single_term = isone(prefactor(a)) && length(a.factors) == 1
    if single_term
        hash(only(a.factors), h)::UInt
    else
        hash(prefactor(a), hash_factors(a.factors, h))::UInt
    end
end
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
catenate(x::NCMul, others...) = NCMul(prefactor(x) * prod(y -> prefactor(y), others), vcat(factors_vector(x.factors), map(y -> factors_vector(y.factors), others)...))

function Base.adjoint(x::NCMul)
    length(x.factors) == 0 && return NCMul(adjoint(prefactor(x)), x.factors)
    ncmul = NCMul(adjoint(prefactor(x)), adjoint_factors(x.factors))
    if autosort()
        return bubble_sort!(ncmul)
    end
    return ncmul
end
adjoint_factors(factors) = collect(Iterators.reverse(Iterators.map(adjoint, factors)))
# keep mixed-type products on Vector{Any} instead of narrowing to the eltype of this particular product
adjoint_factors(factors::Vector{Any}) = Any[adjoint(f) for f in Iterators.reverse(factors)]

isfilterable(x::NCMul) = all(isfilterable, x.factors)