
# The factors of an NCMul may be shared with other products and with the keys of sums, e.g. by `2 * x` or `x + y`.
# They must therefore never be mutated, except by bubble_sort! on a product whose factors were just created.
"""
    NCMul

A product: a coefficient times an ordered list of factors. Create products by multiplying atoms registered with
`@nc`, e.g. `2 * a * b`. Use `coefficient(x)` and `factors(x)` to read the parts of a product.
"""
struct NCMul{C,S,F}
    coeff::C
    factors::F
    function NCMul{C,S,F}(coeff, factors) where {C,S,F}
        S === eltype(F) || throw(ArgumentError("the factor type S = $S of NCMul must be eltype(F) = $(eltype(F))"))
        new{C,S,F}(coeff, factors)
    end
end
NCMul{C,S,F}(ncmul::NCMul{C,S,F}) where {C,S,F} = ncmul
NCMul{C}(coeff, factors) where C = NCMul{C,eltype(factors),typeof(factors)}(coeff, factors)
NCMul(coeff, factors) = NCMul{typeof(coeff)}(coeff, factors)

# Tuple factors are for small products that are cheap to construct. Operations that multiply, sort or mutate collect
# them into a Vector, so their results, as well as zero, one and oneunit, are Vector-backed.
factors_vector(factors::AbstractVector) = factors
factors_vector(factors::Tuple) = collect(eltype(factors), factors)
copy_factors(factors) = copy(factors)
copy_factors(factors::Tuple) = factors_vector(factors)

Base.zero(::Type{NCMul{C,S,F}}) where {C,S,F} = NCMul(zero(C), S[])
Base.one(::Type{NCMul{C,S,F}}) where {C,S,F} = NCMul(one(C), S[])
Base.oneunit(::T) where {T<:NCMul} = oneunit(T)
Base.oneunit(::Type{NCMul{C,S,F}}) where {C,S,F} = NCMul(oneunit(C), S[])
Base.copy(x::NCMul) = NCMul(copy(prefactor(x)), copy_factors(x.factors))
# By value: a zero product is the scalar 0, whatever factors it still carries
isscalar(x::NCMul) = isempty(x.factors) || iszero(x)
"""
    scalar(x)

The number that `x` equals if `x` is a scalar, i.e. a number or a product or sum without factors. Throws an
`ArgumentError` for other products and sums.
"""
scalar(x::NCMul) = isscalar(x) ? prefactor(x) : throw(ArgumentError("NCMul is not a scalar"))
scalar(x::Number) = x
additive_coeff(::NCMul{C}) where C = zero(C)
additive_coeff(::NCMul{Any}) = 0
prefactor(x::NCMul) = x.coeff

function Base.show(io::IO, x::NCMul)
    print_coeff = !isone(prefactor(x)) || isempty(x.factors)
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

# Equality is by value, so all zero products are equal, to each other and to 0, whatever their factors
Base.:(==)(a::NCMul, b::Number) = isscalar(a) && prefactor(a) == b
Base.:(==)(a::Number, b::NCMul) = b == a
Base.:(==)(a::NCMul, b::NCMul) = iszero(a) ? iszero(b) : prefactor(a) == prefactor(b) && factors_equal(a.factors, b.factors)
# Equality and hashing ignore the factors container, so that Tuple-backed and Vector-backed products
# can be used interchangeably as keys of the same Dict
factors_equal(a::F, b::F) where {F} = a == b
factors_equal(a, b) = length(a) == length(b) && all(splat(==), zip(a, b))
hash_factors(factors, h::UInt) = foldl((h, f) -> hash(f, h), factors; init=hash(length(factors), h))
function Base.hash(a::NCMul, h::UInt)
    # by value, consistent with ==
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
        return _with_zero_grades(bubble_sort!(ncmul), () -> _mul_grades(_grades(a), _grades(b)))
    end
    return ncmul
end
catenate(a::NCMul, b::NCMul) = NCMul(prefactor(a) * prefactor(b), vcat(factors_vector(a.factors), factors_vector(b.factors)))

function Base.adjoint(x::NCMul)
    length(x.factors) == 0 && return NCMul(adjoint(prefactor(x)), x.factors)
    ncmul = NCMul(adjoint(prefactor(x)), adjoint_factors(x.factors))
    if autosort()
        return _with_zero_grades(bubble_sort!(ncmul), () -> _adjoint_grades(_grades(x)))
    end
    return ncmul
end
adjoint_factors(factors) = collect(Iterators.reverse(Iterators.map(adjoint, factors)))
# keep mixed-type products on Vector{Any} instead of narrowing to the eltype of this particular product
adjoint_factors(factors::Vector{Any}) = Any[adjoint(f) for f in Iterators.reverse(factors)]
