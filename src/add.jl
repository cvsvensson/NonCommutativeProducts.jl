function filter_scalars!(d::AbstractDict{K,V}) where {K<:NCMul,V}
    coeff = zero(V)
    for (k, v) in d
        if isscalar(k)
            coeff += prefactor(k) * v
            delete!(d, k)
        end
    end
    return coeff
end
function filter_zeros!(d::AbstractDict{K,V}) where {K<:NCMul,V}
    for (k, v) in d
        if iszero(v)
            delete!(d, k)
        end
    end
    return d
end

mutable struct NCAdd{C,K,D}
    coeff::C
    dict::D
    function NCAdd(coeff::C, dict::D; keep_zeros=false) where {C,D<:AbstractDict}
        keep_zeros || filter_zeros!(dict)
        coeff += filter_scalars!(dict)
        new{promote_type(C, valtype(D)),keytype(D),D}(coeff, dict)
    end
end
NCAdd{C,K,D}(ncadd::NCAdd{C,K,D}) where {C,K,D} = ncadd

additive_coeff(x::NCAdd) = x.coeff
add_to_coeff!(a::NCAdd, x::Number) = a.coeff += x
function set_coeff!(a::NCAdd, x::Number)
    a.coeff = x
    return a
end
function set_coeff!!(a::NCAdd, x::Number)
    try
        set_coeff!(a, x)
    catch e
        NCAdd(x, a.dict)
    end
end
Base.convert(::Type{NCAdd{C,K,D}}, x::NCAdd) where {C,K,D} = NCAdd(convert(C, additive_coeff(x)), D(x.dict))
Base.convert(::Type{NCAdd{C,K,D}}, x::Number) where {C,K,D} = NCAdd(x, D())

const MulAdd = Union{NCMul,NCAdd}
function filter_scalars!(x::NCAdd)
    add!!(x, filter_scalars!(x.dict))
end
filter_zeros!(x::NCAdd) = (filter_zeros!(x.dict); return x)
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
Base.copy(x::NCAdd) = NCAdd(copy(additive_coeff(x)), copy(x.dict))

function print_coeff(io, coeff)
    if isreal(coeff)
        print(io, real(coeff), "I")
    else
        print(io, "(", coeff, ")", "I")
    end
end
function Base.show(io::IO, x::NCAdd; max_terms=3)
    compact = get(io, :compact, false)
    print_one = !iszero(additive_coeff(x)) || length(x.dict) == 0
    compact = length(x.dict) > max_terms
    print_sign(s) = compact ? print(io, s) : print(io, " ", s, " ")

    compact && println(io, "Sum with ", length(x.dict) + !iszero(additive_coeff(x)), " terms: ")
    N = min(max_terms, length(x.dict))

    print_one && print_coeff(io, additive_coeff(x))
    for (n, (k, v)) in enumerate(pairs(x.dict))
        n > max_terms && break
        should_print_sign = (n > 1 || print_one)
        if isreal(v)
            v = real(v)
            neg = v < 0
            if neg isa Bool
                if neg
                    print_sign("-")
                    print(io, -real(v) * k)
                else
                    should_print_sign && print_sign("+")
                    print(io, real(v) * k)
                end
            else
                should_print_sign && print_sign("+")
                print(io, "(", v, ")*", k)
            end
        else
            should_print_sign && print_sign("+")
            print(io, "(", v, ")*", k)
        end
    end
    if N < length(x.dict)
        print(io, " + ...")
    end
    return nothing
end

Base.:+(a::Number, b::NCAdd) = iszero(a) ? b : NCAdd(a + additive_coeff(b), b.dict)
Base.:+(a::UniformScaling, b::NCAdd) = iszero(a) ? b : NCAdd(a.λ + additive_coeff(b), b.dict)
Base.:+(a::NCAdd, b::B) where B<:Union{Number,UniformScaling} = b + a
Base.:+(a::NCAdd, b::B) where B<:NCMul = b + a
Base.:/(a::MulAdd, b::Number) = inv(b) * a
Base.:-(a::Union{Number,UniformScaling}, b::MulAdd) = a + (-b)
Base.:-(a::MulAdd, b::Union{Number,MulAdd,UniformScaling}) = a + (-b)
Base.:-(a::NCAdd) = (-1) * a
function Base.:+(a::NCAdd, b::NCAdd)
    coeff = additive_coeff(a) + additive_coeff(b)
    dict = mergewith(+, a.dict, b.dict)
    NCAdd(coeff, dict)
end
add!!(a::NCMul, b::MulAdd, α::Number=One(), β::Number=One()) = add!!(a + 0, b, α, β)

function add!!(a::NCAdd, b::NCMul, α::Number=One(), β::Number=One())
    # compute β * a + α * b
    key = NCMul(1, b.factors)
    coeff = α * prefactor(b)
    newdict, ret = modify!!(a.dict, key) do val
        isnothing(val) && return coeff
        return something(val, 0) * β + coeff
    end
    newcoeff = additive_coeff(a) * β
    if newdict === a.dict
        return set_coeff!!(a, newcoeff)
    end
    return NCAdd(newcoeff, newdict)
end
function add!!(a::NCAdd, b::NCAdd, α::Number=One(), β::Number=One())
    newdict = a.dict
    for (k, v) in b.dict
        coeff = α * v
        newdict, ret = modify!!(a.dict, k) do val
            isnothing(val) && return coeff
            return something(val, 0) * β + coeff
        end
    end
    newcoeff = additive_coeff(a) * β + additive_coeff(b) * α
    if newdict === a.dict
        return set_coeff!!(a, newcoeff)
    end
    return NCAdd(newcoeff, newdict)
end
function add!!(_a::NCAdd, b::Number, α::Number=One(), β::Number=One())
    a = scale!!(_a, β)
    set_coeff!!(a, additive_coeff(a) + α * b)
end
function add!!(_a::NCAdd, b::UniformScaling, α::Number=One(), β::Number=One())
    a = scale!!(_a, β)
    set_coeff!!(a, additive_coeff(a) + α * b.λ)
end

function scale!(x::NCAdd, α::Number)
    x.coeff *= α
    for (k, v) in x.dict
        x.dict[k] = v * α
    end
    return x
end
function scale!!(x::NCAdd, α::Number)
    x = try
        scale!(x, α)
    catch e
        scale(x, α)
    end
end
scale(x::NCAdd, α::Number) = α * x

function NCterms(a::NCAdd)
    (v * k for (k, v) in pairs(a.dict))
end

function Base.:*(x::Number, a::NCAdd)
    acoeff = additive_coeff(a)
    return NCAdd(x * acoeff, Dict(k => v * x for (k, v) in a.dict))
end
Base.:*(a::NCAdd, x::Number) = x * a

function Base.:*(a::NCAdd, b::NCMul)
    c = zero(a)
    filter_scalars!(filter_zeros!(mul!!(c, a, b)))
end
function Base.:*(a::NCMul, b::NCAdd)
    c = zero(b)
    filter_scalars!(filter_zeros!(mul!!(c, a, b)))
end
function Base.:*(a::NCAdd, b::NCAdd)
    c = zero(a)
    filter_scalars!(filter_zeros!(mul!!(c, a, b)))
end
function mul!!(c::NCAdd, a::MulAdd, b::MulAdd)
    acoeff = additive_coeff(a)
    bcoeff = additive_coeff(b)
    if !iszero(acoeff)
        c = add!!(c, b, acoeff, One())
    end
    if !iszero(bcoeff)
        c = add!!(c, a, bcoeff, One())
        c = add!!(c, -acoeff * bcoeff) # We've double counted this term so subtract it
    end
    for bterm in NCterms(b)
        for aterm in NCterms(a)
            newterm = catenate(aterm, bterm)
            c = add!!(c, newterm)
        end
    end
    if autosort()
        return sort!(c)
    end
    return c
end
add!!(c::NCMul, term) = c + term
add!!(c::Number, term) = c + term

function Base.adjoint(x::NCAdd)
    newx = zero(x)
    set_coeff!(newx, adjoint(additive_coeff(x)))
    for (f, v) in x.dict
        newx = add!!(newx, v' * f')
    end
    newx
end

Base.zero(::Type{NCAdd{C,K,D}}) where {C,K,D} = NCAdd(zero(C), D())
Base.one(::Type{NCAdd{C,K,D}}) where {C,K,D} = NCAdd(one(C), D())

@testitem "Consistency between + and add!!" setup = [Fermions] begin
    import NonCommutativeProducts: add!!
    NonCommutativeProducts.disable_autosort!()
    f = Fermion.(1:2)
    a = 1.0 * f[2] * f[1] + 1 + f[1]
    for b in [1.0, 1, f[1], 1.0 * f[1], f[2] * f[1], a]
        a2 = copy(a)
        a3 = add!!(a2, b) # Should mutate
        @test a + b == a3
        @test a2 == a3
        anew = add!!(a, 1im * b) #Should not mutate
        @test a2 !== anew
    end
    @test a == 1.0 * f[2] * f[1] + 1 + f[1]

    NonCommutativeProducts.disable_autosort!()
    @test add!!(f[1] + 1, 1) == f[1] + 2
    @test add!!(f[1] + 1, f[1]) == 2f[1] + 1
    @test add!!(f[1] + 1, 5 * f[1]) == 6f[1] + 1
    @test add!!(f[1] + 1, 1, 2, 5) == 5f[1] + 7
    @test add!!(f[1] + 1, f[1], 2, 5) == 7f[1] + 5
    @test add!!(f[1] + 1, 5 * f[1], 2, 5) == 15f[1] + 5

end
