function filter_ncadd_dict!(d::AbstractDict{K,V}; filter_zeros=true, filter_scalars=true) where {K<:NCMul,V}
    coeff = zero(V)
    !filter_zeros && !filter_scalars && return d, coeff
    for (k, v) in d
        if filter_zeros && dropzero(k, v)
            delete!(d, k)
            continue
        end
        if filter_scalars && isscalar(k)
            delete!(d, k)
            coeff += prefactor(k) * v
        end
    end
    return d, coeff
end

# Zero terms are dropped unless their key is not `isfilterable`, the same rule `filter_ncadd_dict!` uses.
# Every path that mutates the terms of an NCAdd in place goes through this, so that in-place and out-of-place
# results compare and hash equal.
dropzero(k, v) = iszero(v) && isfilterable(k)
filter_zeros!(d::AbstractDict) = filter!(kv -> !dropzero(first(kv), last(kv)), d)

# D is always Dict{K,C}: the terms share the coefficient type of the sum, and Dict is the only container, since
# equality and key lookup assume keys are compared by value. D is kept as a parameter for backwards compatibility.
mutable struct NCAdd{C,K,D<:Dict{K,C}}
    coeff::C
    dict::D
    function NCAdd(coeff::C, dict::D; kwargs...) where {C,D<:AbstractDict}
        _, addcoeff = filter_ncadd_dict!(dict; kwargs...)
        newcoeff = coeff + addcoeff
        T = promote_type(typeof(newcoeff), valtype(D))
        newdict = dict isa Dict{keytype(D),T} ? dict : Dict{keytype(D),T}(dict)
        new{T,keytype(D),typeof(newdict)}(newcoeff, newdict)
    end
end
NCAdd{C,K,D}(ncadd::NCAdd{C,K,D}) where {C,K,D} = ncadd
NCAdd(ncmul::NCMul{C}) where C = NCAdd(zero(C), to_add_dict(ncmul))

additive_coeff(x::NCAdd) = x.coeff
add_to_coeff!(a::NCAdd, x::Number) = a.coeff += x
function set_coeff!(a::NCAdd, x::Number)
    a.coeff = x
    return a
end
function set_coeff!!(a::NCAdd{C}, x::Number) where {C}
    # decide up front, like scale!!, so that set_coeff! is only used where it can't fail
    promote_type(typeof(x), C) <: C ? set_coeff!(a, x) : NCAdd(x, a.dict)
end
# anyadd converts an NCAdd with any key type to an NCAdd whose keys hold their factors in a Vector{Any}.
# Every key type promotes to that one, so it is a closed type for KrylovKit.
function anyadd(x::NCAdd{C}) where {C}
    d = Dict{NCMul{Int,Any,Vector{Any}},C}()
    for (k, v) in x.dict
        d[NCMul(1, Vector{Any}(k.factors))] = v
    end
    NCAdd(x.coeff, d)
end
anyadd(x::NCMul) = anyadd(NCAdd(x))

const MulAdd = Union{NCMul,NCAdd}
function filter_ncadd!!(x::NCAdd; kwargs...)
    dict, coeff = filter_ncadd_dict!(x.dict; kwargs...)
    add!!(x, coeff)
end
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

Base.:+(a::Number, b::NCAdd) = NCAdd(a + additive_coeff(b), copy(b.dict))
Base.:+(a::UniformScaling, b::NCAdd) = a.λ + b
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

# The new value of the term with key `key` when `coeff` is added to it, or nothing (meaning delete the term) if the
# result is a zero term that can be dropped. For use with modify!!.
function _add_to_term(val, key, coeff)
    newval = isnothing(val) ? coeff : something(val) + coeff
    return dropzero(key, newval) ? nothing : newval
end

function add!!(_a::NCAdd, b::NCMul, α::Number=One(), β::Number=One())
    # compute β * a + α * b
    # a scalar product belongs in the coefficient, as in the NCAdd constructor
    isscalar(b) && return add!!(_a, prefactor(b), α, β)
    a = scale!!(_a, β)
    key = term_key(b)
    coeff = α * prefactor(b)
    newdict, _ = modify!!(val -> _add_to_term(val, key, coeff), a.dict, key)
    newcoeff = additive_coeff(a)
    if newdict === a.dict
        return set_coeff!!(a, newcoeff)
    end
    return NCAdd(newcoeff, newdict)
end
function add!!(_a::NCAdd, b::NCAdd, α::Number=One(), β::Number=One())
    # scaling _a in place would also scale b if they alias, e.g. in scale!!(x, x, α)
    b = _a === b ? copy(b) : b
    a = scale!!(_a, β)
    newdict = a.dict
    for (k, v) in b.dict
        coeff = α * v
        newdict, _ = modify!!(val -> _add_to_term(val, k, coeff), newdict, k)
    end
    newcoeff = additive_coeff(a) + additive_coeff(b) * α
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

# Convert, then mutate: sets `y` to the coefficient `coeff` plus the terms `terms`, an iterator of key => value pairs
# that may read from `y.dict`. Everything is converted to the types of `y` before `y` is touched, so that a failed
# conversion (e.g. an InexactError) leaves `y` unchanged. Zero terms are dropped as in the NCAdd constructor.
function _set_ncadd!(y::NCAdd{C,K}, coeff, terms) where {C,K}
    newcoeff = convert(C, coeff)
    newterms = [convert(K, k) => convert(C, v) for (k, v) in terms]
    y.coeff = newcoeff
    empty!(y.dict)
    for (k, v) in newterms
        dropzero(k, v) || (y.dict[k] = v)
    end
    return y
end
_set_ncadd!(y::NCAdd, x::NCAdd) = _set_ncadd!(y, additive_coeff(x), x.dict)
_set_ncadd!(y::NCAdd, x::NCMul) = _set_ncadd!(y, NCAdd(x))

# add! and scale! either update all of `a` or, if the result doesn't fit in the coefficient type, throw and leave
# `a` unchanged
function add!(a::NCAdd{C}, b::Number, α::Number=One(), β::Number=One()) where {C}
    # computes β * a + α * b
    newcoeff = additive_coeff(a) * β + α * b
    promote_type(typeof(β), C) <: C || return _set_ncadd!(a, newcoeff, (k => v * β for (k, v) in a.dict))
    # the scaled terms fit in C, so only the coefficient can fail to convert. Convert it first and scale in place.
    set_coeff!(a, convert(C, newcoeff))
    β isa One && return a
    map!(v -> v * β, values(a.dict))
    filter_zeros!(a.dict)
    return a
end
scale!(x::NCAdd, α::Number) = add!(x, false, One(), α)

function scale!!(x::NCAdd{CS}, α::C) where {CS,C<:Number}
    # decide up front, so that scale! is only used where it can't fail
    if promote_type(C, CS) <: CS
        scale!(x, α)
    else
        scale(x, α)
    end
end
function scale!!(y::NCAdd, x::MulAdd, α::C) where C<:Number
    return add!!(y, x, α, VectorInterface.Zero())
end
scale(x::NCAdd, α::Number) = α * x

function NCterms(a::NCAdd)
    (v * k for (k, v) in pairs(a.dict))
end

function Base.:*(x::Number, a::NCAdd{C,K}) where {C,K}
    dict = copy_dict(a, K, promote_type(typeof(x), C))
    map!(v -> x * v, values(dict))
    return NCAdd(x * additive_coeff(a), dict)
end
Base.:*(a::NCAdd, x::Number) = x * a

function Base.:*(a::NCAdd, b::NCMul)
    c = zero(a)
    filter_ncadd!!(mul!!(c, a, b); filter_zeros=true, filter_scalars=true)
end
function Base.:*(a::NCMul, b::NCAdd)
    c = zero(b)
    filter_ncadd!!(mul!!(c, a, b); filter_zeros=true, filter_scalars=true)
end
function Base.:*(a::NCAdd, b::NCAdd)
    c = zero(a)
    filter_ncadd!!(mul!!(c, a, b); filter_zeros=true, filter_scalars=true)
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
    # adjoint the terms unsorted, so each term stays an NCMul, then sort the sum once
    newx = with(_autosort => false) do
        _adjoint_terms(x)
    end
    autosort() ? sort!(newx) : newx
end
function _adjoint_terms(x::NCAdd)
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
