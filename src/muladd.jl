abstract type AbstractOrdering end
struct NormalOrdering <: AbstractOrdering end
struct NaiveOrdering <: AbstractOrdering end

ordered_product(a, b, ::NaiveOrdering) = NCMul(1, [a, b])

struct NCMul{C,S}
    coeff::C
    factors::Vector{S}
    function NCMul(coeff::C, factors) where {C}
        new{C,eltype(factors)}(coeff, factors)
    end
end
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
mutable struct NCAdd{C,D}
    coeff::C
    dict::D
end
Base.:(==)(a::NCAdd, b::NCAdd) = a.coeff == b.coeff && a.dict == b.dict
Base.hash(a::NCAdd, h::UInt) = hash(a.coeff, hash(a.dict, h))
function filter_scalars!(f::NCAdd)
    for (k, v) in f.dict
        if isscalar(k)
            f.coeff += k.coeff * v
            delete!(f.dict, k)
        end
    end
    f
end
NCMul(f::NCAdd) = (length(f.dict) == 0 && iszero(f.coeff) && return prod(only(f.dict))) || throw(ArgumentError("Cannot convert NCAdd to NCMul: $f"))
function canonicalize!(_f::NCAdd{_C,D}, filter_scalars=true) where {_C,D}
    C = promote_type(_C, valtype(D))
    f = NCAdd(C(_f.coeff), _f.dict)
    filter_scalars && filter_scalars!(f)
    if length(f.dict) == 0
        return f.coeff
    elseif length(f.dict) == 1 && iszero(f.coeff)
        k, v = only(f.dict)
        return v * k
    else
        return NCAdd{C,D}(f.coeff, f.dict)
    end
end
canonicalize!(a::Number) = a
const SM = Union{NCMul}
const SMA = Union{NCMul,NCAdd}

function show_compact_sum(io, x::NCAdd, max_terms=3)
    println(io, "Sum with ", length(x.dict), " terms: ")
    N = min(max_terms, length(x.dict))
    args = sum(v * k for (k, v) in Iterators.take(pairs(x.dict), N))
    show(io, args)
    if N < length(x.dict)
        print(io, " + ...")
    end
    return nothing
end
function Base.show(io::IO, x::NCAdd)
    if length(x.dict) > 3
        return show_compact_sum(io, x)
    end
    compact = get(io, :compact, false)
    args = arguments(x)
    print_one = !iszero(x.coeff)
    if print_one
        if isreal(x.coeff)
            print(io, real(x.coeff), "I")
        else
            print(io, "(", x.coeff, ")", "I")
        end
        args = args[2:end]
    end
    print_sign(s) = compact ? print(io, s) : print(io, " ", s, " ")
    for (n, arg) in enumerate(args)
        k = prod(arg.factors)
        v = arg.coeff
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
    return nothing
end
print_num(io::IO, x) = isreal(x) ? print(io, real(x)) : print(io, "(", x, ")")

Base.:+(a::Number, b::SM) = iszero(a) ? b : NCAdd(a, to_add(b))
Base.:+(a::UniformScaling, b::SM) = iszero(a) ? b : NCAdd(a.λ, to_add(b))
Base.:+(a::SM, b::Union{Number,UniformScaling}) = b + a
Base.:+(a::NCMul, b) = a + NCMul(b)
Base.:+(a, b::NCMul) = 1 * a + b
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
Base.:+(a::SM, b::NCAdd) = canonicalize!(NCAdd(b.coeff, (_merge(+, to_add(a), b.dict; filter=iszero))))
function add!(a::NCAdd, b::NCAdd)
    a.coeff += b.coeff
    a.dict = __merge!(+, a.dict, b.dict; filter=iszero)
    return a
end
function add!(a::NCAdd, b::SM)
    a.dict = __merge!(+, a.dict, to_add_tuple(b); filter=iszero)
    return a
end
add!(a::NCAdd, b::Number) = (a.coeff += b; return a)
add!(a::NCAdd, b::UniformScaling) = (a.coeff += b.λ; return a)
to_add(a::NCMul, coeff=1) = Dict(NCMul(1, a.factors) => a.coeff * coeff)
to_add(a, coeff=1) = Dict(NCMul(a) => coeff)
to_add_tuple(a::NCMul, coeff=1) = (NCMul(1, a.factors) => a.coeff * coeff,)
to_add_tuple(a, coeff=1) = (NCMul(a) => coeff,)

Base.:+(a::Number, b::NCAdd) = iszero(a) ? b : NCAdd(a + b.coeff, b.dict)
Base.:+(a::UniformScaling, b::NCAdd) = iszero(a) ? b : NCAdd(a.λ + b.coeff, b.dict)
Base.:+(a::NCAdd, b::Union{Number,SM,UniformScaling}) = b + a
Base.:/(a::SMA, b::Number) = inv(b) * a
Base.:-(a::Union{Number,UniformScaling}, b::SMA) = a + (-b)
Base.:-(a::SMA, b::Union{Number,SMA,UniformScaling}) = a + (-b)
Base.:-(a::SMA) = -1 * a
function NCterms(a::NCAdd)
    (v * k for (k, v) in pairs(a.dict))
end
NCterms(a::NCMul) = (a,)
function allterms(a::NCAdd)
    [a.coeff, [v * k for (k, v) in pairs(a.dict)]...]
end
function Base.:+(a::NCAdd, b::NCAdd)
    coeff = a.coeff + b.coeff
    dict = _merge(+, a.dict, b.dict; filter=iszero)
    NCAdd(coeff, dict)
end
Base.:^(a::Union{NCMul,NCAdd}, b) = Base.power_by_squaring(a, b)

ordered_product(x::Number, y::Number, ordering) = x * y
ordered_product(x::Number, a::NCMul, ordering) = NCMul(x * a.coeff, a.factors)
ordered_product(x::Number, a::NCAdd, ordering) = NCAdd(x * a.coeff, Dict(k => v * x for (k, v) in a.dict))
ordered_product(a::NCMul, b::NCMul, ::NaiveOrdering) = NCMul(a.coeff * b.coeff, vcat(a.factors, b.factors))
ordered_product(a::SMA, x::Number, ordering) = ordered_product(x, a, ordering)

ordered_product(bs::NCMul, ordering) = ordered_product(NCMul(a), bs, ordering)
ordered_product(as::NCMul, b, ordering) = ordered_product(as, NCMul(b), ordering)
ordered_product(as::NCMul, bs::NCMul, ::NormalOrdering) = canonicalize!(normal_order(ordered_product(as, bs, NaiveOrdering())))
additive_coeff(a::NCAdd) = a.coeff
additive_coeff(a::SM) = 0

Base.:*(a::SMA, b::SMA) = ordered_product(a, b, NormalOrdering())
Base.:*(a::SMA, x::Number) = ordered_product(x, a, NaiveOrdering())
Base.:*(x::Number, a::SMA) = ordered_product(x, a, NaiveOrdering())
function ordered_product(a::NCAdd, b::SM, ordering::AbstractOrdering)
    c = zero(a)
    return trymul!(c, a, b, ordering)
end
function ordered_product(a::SM, b::NCAdd, ordering)
    c = zero(b)
    return trymul!(c, a, b, ordering)
end
function ordered_product(a::NCAdd, b::NCAdd, ordering)
    c = zero(a)
    return trymul!(c, a, b, ordering)
end
function trymul!(c::NCAdd, a::SMA, b::SMA, ordering)
    acoeff = additive_coeff(a)
    bcoeff = additive_coeff(b)
    if !iszero(acoeff)
        c = tryadd!(c, ordered_product(acoeff, b, ordering))
    end
    if !iszero(bcoeff)
        c = tryadd!(c, ordered_product(a, bcoeff, ordering))
        c = tryadd!(c, -acoeff * bcoeff) # We've double counted this term so subtract it
    end
    for bterm in NCterms(b)
        for aterm in NCterms(a)
            newterm = ordered_product(aterm, bterm, ordering)
            c = tryadd!(c, newterm)
        end
    end
    return canonicalize!(c)
end

"""
    tryadd!(c::NCAdd, term)

This function tries to add `term` to `c` in place. If it fails, it catches the error and returns the out of place sum.
"""
function tryadd!(c::NCAdd, term)
    try
        return add!(c, term)
    catch e
        @debug e
    end
    c + term
end
tryadd!(c::SM, term) = c + term
tryadd!(c::Number, term) = c + term


Base.adjoint(x::NCMul) = length(x.factors) == 0 ? NCMul(adjoint(x.coeff), x.factors) : adjoint(x.coeff) * foldr(*, Iterators.reverse(Iterators.map(adjoint, x.factors)))
function Base.adjoint(x::NCAdd)
    newx = zero(x)
    newx.coeff = adjoint(x.coeff)
    for (f, v) in x.dict
        add!(newx, v' * f')
    end
    newx
end
Base.zero(::NCAdd{C,D}) where {C,D} = NCAdd(zero(C), D())

function sorted_noduplicates(v)
    I = eachindex(v)
    for i in I[1:end-1]
        isequal(v[i], v[i+1]) && return false
    end
    return true
end

isscalar(x::NCMul) = iszero(x.coeff) || (length(x.factors) == 0)
isscalar(x::NCAdd) = length(x.dict) == 0 || all(isscalar, keys(x.dict)) || all(iszero(values(x.dict)))


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


#From SymbolicUtils
function _merge(f, d, others...; filter=x -> false)
    T = Union{promote_type(valtype(d), valtype.(others)...)}
    K = Union{keytype(d),keytype.(others)...}
    d2 = Dict{K,T}(d)
    return __merge!(f, d2, others...; filter)
end
function __merge!(f::F, d, others...; filter=x -> false) where {F}
    acc = d
    for other in others
        for (k, v) in other
            v = f(v)
            ak = get(acc, k, nothing)
            if ak !== nothing
                v = ak + v
            end
            if filter(v)
                delete!(acc, k)
            else
                acc[k] = v
            end
        end
    end
    acc
end


Base.copy(x::NCAdd) = NCAdd(copy(x.coeff), copy(x.dict))
@testitem "Consistency between + and add!" begin
    import NonCommutativeProducts: add!
    @fermions f
    a = 1.0 * f[2] * f[1] + 1 + f[1]
    for b in [1.0, 1, f[1], 1.0 * f[1], f[2] * f[1], a]
        a2 = copy(a)
        a3 = add!(a2, b)
        @test a + b == a3
        @test a2 == a3
        @test_throws InexactError add!(a, 1im * b)
    end
    @test a == 1.0 * f[2] * f[1] + 1 + f[1]
end
