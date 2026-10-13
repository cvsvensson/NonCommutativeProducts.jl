function filter_ncadd_dict!(d::AbstractDict{K,V}; filter_zeros=true, filter_scalars=true) where {K<:NCMul,V}
    coeff = zero(V)
    !filter_zeros && !filter_scalars && return d, coeff
    for (k, v) in d
        if filter_zeros && iszero(v)
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

filter_zeros!(d::AbstractDict) = filter!(kv -> !iszero(last(kv)), d)

# The terms are always a Dict{K,C}: they share the coefficient type of the sum, and Dict is the only container,
# since equality and key lookup assume keys are compared by value.
# `zero_grades` holds the grades of a zero sum (see grading.jl), and is `nothing` for a nonzero sum and for a zero
# without grades (∅). It is not part of the value, so equality and hashing ignore it.
mutable struct NCAdd{C,K}
    coeff::C
    dict::Dict{K,C}
    zero_grades::Union{Nothing,Set{Any}}
    # The coefficient type follows arithmetic, so absorbing scalar terms may widen it (false + false isa Int)
    function NCAdd(coeff::C, dict::D; zero_grades=nothing, kwargs...) where {C,D<:AbstractDict}
        _, addcoeff = filter_ncadd_dict!(dict; kwargs...)
        newcoeff = coeff + addcoeff
        T = promote_type(typeof(newcoeff), valtype(D))
        newdict = dict isa Dict{keytype(D),T} ? dict : Dict{keytype(D),T}(dict)
        new{T,keytype(D)}(newcoeff, newdict, _as_grades(zero_grades))
    end
    # Exactly NCAdd{C,K}, or an InexactError if the coefficients don't fit in C
    function NCAdd{C,K}(coeff, dict::AbstractDict; zero_grades=nothing, kwargs...) where {C,K}
        newdict = dict isa Dict{K,C} ? dict : Dict{K,C}(dict)
        _, addcoeff = filter_ncadd_dict!(newdict; kwargs...)
        new{C,K}(convert(C, coeff + addcoeff), newdict, _as_grades(zero_grades))
    end
end
NCAdd{C,K}(ncadd::NCAdd{C,K}) where {C,K} = ncadd
NCAdd(ncmul::NCMul{C}) where C = _with_zero_grades(NCAdd(zero(C), to_add_dict(ncmul)), () -> _grades(ncmul))

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
anyadd(x::Union{NCMul{C},NCAdd{C}}) where {C} = convert(NCAdd{C,NCMul{Int,Any,Vector{Any}}}, x)

const MulAdd = Union{NCMul,NCAdd}
# Like the other helpers, this doesn't maintain the grades of a zero result, which is left to the public operation that
# calls it (see _with_zero_grades).
function filter_ncadd!!(x::NCAdd; kwargs...)
    dict, coeff = filter_ncadd_dict!(x.dict; kwargs...)
    _add_scaled!!(x, coeff, One())
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
Base.copy(x::NCAdd) = NCAdd(copy(additive_coeff(x)), copy(x.dict); zero_grades=x.zero_grades)

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

Base.:+(a::Number, b::NCAdd) = _with_zero_grades(NCAdd(a + additive_coeff(b), copy(b.dict)), () -> _union_grades(_grades(a), _grades(b)))
Base.:+(a::UniformScaling, b::NCAdd) = a.λ + b
Base.:+(a::NCAdd, b::B) where B<:Union{Number,UniformScaling} = b + a
Base.:+(a::NCAdd, b::B) where B<:NCMul = b + a
Base.:/(a::MulAdd, b::Number) = inv(b) * a
Base.:-(a::Union{Number,UniformScaling}, b::MulAdd) = a + (-b)
Base.:-(a::MulAdd, b::Union{Number,MulAdd,UniformScaling}) = a + (-b)
Base.:-(a::NCAdd) = (-1) * a
function Base.:+(a::NCAdd{C1,K1}, b::NCAdd{C2,K2}) where {C1,K1,C2,K2}
    dict = copy_dict(a, promote_type(K1, K2), promote_type(C1, C2))
    for (k, v) in b.dict
        dict = add_term!!(dict, k, v)
    end
    sum = NCAdd(additive_coeff(a) + additive_coeff(b), dict)
    return _with_zero_grades(sum, () -> _union_grades(_grades(a), _grades(b)))
end
add!!(a::NCMul, b::MulAdd, α::Number=One(), β::Number=One()) = add!!(a + 0, b, α, β)

# Add `coeff` to the term with key `key` in `dict`, and delete the term if it becomes zero.
# Mutates `dict` if the result fits, and otherwise returns a widened copy (e.g. true + true isa Int).
function add_term!!(dict::AbstractDict, key, coeff)
    newdict, _ = modify!!(dict, key) do val
        newval = isnothing(val) ? coeff : something(val) + coeff
        iszero(newval) ? nothing : newval
    end
    return newdict
end

# Computes β * a + α * b. A zero result has the grades of β * a and b. Those of a nonzero β * a needn't be computed,
# since its terms cancel with terms of b, and so its grades are among those of b. The grades of the destination are
# dropped when β isa Zero, since it is then ignored.
function add!!(_a::NCAdd, b::Union{MulAdd,Number,UniformScaling}, α::Number=One(), β::Number=One())
    # scaling _a in place would also scale b if they alias, e.g. in scale!!(x, x, α)
    b = _a === b ? copy(b) : b
    a = scale!!(_a, β)
    agrades = β isa VectorInterface.Zero ? nothing : _zero_grades(a)
    return _with_zero_grades(_add_scaled!!(a, b, α), () -> _union_grades(agrades, _grades(b)))
end
# a + α * b without the grades of a zero result. Loops that build a sum term by term call this directly, and leave the
# grades to the operation that calls them.
function _add_scaled!!(a::NCAdd, b::NCMul, α)
    # a scalar product belongs in the coefficient, as in the NCAdd constructor
    isscalar(b) && return _add_scaled!!(a, prefactor(b), α)
    newdict = add_term!!(a.dict, term_key(b), α * prefactor(b))
    newcoeff = additive_coeff(a)
    if newdict === a.dict
        return set_coeff!!(a, newcoeff)
    end
    return NCAdd(newcoeff, newdict)
end
function _add_scaled!!(a::NCAdd, b::NCAdd, α)
    newdict = a.dict
    for (k, v) in b.dict
        newdict = add_term!!(newdict, k, α * v)
    end
    newcoeff = additive_coeff(a) + additive_coeff(b) * α
    if newdict === a.dict
        return set_coeff!!(a, newcoeff)
    end
    return NCAdd(newcoeff, newdict)
end
_add_scaled!!(a::NCAdd, b::Number, α) = set_coeff!!(a, additive_coeff(a) + α * b)
_add_scaled!!(a::NCAdd, b::UniformScaling, α) = set_coeff!!(a, additive_coeff(a) + α * b.λ)

# Convert, then mutate: sets `y` to the coefficient `coeff` plus the terms `terms`, an iterator of key => value pairs
# that may read from `y.dict`. Everything is converted to the types of `y` before `y` is touched, so that a failed
# conversion (e.g. an InexactError) leaves `y` unchanged. Zero terms are dropped as in the NCAdd constructor.
function _set_ncadd!(y::NCAdd{C,K}, coeff, terms) where {C,K}
    newcoeff = convert(C, coeff)
    newterms = [convert(K, k) => convert(C, v) for (k, v) in terms]
    y.coeff = newcoeff
    empty!(y.dict)
    for (k, v) in newterms
        iszero(v) || (y.dict[k] = v)
    end
    return y
end
function _set_ncadd!(y::NCAdd, x::NCAdd)
    _set_ncadd!(y, additive_coeff(x), x.dict)
    y.zero_grades = x.zero_grades
    return y
end
# Fast path: the terms of `x` already have the types of `y`, so only the coefficient can fail to convert, and the
# terms are copied directly instead of through a temporary vector. If `x` shares its terms with `y` (e.g. x === y),
# they are already in place.
function _set_ncadd!(y::NCAdd{C,K}, x::NCAdd{C,K}) where {C,K}
    newcoeff = convert(C, additive_coeff(x))
    if x.dict !== y.dict
        empty!(y.dict)
        for (k, v) in x.dict
            iszero(v) || (y.dict[k] = v)
        end
    end
    y.coeff = newcoeff
    y.zero_grades = x.zero_grades
    return y
end
_set_ncadd!(y::NCAdd, x::NCMul) = _set_ncadd!(y, NCAdd(x))

# add! and scale! either update all of `a` or, if the result doesn't fit in the coefficient type, throw and leave
# `a` unchanged. An error from the arithmetic itself (e.g. an OverflowError with Rational{Int}) may leave `a` partly
# updated, since checking for that would need a temporary copy of the terms.
function add!(a::NCAdd{C}, b::Number, α::Number=One(), β::Number=One()) where {C}
    # computes β * a + α * b. The grades of a zero result are as in add!!, except that if a is scaled by zero they
    # are those of a, which are computed before its terms are dropped.
    agrades = β isa VectorInterface.Zero ? nothing : iszero(β) ? _grades(a) : _zero_grades(a)
    newgrades = () -> _union_grades(agrades, _grades(b))
    newcoeff = additive_coeff(a) * β + α * b
    if !(promote_type(typeof(β), C) <: C)
        return _with_zero_grades(_set_ncadd!(a, newcoeff, (k => v * β for (k, v) in a.dict)), newgrades)
    end
    # the scaled terms fit in C, so only the coefficient can fail to convert. Convert it first and scale in place.
    set_coeff!(a, convert(C, newcoeff))
    if !(β isa One)
        map!(v -> v * β, values(a.dict))
        iszero(β) && filter_zeros!(a.dict)
    end
    return _with_zero_grades(a, newgrades)
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

NCterms(a::NCAdd, copy=identity) = (NCMul(v, copy(k.factors)) for (k, v) in pairs(a.dict))

function Base.:*(x::Number, a::NCAdd{C,K}) where {C,K}
    dict = copy_dict(a, K, promote_type(typeof(x), C))
    map!(v -> x * v, values(dict))
    return _with_zero_grades(NCAdd(x * additive_coeff(a), dict), () -> _grades(a))
end
Base.:*(a::NCAdd, x::Number) = x * a

# zero(typeof(a)) rather than zero(a), which would compute the grades of a
Base.:*(a::NCAdd, b::NCMul) = _mul_sums(zero(typeof(a)), a, b)
Base.:*(a::NCMul, b::NCAdd) = _mul_sums(zero(typeof(b)), a, b)
Base.:*(a::NCAdd, b::NCAdd) = _mul_sums(zero(typeof(a)), a, b)
# The product a * b, accumulated in the fresh zero `c`. `c` is not an operand, so a zero result has the grades of a
# and b only.
function _mul_sums(c::NCAdd, a::MulAdd, b::MulAdd)
    product = filter_ncadd!!(_mul!!(c, a, b); filter_zeros=true, filter_scalars=true)
    return _with_zero_grades(product, () -> _mul_grades(_grades(a), _grades(b)))
end
# Computes c + a * b. A zero result has the grades of c and the products of the grades of a and b. As in add!!, those
# of a nonzero c needn't be computed, since they are among the latter.
function mul!!(c::NCAdd, a::MulAdd, b::MulAdd)
    cgrades = _zero_grades(c)
    return _with_zero_grades(_mul!!(c, a, b), () -> _union_grades(cgrades, _mul_grades(_grades(a), _grades(b))))
end
function _mul!!(c::NCAdd, a::MulAdd, b::MulAdd)
    c = _add_coeff_terms!!(c, a, b)
    autosort() && return _add_sorted_products!!(c, a, b)
    for bterm in NCterms(b)
        for aterm in NCterms(a)
            c = _add_scaled!!(c, catenate(aterm, bterm), One())
        end
    end
    return c
end
# c + acoeff * b + a * bcoeff - acoeff * bcoeff, the part of c + a * b that involves the additive coefficients
function _add_coeff_terms!!(c::NCAdd, a::MulAdd, b::MulAdd)
    acoeff = additive_coeff(a)
    bcoeff = additive_coeff(b)
    if !iszero(acoeff)
        c = _add_scaled!!(c, b, acoeff)
    end
    if !iszero(bcoeff)
        c = _add_scaled!!(c, a, bcoeff)
        c = _add_scaled!!(c, -acoeff * bcoeff, One()) # We've double counted this term so subtract it
    end
    return c
end
# The sorted sum of c and the products of the terms of a and b. The products are sorted together with the terms of c
# as one list before they are added up, which avoids hashing every unsorted product. The result is the same as
# adding the products to c and then sorting it.
function _add_sorted_products!!(c::NCAdd{C}, a::MulAdd, b::MulAdd) where {C}
    products = [catenate(aterm, bterm) for bterm in NCterms(b) for aterm in NCterms(a)]
    isempty(products) && return bubble_sort!(c)
    terms = isempty(c.dict) ? products : vcat(products, collect(NCterms(c, copy_factors)))
    # the coefficient type of c with the products added, also if they all vanish when sorted
    return _add_scaled!!(_bubble_sort!(terms, promote_type(C, _coeff_type(products))), additive_coeff(c), One())
end
_coeff_type(::Vector{<:NCMul{C}}) where {C} = C
_coeff_type(terms::Vector{<:NCMul}) = mapreduce(t -> typeof(prefactor(t)), promote_type, terms)
add!!(c::NCMul, term) = c + term
add!!(c::Number, term) = c + term

function Base.adjoint(x::NCAdd)
    # adjoint the terms unsorted, so each term stays an NCMul, then sort the sum once
    newx = with(_autosort => false) do
        _adjoint_terms(x)
    end
    newx = autosort() ? bubble_sort!(newx) : newx
    return _with_zero_grades(newx, () -> _adjoint_grades(_grades(x)))
end
function _adjoint_terms(x::NCAdd)
    newx = zero(typeof(x))
    set_coeff!(newx, adjoint(additive_coeff(x)))
    for (f, v) in x.dict
        newx = _add_scaled!!(newx, v' * f', One())
    end
    newx
end

Base.zero(::Type{NCAdd{C,K}}) where {C,K} = NCAdd{C,K}(zero(C), Dict{K,C}())
Base.one(::Type{NCAdd{C,K}}) where {C,K} = NCAdd{C,K}(one(C), Dict{K,C}())

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
