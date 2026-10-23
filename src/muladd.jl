function NCMul(f::NCAdd)
    if length(f.dict) == 1 && iszero(additive_coeff(f))
        key, coeff = only(f.dict)
        return key * coeff
    end
    throw(ArgumentError("Cannot convert NCAdd to NCMul, since the sum is not a single product. Use `terms` to get its products, or `scalar` if it has no factors: $f"))
end

# zero(x) keeps the grades of x: the zero of a product keeps its factors, with the same type as zero(typeof(x)).
# zero of a type has no grades.
Base.zero(nc::NCMul{C,S}) where {C,S} = NCMul(zero(C), convert(Vector{S}, factors_vector(nc.factors)))
Base.zero(nc::NCAdd) = _set_zero_grades!(zero(typeof(nc)), _grades(nc))
Base.one(nc::Union{<:NCAdd,<:NCMul}) = one(typeof(nc))
function Base.:(==)(a::NCAdd, b::NCMul)
    isscalar(b) && return isempty(a.dict) && additive_coeff(a) == prefactor(b)
    iszero(additive_coeff(a)) || return false
    length(a.dict) == 1 || return false
    ncmul, coeff = only(a.dict)
    NCMul(coeff, ncmul.factors) == b
end
Base.:(==)(a::NCMul, b::NCAdd) = b == a

function Base.:+(a::NCMul{C1}, b::NCMul{C2}) where {C1,C2}
    C = promote_type(C1, C2)
    K = promote_type(to_add_dict_type(typeof(a)), to_add_dict_type(typeof(b)))
    # convert the keys first so that the Dict has a single concrete key type, e.g. for Tuple- and Vector-backed terms
    dict = add_term!!(Dict{K,C}(term_key(K, a) => prefactor(a)), term_key(K, b), prefactor(b))
    return _with_zero_grades(NCAdd(zero(C), dict), () -> _union_grades(_grades(a), _grades(b)))
end

function Base.:+(a::A, b::NCMul{C}) where {A<:Number,C}
    return _with_zero_grades(NCAdd(a, to_add_dict(promote_type(A, C), b)), () -> _union_grades(_grades(a), _grades(b)))
end
Base.:+(a::UniformScaling, b::NCMul) = a.λ + b
Base.:+(a::NCMul, b::Union{Number,UniformScaling}) = b + a
Base.:*(x::MulAdd, J::UniformScaling) = x * J.λ
Base.:*(J::UniformScaling, x::MulAdd) = J.λ * x
function Base.:+(a::NCMul{C1}, b::NCAdd{C2,K2}) where {C1,C2,K2}
    C = promote_type(C1, C2)
    K = promote_type(to_add_dict_type(typeof(a)), K2)
    newdict = add_term!!(copy_dict(b, K, C), term_key(K, a), prefactor(a))
    return _with_zero_grades(NCAdd(additive_coeff(b), newdict), () -> _union_grades(_grades(a), _grades(b)))
end

to_add_dict(a::NCMul{C}) where C = to_add_dict(C, a)
to_add_dict(::Type{T}, a::NCMul) where {T<:Number} = Dict(term_key(a) => convert(T, prefactor(a)))
# The key of the term `a` in a sum, which may share its factors with `a` (see the comment on NCMul)
term_key(a::NCMul) = term_key(to_add_dict_type(typeof(a)), a)
term_key(::Type{K}, a::NCMul) where {K} = convert(K, NCMul(1, a.factors))
# The key type of a sum with the term `a`, or of the sum `a`. Keys of Tuple-backed products are Vector-backed, so that
# a sum can hold products of any length.
to_add_dict_type(::Type{NCMul{C,S,F}}) where {C,S,F} = NCMul{Int,S,F}
to_add_dict_type(::Type{NCMul{C,S,F}}) where {C,S,F<:Tuple} = NCMul{Int,S,Vector{S}}
to_add_dict_type(::Type{NCAdd{C,K}}) where {C,K} = K
to_add_dict_type(::Type{NCMul{C,S}}) where {C,S} = NCMul{Int,S}
to_add_dict_type(::Type{NCMul{C}}) where C = NCMul{Int}
to_add_dict_type(::Type{NCMul}) = NCMul{Int}
function Base.:^(a::Union{NCAdd,NCMul}, b::Integer)
    b < 0 && throw(DomainError(b, "products and sums have no inverse, so the exponent must be non-negative"))
    ret = Base.power_by_squaring(a, b)
    # for b >= 2, ret is a product computed by `*`, which already sorts
    autosort() && b < 2 && return bubble_sort(ret)
    return ret
end
# Base rewrites a negative literal power `x^-n` as `inv(x)^n`. Forward it to `^` instead.
Base.literal_pow(::typeof(^), a::Union{NCAdd,NCMul}, ::Val{p}) where {p} = a^p

# The product `ncmul` of the atoms or products `x` and `y`, sorted if autosort is enabled. Used by @nc.
function _autosort_product(ncmul::NCMul, x, y)
    autosort() || return ncmul
    return _with_zero_grades(bubble_sort!(ncmul), () -> _mul_grades(_grades(x), _grades(y)))
end

macro nc_common(T)
    quote
        NonCommutativeProducts._is_atom(::Type{<:$(esc(T))}) = true

        NonCommutativeProducts.NCMul(f::$(esc(T))) = NCMul(1, [f])
        NonCommutativeProducts.NCAdd(f::$(esc(T))) = NCAdd(NCMul(f))
        NonCommutativeProducts.ncmapreduce(f, ops::Tuple, x::$(esc(T)); scalarmap=identity) = f(x)
        NonCommutativeProducts.terms(x::$(esc(T))) = terms(NCMul(x))
        NonCommutativeProducts.factors(x::$(esc(T))) = [x]
        NonCommutativeProducts.coefficient(x::$(esc(T))) = 1
        NonCommutativeProducts.nterms(::$(esc(T))) = 1

        Base.:+(x::$(esc(T)), y::$(esc(T))) = NCMul(x) + NCMul(y)
        Base.:+(x::$(esc(T)), y::Union{Number,UniformScaling,NCMul,NCAdd}) = NCMul(x) + y
        Base.:+(x::Union{Number,UniformScaling,NCMul,NCAdd}, y::$(esc(T))) = y + x

        NonCommutativeProducts.add!!(x::NCAdd, y::$(esc(T))) = add!!(x, NCMul(y))

        Base.:-(x::$(esc(T)), y::$(esc(T))) = NCMul(x) - NCMul(y)
        Base.:-(x::$(esc(T))) = NCMul(-1, [x])
        Base.:-(x::Union{Number,UniformScaling,NCMul,NCAdd}, y::$(esc(T))) = x - NCMul(y)
        Base.:-(x::$(esc(T)), y::Union{Number,UniformScaling,NCMul,NCAdd}) = NCMul(x) - y

        Base.:*(x::$(esc(T)), y::$(esc(T))) = _autosort_product(NCMul(1, [x, y]), x, y)
        Base.:*(x::$(esc(T)), y::NCMul) = _autosort_product(NCMul(prefactor(y), pushfirst!!(copy_factors(y.factors), x)), x, y)
        Base.:*(x::NCMul, y::$(esc(T))) = _autosort_product(NCMul(prefactor(x), push!!(copy_factors(x.factors), y)), x, y)
        Base.:*(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = x * NCMul(y)
        Base.:*(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = NCMul(x) * y
        Base.:/(x::$(esc(T)), y::Number) = NCMul(x) / y
        Base.://(x::$(esc(T)), y::Number) = NCMul(x) // y

        Base.:^(a::$(esc(T)), b) = NCMul(a)^b
        Base.literal_pow(::typeof(^), a::$(esc(T)), ::Val{p}) where {p} = NCMul(a)^p

        Base.:(==)(a::$(esc(T)), b::Union{NCMul,NCAdd}) = NCMul(a) == b
        Base.:(==)(a::Union{NCMul,NCAdd}, b::$(esc(T))) = b == a

        # the zero of an atom keeps it as a factor for its grade, with the same type as zero(T)
        Base.zero(a::$(esc(T))) = NCMul(0, typeof(a)[a])
        function Base.zero(::Type{W}) where W<:$(esc(T))
            NCMul(0, W[])
        end
        Base.one(a::$(esc(T))) = one($(esc(T)))
        function Base.one(::Type{W}) where W<:$(esc(T))
            NCMul(1, W[])
        end
        Base.oneunit(a::$(esc(T))) = oneunit($(esc(T)))
        function Base.oneunit(::Type{W}) where W<:$(esc(T))
            NCMul(1, W[])
        end

        # an atom promotes and converts as the product with itself as the only factor
        Base.promote_rule(::Type{W}, ::Type{NC}) where {W<:$(esc(T)),NC<:MulAdd} = promote_type(NCMul{Int,W,Vector{W}}, NC)
        Base.convert(::Type{NC}, x::$(esc(T))) where {NC<:MulAdd} = convert(NC, NCMul(x))

        NonCommutativeProducts.add!!(x::MulAdd, y::$(esc(T)), α::Number, β::Number) = add!!(x, NCMul(y), α, β)
        NonCommutativeProducts.add!!(x::$(esc(T)), y::$(esc(T)), α::Number, β::Number) = add!!(NCMul(x), NCMul(y), α, β)
        NonCommutativeProducts.add!!(x::$(esc(T)), y::MulAdd, α::Number, β::Number) = add!!(NCMul(x), y, α, β)

        @nc_vectorinterface $(esc(T))

        NonCommutativeProducts.anyadd(x::$(esc(T))) = anyadd(NCAdd(x))
    end
end

const _DEFAULT_AUTOSORT = Ref(false)
const _autosort = ScopedValue{Bool}()

"""
    autosort()

Whether products are sorted with `mul_effect` as they are formed. Set the global default with `enable_autosort!` and
`disable_autosort!`, or override it for a block of code with `with_autosort`. Autosorting is off by default.

Returns the current setting as a `Bool`: the value set by the innermost `with_autosort` if one is active, otherwise the global default.
"""
function autosort()
    v = ScopedValues.get(_autosort)
    return v === nothing ? _DEFAULT_AUTOSORT[] : something(v)
end
"""
    enable_autosort!()

Sort products with `mul_effect` as they are formed, unless overridden by `with_autosort`.
"""
enable_autosort!() = _DEFAULT_AUTOSORT[] = true
"""
    disable_autosort!()

Do not sort products as they are formed, unless overridden by `with_autosort`. This is the default.
"""
disable_autosort!() = _DEFAULT_AUTOSORT[] = false

"""
    with_autosort(f, enabled::Bool=true)

Call `f()` with autosorting set to `enabled`, overriding the global setting from `enable_autosort!` and `disable_autosort!` for the duration of the call. Returns the
value of `f()`. Works with a `do` block:

    with_autosort(false) do
        a * b
    end
"""
with_autosort(f, enabled::Bool=true) = with(f, _autosort => enabled)

"""
    @nc T1 T2 ...

Define addition and multiplication for the given types, both within each type and between them, and make them
implement the VectorInterface.jl API.
"""
macro nc(types...)
    nc_common_calls = [:(@nc_common $(esc(T))) for T in types]
    quote
        $(nc_common_calls...)
        @nc_pairs $(esc.(types)...)
    end

end
macro nc_pairs(types...)
    mul_pairs = Expr[]
    for T1 in types
        for T2 in types
            T1 == T2 && continue
            push!(mul_pairs, :(Base.:*(x::$(esc(T1)), y::$(esc(T2))) = _autosort_product(NCMul(1, [x, y]), x, y)))
        end
    end

    add_pairs = Expr[]
    for T1 in types
        for T2 in types
            T1 == T2 && continue
            push!(add_pairs, :(Base.:+(x::$(esc(T1)), y::$(esc(T2))) = NCMul(x) + NCMul(y)))
            push!(add_pairs, :(Base.:-(x::$(esc(T1)), y::$(esc(T2))) = NCMul(x) - NCMul(y)))
        end
    end

    quote
        $(mul_pairs...)
        $(add_pairs...)
    end

end
"""
    @commutative T1 T2 ...

Define mixed arithmetic for types that commute with each other, and order them as listed: `mul_effect(x::Ti, y::Tj)`
is `nothing` if `i < j` and `y * x` if `i > j`. The types must already be registered with `@nc`.
"""
macro commutative(types...)
    mul_effect = Expr[]
    for (n1, T1) in enumerate(types)
        for (n2, T2) in enumerate(types)
            n1 >= n2 && continue
            push!(mul_effect, :(NonCommutativeProducts.mul_effect(x::$(esc(T1)), y::$(esc(T2))) = nothing))
            push!(mul_effect, :(NonCommutativeProducts.mul_effect(x::$(esc(T2)), y::$(esc(T1))) = y * x))
        end
    end
    quote
        $(mul_effect...)
        @nc_pairs $(esc.(types)...)
    end
end
