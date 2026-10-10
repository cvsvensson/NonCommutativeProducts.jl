NCMul(f::NCAdd) = (length(f.dict) == 1 && iszero(additive_coeff(f)) && return prod(only(f.dict))) || throw(ArgumentError("Cannot convert NCAdd to NCMul: $f"))

Base.zero(nc::Union{<:NCAdd,<:NCMul}) = zero(typeof(nc))
Base.one(nc::Union{<:NCAdd,<:NCMul}) = one(typeof(nc))
function Base.:(==)(a::NCAdd, b::NCMul)
    if isscalar(b)
        additive_coeff(a) == prefactor(b) && length(a.dict) == 0 && return true
    end
    iszero(additive_coeff(a)) || return false
    length(a.dict) == 1 || return false
    ncmul, coeff = only(a.dict)
    NCMul(coeff, ncmul.factors) == b
end
Base.:(==)(a::NCMul, b::NCAdd) = b == a

function Base.:+(a::NCMul{C1}, b::NCMul{C2}) where {C1,C2}
    C = promote_type(C1, C2)
    K = promote_type(to_add_dict_type(typeof(a)), to_add_dict_type(typeof(b)))
    # convert the keys up front so that the Dict has a single concrete key type, e.g. when one term is Tuple-backed and the other Vector-backed
    dict = add_term!!(Dict{K,C}(term_key(K, a) => prefactor(a)), term_key(K, b), prefactor(b))
    return NCAdd(zero(C), dict)
end

Base.:+(a::A, b::NCMul{C}) where {A<:Number,C} = NCAdd(a, to_add_dict(promote_type(A, C), b))
Base.:+(a::UniformScaling, b::NCMul) = a.λ + b
Base.:+(a::NCMul, b::Union{Number,UniformScaling}) = b + a
function Base.:+(a::NCMul{C1}, b::NCAdd{C2,K2}) where {C1,C2,K2}
    C = promote_type(C1, C2)
    K = promote_type(to_add_dict_type(typeof(a)), K2)
    newdict = add_term!!(copy_dict(b, K, C), term_key(K, a), prefactor(a))
    return NCAdd(additive_coeff(b), newdict)
end

to_add_dict(a::NCMul{C}) where C = to_add_dict(C, a)
to_add_dict(::Type{T}, a::NCMul) where {T<:Number} = Dict(term_key(a) => convert(T, prefactor(a)))
# The key of the term `a` in a sum. It may share its factors with `a`, since factors are never mutated once they
# are part of an NCMul that escapes the package (see the comment on NCMul).
term_key(a::NCMul) = term_key(to_add_dict_type(typeof(a)), a)
term_key(::Type{K}, a::NCMul) where {K} = convert(K, NCMul(1, a.factors))
# The key type of a sum with the term `a`, or of the sum `a`. Every key of a sum built from products is derived
# from this. Like zero and one (see mul.jl), keys built from Tuple-backed products are Vector-backed, so that the
# sum can hold products of any length.
to_add_dict_type(::Type{NCMul{C,S,F}}) where {C,S,F} = NCMul{Int,S,F}
to_add_dict_type(::Type{NCMul{C,S,F}}) where {C,S,F<:Tuple} = NCMul{Int,S,Vector{S}}
to_add_dict_type(::Type{NCAdd{C,K}}) where {C,K} = K
to_add_dict_type(::Type{NCMul{C,S}}) where {C,S} = NCMul{Int,S}
to_add_dict_type(::Type{NCMul{C}}) where C = NCMul{Int}
to_add_dict_type(::Type{NCMul}) = NCMul{Int}
function Base.:^(a::Union{NCAdd,NCMul}, b::Int)
    ret = Base.power_by_squaring(a, b)
    autosort() && return bubble_sort!(ret)
    return ret
end

macro nc_common(T)
    quote
        NonCommutativeProducts.NCMul(f::$(esc(T))) = NCMul(1, [f])
        NonCommutativeProducts.NCAdd(f::$(esc(T))) = NCAdd(NCMul(f))
        NonCommutativeProducts.ncmapreduce(f, ops::Tuple, x::$(esc(T)); scalarmap=identity) = f(x)

        Base.:+(x::$(esc(T)), y::$(esc(T))) = NCMul(x) + NCMul(y)
        Base.:+(x::$(esc(T)), y::Union{Number,UniformScaling,NCMul,NCAdd}) = NCMul(x) + y
        Base.:+(x::Union{Number,UniformScaling,NCMul,NCAdd}, y::$(esc(T))) = y + x

        NonCommutativeProducts.add!!(x::NCAdd, y::$(esc(T))) = add!!(x, NCMul(y))

        Base.:-(x::$(esc(T)), y::$(esc(T))) = NCMul(x) - NCMul(y)
        Base.:-(x::$(esc(T))) = NCMul(-1, [x])
        Base.:-(x::Union{Number,UniformScaling,NCMul,NCAdd}, y::$(esc(T))) = x - NCMul(y)
        Base.:-(x::$(esc(T)), y::Union{Number,UniformScaling,NCMul,NCAdd}) = NCMul(x) - y

        Base.:*(x::$(esc(T)), y::$(esc(T))) = autosort() ? bubble_sort!(NCMul(1, [x, y])) : NCMul(1, [x, y])
        function Base.:*(x::$(esc(T)), y::NCMul)
            ncmul = NCMul(prefactor(y), pushfirst!!(copy_factors(y.factors), x))
            autosort() ? bubble_sort!(ncmul) : ncmul
        end
        function Base.:*(x::NCMul, y::$(esc(T)))
            ncmul = NCMul(prefactor(x), push!!(copy_factors(x.factors), y))
            autosort() ? bubble_sort!(ncmul) : ncmul
        end
        Base.:*(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = x * NCMul(y)
        Base.:*(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = NCMul(x) * y
        Base.:/(x::$(esc(T)), y::Number) = NCMul(x) / y

        Base.:^(a::$(esc(T)), b) = NCMul(a)^b

        Base.:(==)(a::$(esc(T)), b::Union{NCMul,NCAdd}) = NCMul(a) == b
        Base.:(==)(a::Union{NCMul,NCAdd}, b::$(esc(T))) = b == a

        Base.zero(a::$(esc(T))) = zero($(esc(T)))
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

        VectorInterface.inner(x::MulAdd, y::$(esc(T))) = _inner(x, y)
        VectorInterface.inner(x::$(esc(T)), y::MulAdd) = _inner(x, y)
        VectorInterface.inner(x::$(esc(T)), y::$(esc(T))) = _inner(x, y)
        LinearAlgebra.norm(x::$(esc(T))) = sqrt(VectorInterface.inner(x, x))

        NonCommutativeProducts.add!!(x::MulAdd, y::$(esc(T)), α::Number, β::Number) = add!!(x, NCMul(y), α, β)
        NonCommutativeProducts.add!!(x::$(esc(T)), y::$(esc(T)), α::Number, β::Number) = add!!(NCMul(x), NCMul(y), α, β)
        NonCommutativeProducts.add!!(x::$(esc(T)), y::MulAdd, α::Number, β::Number) = add!!(NCMul(x), y, α, β)

        VectorInterface.scale(x::$(esc(T)), α::Number) = α * x
        VectorInterface.scale!!(x::$(esc(T)), α::Number) = α * x
        VectorInterface.scale!!(a::NCAdd, x::$(esc(T)), α::Number) = add!!(a, NCMul(x), α, VectorInterface.Zero())
        VectorInterface.zerovector(x::$(esc(T)), ::Type{S}) where {S<:Number} = VectorInterface.zerovector(NCMul(x), S)

        NonCommutativeProducts.anyadd(x::$(esc(T))) = anyadd(NCAdd(x))
    end
end

isfilterable(x) = true
const _DEFAULT_AUTOSORT = Ref(false)
const _autosort = ScopedValue{Bool}()

function autosort()
    Base.isassigned(_autosort) && return _autosort[]
    return _DEFAULT_AUTOSORT[]
end
enable_autosort!() = _DEFAULT_AUTOSORT[] = true
disable_autosort!() = _DEFAULT_AUTOSORT[] = false


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
            push!(mul_pairs, :(Base.:*(x::$(esc(T1)), y::$(esc(T2))) = autosort() ? bubble_sort!(NCMul(1, [x, y])) : NCMul(1, [x, y])))
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

@testitem "Addition and multiplication with different number types" setup = [Fermions] begin
    # This tests that addition and multiplication with different number types does not throw errors due to promotion issues and in place operations
    f1 = Fermion(1)
    f2 = Fermion(2)
    nums = (1, 1/2, 1//3, 0.1+1im, big(4))
    for x in nums
        @test x + f1 == f1 + x
        for y in nums
            @test x + y*f1 == y*f1 + x
            for z in nums
                @test z*(x+y*f1) isa NonCommutativeProducts.NCAdd #throws error on v0.4.4
            end
        end
    end
end
