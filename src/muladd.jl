NCMul(f::NCAdd) = (length(f.dict) == 1 && iszero(additive_coeff(f)) && return prod(only(f.dict))) || throw(ArgumentError("Cannot convert NCAdd to NCMul: $f"))

Base.zero(nc::Union{<:NCAdd,<:NCMul}) = zero(typeof(nc))
Base.one(nc::Union{<:NCAdd,<:NCMul}) = one(typeof(nc))
function Base.promote_rule(::Type{<:NCAdd{C1,NCMul{Int,S,VS},D1}}, ::Type{<:NCAdd{C2,NCMul{Int,S,VS},D2}}) where {C1,C2,D1,D2,S,VS}
    C = promote_type(C1, C2)
    D = promote_type(D1, D2)
    NCMUL = NCMul{Int,S,VS}
    return NCAdd{C,NCMUL,D}
end
Base.promote_rule(::Type{<:NCMul{C1,S1,VS1}}, ::Type{<:NCMul{C2,S2,VS2}}) where {C1,C2,S1,S2,VS1,VS2} = NCMul{promote_type(C1, C2),promote_type(S1, S2),promote_type(VS1, VS2)}

function Base.promote_rule(::Type{<:NCMul{C1,S,VS1}}, ::Type{<:NCAdd{C2,NCMul{Int,S,VS2},D}}) where {C1,C2,S,VS1,VS2,D}
    C = promote_type(C1, C2)
    VS = promote_type(VS1, VS2)
    NCMUL = NCMul{Int,S,VS}
    DD = promote_type(Dict{NCMUL,C}, D)
    return NCAdd{C,NCMUL,DD}
end
function Base.promote_rule(::Type{A}, ::Type{M}) where {A<:NCAdd,M<:NCMul}
    promote_rule(M, A)
end

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

function Base.:+(a::NCMul{C1,F1,S}, b::NCMul{C2,F2,S}) where {C1,C2,F1,F2,S}
    if a.factors == b.factors
        return NCAdd(0, Dict(NCMul(1, a.factors) => prefactor(a) + prefactor(b)))
    end
    F = Union{F1,F2}
    C = promote_type(C1, C2)
    return NCAdd(0, Dict{NCMul{Int,F,S},C}(NCMul{Int,F,S}(1, a.factors) => prefactor(a), NCMul{Int,F,S}(1, b.factors) => prefactor(b)))
end
function Base.:+(a::NCMul{C1,F1,S1}, b::NCMul{C2,F2,S2}) where {C1,C2,F1,F2,S1,S2}
    if a.factors == b.factors
        return NCAdd(0, Dict(NCMul(1, a.factors) => prefactor(a) + prefactor(b)))
    end
    F = promote_type(F1, F2)
    C = promote_type(C1, C2)
    S = promote_type(S1, S2)
    return NCAdd(0, Dict{NCMul{Int,F,S},C}(NCMul{Int,F,S}(1, a.factors) => prefactor(a), NCMul{Int,F,S}(1, b.factors) => prefactor(b)))
end

Base.:+(a::Number, b::NCMul) = NCAdd(a, to_add_dict(b))
Base.:+(a::UniformScaling, b::NCMul) = NCAdd(a.λ, to_add_dict(b))
Base.:+(a::NCMul, b::Union{Number,UniformScaling}) = b + a
function Base.:+(a::NCMul, b::NCAdd)
    newdict = copy(b.dict)
    for (k, v) in b.dict
        if k.factors == a.factors
            newdict[k] = v + prefactor(a)
            return NCAdd(additive_coeff(b), newdict)
        end
    end
    newdict2 = setindex!!(newdict, prefactor(a), NCMul(1, a.factors))
    nc = NCAdd(additive_coeff(b), newdict2)
    return nc
end
function Base.convert(::Type{NCAdd{C,NCMul{Int,S,F},_D}}, x::NCMul{C2,S,F}) where {C,C2,S,F,_D}
    D = Dict{NCMul{Int,S,F},C2}
    NCAdd(zero(C), D(to_add_dict(x)))
end

to_add_dict(a::NCMul) = Dict(NCMul(1, a.factors) => prefactor(a))
to_add_dict_type(::Type{NCMul{C,S,F}}) where {C,S,F} = NCMul{Int,S,F}
to_add_dict_type(::Type{NCMul}) = NCMul{Int}
function Base.:^(a::Union{NCAdd,NCMul}, b::Int)
    ret = Base.power_by_squaring(a, b)
    autosort() && return sort!(ret)
    return ret
end

macro nc_common(T)
    quote
        NonCommutativeProducts.NCMul(f::$(esc(T))) = NCMul(1, [f])
        NonCommutativeProducts.ncmapreduce(f, ops::Tuple, x::$(esc(T)); scalarmap=identity) = f(x)

        Base.:+(x::$(esc(T)), y::$(esc(T))) = NCMul(x) + NCMul(y)
        Base.:+(x::$(esc(T)), y::Union{Number,UniformScaling,NCMul,NCAdd}) = NCMul(x) + y
        Base.:+(x::Union{Number,UniformScaling,NCMul,NCAdd}, y::$(esc(T))) = y + x

        NonCommutativeProducts.add!!(x::NCAdd, y::$(esc(T))) = add!!(x, NCMul(y))

        Base.:-(x::$(esc(T)), y::$(esc(T))) = NCMul(x) - NCMul(y)
        Base.:-(x::$(esc(T))) = NCMul(-1, [x])
        Base.:-(x::Union{Number,UniformScaling,NCMul,NCAdd}, y::$(esc(T))) = x - NCMul(y)
        Base.:-(x::$(esc(T)), y::Union{Number,UniformScaling,NCMul,NCAdd}) = NCMul(x) - y

        Base.:*(x::$(esc(T)), y::$(esc(T))) = autosort() ? sort!(NCMul(1, [x, y])) : NCMul(1, [x, y])
        function Base.:*(x::$(esc(T)), y::NCMul)
            ncmul = NCMul(prefactor(y), pushfirst!!(copy(y.factors), x))
            autosort() ? sort!(ncmul) : ncmul
        end
        function Base.:*(x::NCMul, y::$(esc(T)))
            ncmul = NCMul(prefactor(x), push!!(copy(x.factors), y))
            autosort() ? sort!(ncmul) : ncmul
        end
        Base.:*(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = x * NCMul(y)
        Base.:*(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = NCMul(x) * y
        Base.:/(x::$(esc(T)), y::Number) = NCMul(x) / y

        Base.:^(a::$(esc(T)), b) = NCMul(a)^b
        Base.convert(::Type{NCMul{C,S,F}}, x::$(esc(T))) where {C,S<:$(esc(T)),F} = NCMul(one(C), S[x])

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

        Base.promote_rule(::Type{W}, ::Type{W}) where W<:$(esc(T)) = return W
        function Base.promote_rule(::Type{W}, ::Type{<:NCAdd{C,NCMul{Int,W,VS},D}}) where {C,VS,D,W<:$(esc(T))}
            return NCAdd{C,NCMul{Int,W,VS},D}
        end
        function Base.promote_rule(::Type{W}, ::Type{<:NCMul{C,W,VS}}) where {C,VS,W<:$(esc(T))}
            return NCMul{C,W,VS}
        end
        function Base.promote_rule(::Type{NC}, ::Type{W}) where {NC<:MulAdd,W<:$(esc(T))}
            return promote_rule(W, NC)
        end

        Base.convert(::Type{NCMul{C,W,V}}, x::W) where {C,V,W<:$(esc(T))} = one(C) * x
        Base.convert(::Type{NCAdd{C,NCMul{Int,W,V},D}}, x::W) where {C,V,D,W<:$(esc(T))} = NCAdd(zero(C), D(NCMul(1, [x]) => one(C)))

        VectorInterface.inner(x::MulAdd, y::$(esc(T))) = _inner(x, y)
        VectorInterface.inner(x::$(esc(T)), y::MulAdd) = _inner(x, y)
        VectorInterface.inner(x::$(esc(T)), y::$(esc(T))) = _inner(x, y)
        LinearAlgebra.norm(x::$(esc(T))) = sqrt(VectorInterface.inner(x, x))

        VectorInterface.add!!(x::MulAdd, y::$(esc(T)), α::Number, β::Number) = add!!(x, y + 0, α, β)
        VectorInterface.scale(x::$(esc(T)), α::Number) = α * x
        VectorInterface.scale!!(x::$(esc(T)), α::Number) = α * x
        VectorInterface.scale!!(a::NCAdd, x::$(esc(T)), α::Number) = α * x
        VectorInterface.zerovector(x::$(esc(T)), ::Type{S}) where {S<:Number} = VectorInterface.zerovector(NCMul(x), S)
    end
end
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
            push!(mul_pairs, :(Base.:*(x::$(esc(T1)), y::$(esc(T2))) = autosort() ? sort!(NCMul(1, [x, y])) : NCMul(1, [x, y])))
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
