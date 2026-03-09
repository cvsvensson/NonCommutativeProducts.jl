NCMul(f::NCAdd) = (length(f.dict) == 1 && iszero(f.coeff) && return prod(only(f.dict))) || throw(ArgumentError("Cannot convert NCAdd to NCMul: $f"))


function Base.:(==)(a::NCAdd, b::NCMul)
    iszero(a.coeff) || return false
    length(a.dict) == 1 || return false
    ncmul, coeff = only(a.dict)
    NCMul(coeff, ncmul.factors) == b
end
Base.:(==)(a::NCMul, b::NCAdd) = b == a

function Base.:+(a::NCMul{C1,F1,S}, b::NCMul{C2,F2,S}) where {C1,C2,F1,F2,S}
    if a.factors == b.factors
        return NCAdd(0, Dict(NCMul(1, a.factors) => a.coeff + b.coeff))
    end
    F = Union{F1,F2}
    C = promote_type(C1, C2)
    return NCAdd(0, Dict{NCMul{C,F,S},C}(NCMul{C,F,S}(1, a.factors) => a.coeff, NCMul{C,F,S}(1, b.factors) => b.coeff))
end
function Base.:+(a::NCMul{C1,F1,S1}, b::NCMul{C2,F2,S2}) where {C1,C2,F1,F2,S1,S2}
    if a.factors == b.factors
        return NCAdd(0, Dict(NCMul(1, a.factors) => a.coeff + b.coeff))
    end
    F = promote_type(F1, F2)
    C = promote_type(C1, C2)
    S = promote_type(S1, S2)
    return NCAdd(0, Dict{NCMul{C,F,S},C}(NCMul{C,F,S}(1, a.factors) => a.coeff, NCMul{C,F,S}(1, b.factors) => b.coeff))
end

Base.:+(a::Number, b::NCMul) = NCAdd(a, to_add_dict(b))
Base.:+(a::UniformScaling, b::NCMul) = NCAdd(a.λ, to_add_dict(b))
Base.:+(a::NCMul, b::Union{Number,UniformScaling}) = b + a
# Base.:+(a::NCMul, b::NCAdd) = NCAdd(b.coeff, mergewith!!(+, to_add_dict(a), b.dict))
function Base.:+(a::NCMul, b::NCAdd)
    newdict = copy(b.dict)
    for (k, v) in b.dict
        if k.factors == a.factors
            newdict[k] = v + a.coeff
            return NCAdd(b.coeff, newdict)
            # return NCAdd(b.coeff, setindex!!(newdict, k => v + ))
        end
    end
    # println(newdict)
    # println(a.factors)
    # println(a.coeff)
    # println(setindex!!(newdict, NCMul(1, a.factors), a.coeff))
    nc = NCAdd(b.coeff, setindex!!(newdict, a.coeff, NCMul(1, a.factors)))
    return nc

    # NCAdd(b.coeff, mergewith!!(+, to_add_dict(a), b.dict))
end
Base.convert(::Type{NCAdd{C,NCMul{C,S,F},D}}, x::NCMul{C,S,F}) where {C,S,F,D} = NCAdd(zero(C), D(to_add_dict(x)))

to_add_dict(a::NCMul) = Dict(NCMul(1, a.factors) => a.coeff)

function Base.:^(a::Union{NCAdd,NCMul}, b::Int)
    ret = Base.power_by_squaring(a, b)
    autosort() && return bubble_sort(ret)
    return ret
end

macro nc_common(T)
    quote
        NonCommutativeProducts.NCMul(f::$(esc(T))) = NCMul(1, [f])

        Base.:+(x::$(esc(T)), y::$(esc(T))) = NCMul(x) + NCMul(y)
        Base.:+(x::$(esc(T)), y::Union{Number,UniformScaling,NCMul,NCAdd}) = NCMul(x) + y
        Base.:+(x::Union{Number,UniformScaling,NCMul,NCAdd}, y::$(esc(T))) = y + x

        NonCommutativeProducts.add!!(x::NCAdd, y::$(esc(T))) = add!!(x, NCMul(y))

        Base.:-(x::$(esc(T)), y::$(esc(T))) = NCMul(x) - NCMul(y)
        Base.:-(x::$(esc(T))) = NCMul(-1, [x])
        Base.:-(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = x - NCMul(y)
        Base.:-(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = NCMul(x) - y

        Base.:*(x::$(esc(T)), y::$(esc(T))) = autosort() ? bubble_sort!(NCMul(1, [x, y])) : NCMul(1, [x, y])
        Base.:*(x::$(esc(T)), y::NCMul) = autosort() ? bubble_sort!(NCMul(y.coeff, pushfirst!!(copy(y.factors), x))) : NCMul(y.coeff, pushfirst!!(copy(y.factors), x))
        Base.:*(x::NCMul, y::$(esc(T))) = autosort() ? bubble_sort!(NCMul(x.coeff, push!!(copy(x.factors), y))) : NCMul(x.coeff, push!!(copy(x.factors), y))
        Base.:*(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = x * NCMul(y)
        Base.:*(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = NCMul(x) * y

        Base.:^(a::$(esc(T)), b) = NCMul(a)^b
        Base.convert(::Type{NCMul{C,S,F}}, x::$(esc(T))) where {C,S<:$(esc(T)),F} = NCMul(one(C), S[x])

        Base.:(==)(a::$(esc(T)), b::Union{NCMul,NCAdd}) = NCMul(a) == b
        Base.:(==)(a::Union{NCMul,NCAdd}, b::$(esc(T))) = b == a
    end
end


const _autosort = Ref(false)

autosort() = _autosort[]
enable_autosort!() = _autosort[] = true
disable_autosort!() = _autosort[] = false


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
            push!(mul_effect, :(NonCommutativeProducts.mul_effect(x::$(esc(T2)), y::$(esc(T1))) = Swap(1)))
        end
    end
    quote
        $(mul_effect...)
        @nc_pairs $(esc.(types)...)
    end
end
