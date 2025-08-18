NCMul(f::NCAdd) = (length(f.dict) == 1 && iszero(f.coeff) && return prod(only(f.dict))) || throw(ArgumentError("Cannot convert NCAdd to NCMul: $f"))


function Base.:(==)(a::NCAdd, b::NCMul)
    iszero(a.coeff) || return false
    length(a.dict) == 1 || return false
    ncmul, coeff = only(a.dict)
    NCMul(coeff, ncmul.factors) == b
end
Base.:(==)(a::NCMul, b::NCAdd) = b == a

function Base.:+(_a::NCMul, _b::NCMul)
    a, b = promote(_a, _b)
    if a.factors == b.factors
        return NCAdd(0, Dict(NCMul(1, a.factors) => a.coeff + b.coeff))
    end
    return NCAdd(0, Dict(NCMul(1, a.factors) => a.coeff, NCMul(1, b.factors) => b.coeff))
end

Base.:+(a::Number, b::NCMul) = NCAdd(a, to_add(b))
Base.:+(a::UniformScaling, b::NCMul) = NCAdd(a.λ, to_add(b))
Base.:+(a::NCMul, b::Union{Number,UniformScaling}) = b + a
Base.:+(a::NCMul, b::NCAdd) = NCAdd(b.coeff, mergewith!!(+, to_add(a), b.dict))

# function Base.:+(a::NCMul{CA,KA,FA}, b::NCMul{CB,KB,FB}) where {CA,CB,KA,KB,FA,FB}
#     C = promote_type(CA, CB)
#     K = Union{NCMul{C,KA,FA},NCMul{C,KB,FB}}
#     D = Dict{K,C}
#     if a.factors == b.factors
#         coeff = a.coeff + b.coeff
#         return (NCAdd(0, D(NCMul(1, a.factors) => coeff)))
#     end
#     at, bt = to_add_tuple(a), to_add_tuple(b)
#     return (NCAdd(0, D(at..., bt...)))
# end

to_add(a::NCMul, coeff=1) = Dict(NCMul(1, a.factors) => a.coeff * coeff)
to_add(a, coeff=1) = Dict(NCMul(a) => coeff)
to_add_tuple(a::NCMul, coeff=1) = (NCMul(1, a.factors) => a.coeff * coeff,)
to_add_tuple(a, coeff=1) = (NCMul(a) => coeff,)

Base.:^(a::Union{NCMul,NCAdd}, b) = Base.power_by_squaring(a, b)

macro nc(T)
    quote
        Base.:+(x::$(esc(T)), y::$(esc(T))) = NCMul(1, [x]) + NCMul(1, [y])
        Base.:+(x::$(esc(T)), y::Union{Number,UniformScaling,NCMul,NCAdd}) = NCMul(1, [x]) + y
        Base.:+(x::Union{Number,UniformScaling,NCMul,NCAdd}, y::$(esc(T))) = y + x

        NonCommutativeProducts.add!!(x::NCAdd, y::$(esc(T))) = add!!(x, NCMul(1, [y]))
        NonCommutativeProducts.eager(::T) where T<:$(esc(T)) = eager(T)
        NonCommutativeProducts.eager(::Type{<:$(esc(T))}) = false

        Base.:-(x::$(esc(T)), y::$(esc(T))) = NCMul(1, [x]) - NCMul(1, [y])
        Base.:-(x::$(esc(T))) = NCMul(-1, [x])
        Base.:-(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = x - NCMul(1, [y])
        Base.:-(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = NCMul(1, [x]) - y

        Base.:*(x::$(esc(T)), y::$(esc(T))) = NCMul(1, [x, y])
        Base.:*(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = x * NCMul(1, [y])
        Base.:*(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = NCMul(1, [x]) * y
        Base.:*(x::$(esc(T)), y::NCMul) = NCMul(y.coeff, pushfirst!!(copy(y.factors), x))
        Base.:*(x::NCMul, y::$(esc(T))) = NCMul(x.coeff, push!!(copy(x.factors), y))
        Base.:^(a::$(esc(T)), b) = Base.power_by_squaring(a, b)
        Base.convert(::Type{NCMul{C,S,F}}, x::$(esc(T))) where {C,S<:$(esc(T)),F} = NCMul(one(C), S[x])
    end
end
eager(::NCMul{C,T}) where {C,T} = eager(T)
Ordering(::NCMul{C,T}) where {C,T} = Ordering(T)
# const ADD{T} = NCAdd{C,MUL{T},D} where {C,D}
macro nc_eager(T, ordering)
    quote
        # const _ordering = $(esc(ordering))
        # MUL = NCMul{C,<:$(esc(T)),F} where {C,F}
        # ADD = NCAdd{T1,NCMul{T2,<:$(esc(T)),T3},T4} where {T1,T2,T3,T4}
        Base.:+(x::$(esc(T)), y::$(esc(T))) = NCMul(1, [x]) + NCMul(1, [y])
        Base.:+(x::$(esc(T)), y::Union{Number,UniformScaling,NCMul,NCAdd}) = NCMul(1, [x]) + y
        Base.:+(x::Union{Number,UniformScaling,NCMul,NCAdd}, y::$(esc(T))) = y + x

        NonCommutativeProducts.add!!(x::NCAdd, y::$(esc(T))) = add!!(x, NCMul(1, [y]))

        # NonCommutativeProducts.Ordering(::<:$(esc(T))) = $(esc(ordering))
        NonCommutativeProducts.Ordering(::T) where T<:$(esc(T)) = Ordering(T)
        NonCommutativeProducts.Ordering(::Type{<:$(esc(T))}) = $(esc(ordering))
        NonCommutativeProducts.eager(::T) where T<:$(esc(T)) = eager(T)
        NonCommutativeProducts.eager(::Type{<:$(esc(T))}) = true
        # NonCommutativeProducts.eager(::NCMul{C,$(esc(T))}) where C = true
        # function Base.:*(a::NCAdd{T1,<:NCMul{T2,<:$(esc(T)),T3},T4}, b::NCMul{C,<:$(esc(T)),F}) where {T1,T2,T3,T4,C,F}
        #     c = zero(a)
        #     return bubble_sort(mul!!(c, a, b), _ordering)
        # end
        # function Base.:*(a::NCMul{C,<:$(esc(T)),F}, b::NCAdd{T1,<:NCMul{T2,<:$(esc(T)),T3},T4}) where {C,F,T1,T2,T3,T4}
        #     c = zero(b)
        #     return bubble_sort(mul!!(c, a, b), $(esc(ordering)))
        # end

        # function Base.:*(a::$(esc(ADD)), b::$(esc(ADD)))
        # function Base.:*(a::NCAdd{T1,<:NCMul{T2,<:$(esc(T)),T3},T4}, b::NCAdd{T12,<:NCMul{T22,<:$(esc(T)),T32},T42}) where {T1,T2,T3,T4,T12,T22,T32,T42}
        #     c = zero(a)
        #     bubble_sort(mul!!(c, a, b), $(esc(ordering)))
        # end

        Base.:-(x::$(esc(T)), y::$(esc(T))) = NCMul(1, [x]) - NCMul(1, [y])
        Base.:-(x::$(esc(T))) = NCMul(-1, [x])
        Base.:-(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = x - NCMul(1, [y])
        Base.:-(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = NCMul(1, [x]) - y

        Base.:*(x::$(esc(T)), y::$(esc(T))) = bubble_sort(NCMul(1, [x, y]), $(esc(ordering)))
        Base.:*(x::Union{Number,UniformScaling,NCAdd}, y::$(esc(T))) = bubble_sort(x * NCMul(1, [y]), $(esc(ordering)))
        Base.:*(x::$(esc(T)), y::Union{Number,UniformScaling,NCAdd}) = bubble_sort(NCMul(1, [x]) * y, $(esc(ordering)))
        Base.:*(x::$(esc(T)), y::NCMul) = bubble_sort(NCMul(y.coeff, pushfirst!!(copy(y.factors), x)), $(esc(ordering)))
        Base.:*(x::NCMul, y::$(esc(T))) = bubble_sort(NCMul(x.coeff, push!!(copy(x.factors), y)), $(esc(ordering)))
        Base.:^(a::$(esc(T)), b) = bubble_sort(Base.power_by_squaring(a, b), $(esc(ordering)))
        Base.convert(::Type{NCMul{C,S,F}}, x::$(esc(T))) where {C,S<:$(esc(T)),F} = NCMul(one(C), S[x])
    end
end

