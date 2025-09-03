TermInterface.isexpr(::NCAdd) = true
TermInterface.iscall(::NCAdd) = true
TermInterface.head(::NCAdd) = :call
TermInterface.operation(::NCAdd) = +
TermInterface.children(x::NCAdd) = arguments(x)

TermInterface.isexpr(::NCMul) = true
TermInterface.iscall(::NCMul) = true
TermInterface.head(::NCMul) = :call
TermInterface.operation(::NCMul) = *
TermInterface.children(x::NCMul) = arguments(x)

function TermInterface.arguments(x::NCAdd{C1,NCMul{C2,T,F}}) where {C1,C2,T,F}
    C = promote_type(C1, C2)
    args = NCMul{C,T,F}[NCMul(C(x.coeff), T[])]
    for (k, v) in x.dict
        args = push!(args, v * k)
    end
    return args
end
function TermInterface.maketerm(::Type{NCAdd}, head, args, metadata)
    @assert head == :+ "Head must be :+ for NCAdd"
    sum(args)
end

function TermInterface.arguments(x::NCMul{C,T,F}) where {C,T,F}
    args = Union{C,T}[x.coeff]
    append!(args, x.factors)
    args
end
function TermInterface.maketerm(::Type{NCMul}, head, args, metadata)
    @assert head == :* "Head must be :* for NCMul"
    prod(args)
end
