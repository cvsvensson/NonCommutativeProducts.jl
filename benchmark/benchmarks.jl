using BenchmarkTools
using NonCommutativeProducts
using Random
const SUITE = BenchmarkGroup()
Random.seed!(1)

import NonCommutativeProducts: AddTerms, Swap, @nc
struct Fermion{L}
    label::L
    creation::Bool
end
Base.adjoint(x::Fermion) = Fermion(x.label, !x.creation)
Fermion(k) = Fermion(k, false)
Base.show(io::IO, x::Fermion) = print(io, "c", x.creation ? "†" : "", "[", x.label, "]")

@nc Fermion
NonCommutativeProducts.enable_autosort!()
function should_swap(a::Fermion, b::Fermion)
    if a.creation == b.creation
        return a.label > b.label
    else
        return a.creation < b.creation
    end
end

function NonCommutativeProducts.mul_effect(a::Fermion, b::Fermion)
    a == b && return 0 # a*b => 0
    if should_swap(a, b)
        a.label == b.label && xor(a.creation, b.creation) && return AddTerms((Swap(-1), 1)) #  a*b => -b*a + 1
        return Swap(-1) # a*b => -b*a
    else
        return nothing
    end
end

SUITE["symbolic_sum"] = @benchmarkable sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' for n in 1:1000)
SUITE["symbolic_sum_square"] = @benchmarkable sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' for n in 1:10)^3

labels = shuffle(1:10)
SUITE["symbolic_deep_product"] = @benchmarkable prod(Fermion(l) for l in labels) * prod(Fermion(l)' for l in labels)
