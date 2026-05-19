using BenchmarkTools
using NonCommutativeProducts
using Random
const SUITE = BenchmarkGroup()
Random.seed!(1)

include("../test/rules/fermions.jl")
NonCommutativeProducts.enable_autosort!()
SUITE["symbolic_sum"] = @benchmarkable sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' for n in 1:1000)
SUITE["symbolic_sum_square"] = @benchmarkable sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' for n in 1:10)^3

labels = shuffle(1:10)
SUITE["symbolic_deep_product"] = @benchmarkable prod(Fermion(l) for l in labels) * prod(Fermion(l)' for l in labels)

include("../test/rules/bosons.jl")
NonCommutativeProducts.@commutative Boson Fermion
SUITE["symbolic_sum_mixed"] = @benchmarkable sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' + Boson(n)' * Boson(n + 1) for n in 1:1000)
SUITE["symbolic_sum_square_mixed"] = @benchmarkable sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' + Boson(n)' * Boson(n + 1) for n in 1:10)^3

labels = shuffle(1:6)
SUITE["symbolic_deep_product_mixed"] = @benchmarkable prod(Fermion(l) for l in labels)* prod(Boson(l) for l in labels) * prod(Fermion(l)' for l in labels)