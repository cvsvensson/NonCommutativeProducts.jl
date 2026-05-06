using BenchmarkTools
using NonCommutativeProducts
using Random
const SUITE = BenchmarkGroup()
Random.seed!(1)

include("../test/rules/fermions.jl")

SUITE["symbolic_sum"] = @benchmarkable sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' for n in 1:1000)
SUITE["symbolic_sum_square"] = @benchmarkable sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' for n in 1:10)^3

labels = shuffle(1:10)
SUITE["symbolic_deep_product"] = @benchmarkable prod(Fermion(l) for l in labels) * prod(Fermion(l)' for l in labels)
