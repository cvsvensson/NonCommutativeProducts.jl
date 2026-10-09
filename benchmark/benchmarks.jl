using BenchmarkTools
using NonCommutativeProducts
using Random
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1
const SUITE = BenchmarkGroup()
Random.seed!(1)

include("../test/rules/fermions.jl")
NonCommutativeProducts.enable_autosort!()
SUITE["symbolic_sum"] = @benchmarkable sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' for n in 1:1000)
SUITE["symbolic_sum_square"] = @benchmarkable sum(Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' for n in 1:10)^3

labels = shuffle(1:10)
SUITE["symbolic_deep_product"] = @benchmarkable prod(Fermion(l) for l in $labels) * prod(Fermion(l)' for l in $labels)

include("../test/rules/bosons.jl")
NonCommutativeProducts.@commutative Boson Fermion
# Boson is a single mode, and Boson() * Boson()' is the only boson pair that adds a term when reordered (b b† = b† b + 1).
# The Boson in the middle of the last term is moved to the left of the Fermion, which is a plain swap.
mixed_term(n) = Fermion(n)' * Fermion(n) + Fermion(n) * Fermion(n + 1)' + Boson()' * Fermion(n) * Boson()
SUITE["symbolic_sum_mixed"] = @benchmarkable sum(mixed_term(n) for n in 1:1000)
SUITE["symbolic_sum_square_mixed"] = @benchmarkable sum(mixed_term(n) for n in 1:8)^3

# A single product with lots of reordering but few additions: fermions with distinct labels only swap, and bosons
# only swap past fermions, so the extra terms come from the three b b† pairs.
mixed_ops = shuffle([Fermion.(1:6); adjoint.(Fermion.(7:12)); fill(Boson(), 3); fill(Boson()', 3)])
SUITE["symbolic_deep_product_mixed"] = @benchmarkable prod($mixed_ops)
