# NonCommutativeProducts.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cvsvensson.github.io/NonCommutativeProducts.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cvsvensson.github.io/NonCommutativeProducts.jl/dev/)
[![Build Status](https://github.com/cvsvensson/NonCommutativeProducts.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cvsvensson/NonCommutativeProducts.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cvsvensson/NonCommutativeProducts.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cvsvensson/NonCommutativeProducts.jl)
[![Aqua QA](https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Version](https://img.shields.io/github/v/tag/cvsvensson/NonCommutativeProducts.jl)](https://github.com/cvsvensson/NonCommutativeProducts.jl/tags)

**NonCommutativeProducts.jl** is a Julia package for sorting non-commuting objects, such as operators in quantum mechanics. Users must specify custom commutation relations and sorting orders, as there are no inbuilt ones in this package. With a list of types `Ts` that represents your non-commuting objects, call `@nc Ts...` to define addition and multiplication for these types, and define `mul_effect(a, b)` for each pair of types to say what `a*b` is replaced by while sorting.

## Installation

```julia
using Pkg
Pkg.add("NonCommutativeProducts")
```

## Example: Fermions

Fermions sorted in normal order, with all creation operators before annihilation operators:
```julia
struct Fermion
    label::Int
    dagger::Bool
end
Base.adjoint(x::Fermion) = Fermion(x.label, !x.dagger)
Fermion(k) = Fermion(k, false)
Base.show(io::IO, x::Fermion) = print(io, "c", x.dagger ? "†" : "", "[", x.label, "]")

using NonCommutativeProducts
import NonCommutativeProducts: @nc, mul_effect
@nc Fermion

function mul_effect(a::Fermion, b::Fermion)
    # If the fermion is multiplied with itself, replace it by zero.
    (a.dagger, a.label) == (b.dagger, b.label) && return 0
    # If already sorted, return nothing
    (!a.dagger, a.label) < (!b.dagger, b.label) && return nothing
    # Non-trivial anti-commutation relation
    (!a.dagger, a.label) == (b.dagger, b.label) && return -b*a + 1
    # Trivial anti-commutation relation
    return -b*a
end

Fermion(1)'*Fermion(1) |> sort
#c†[1]*c[1]
Fermion(1)*Fermion(1)' |> sort
#1I - c†[1]*c[1]
```

## Documentation

See the [stable](https://cvsvensson.github.io/NonCommutativeProducts.jl/stable/) or [development](https://cvsvensson.github.io/NonCommutativeProducts.jl/dev/) documentation for automatic sorting, multiple types, traversing expressions, linear algebra with [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl), gradings, performance tips and the API reference.
