# NonCommutativeProducts.jl

<!--  [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cvsvensson.github.io/NonCommutativeProducts.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cvsvensson.github.io/NonCommutativeProducts.jl/dev/) -->
[![Build Status](https://github.com/cvsvensson/NonCommutativeProducts.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cvsvensson/NonCommutativeProducts.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cvsvensson/NonCommutativeProducts.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cvsvensson/NonCommutativeProducts.jl)

**NonCommutativeProducts.jl** is a Julia package for sorting non-commuting objects, such as operators in quantum mechanics. Users must specify custom commutation relations and sorting orders, as there are no inbuilt ones in this package. 

## How to do it
With a list of types `Ts` that represents your non-commuting objects, call `@nc Ts...` which defines addition and multiplication for these types. In order to sort them, define `mul_effect(a::Tj, b::Tk)` for each pair of types to define the behaviour of `a*b`.

The function `mul_effect(a::Tj, b::Tk)` can have return values 
* `nothing`: Keeps `a*b`. This return value is important as the sorting only terminates when this is the return value for each neighbouring pair product. If you don't have this, you'll get stuck in an infinite loop.
* `λ::Number`: replaces `a*b` by `λ`.
* `x::T`: replaces `a*b` by `x`
* `NCMul(λ::Number, xs::Vector{T})`: replaces `a*b` by `λ*x[1]*x[2]...`
* `Swap(λ::Number)`: Replaces `a*b` by `λ*b*a`. This is an instance of the rule above but more convenient.
* `AddTerms(terms)`: `a*b` should be replaced by a sum of terms. `terms` should be an iterable such as a vector or a tuple, and the elements can be of the four types above.

This is a young package and not every combination has been tested thoroughly. Do write tests for your specific use case to verify that it works and please report any bugs here.

No names are currently exported from the package, so you'll need to import the names you need.

Let's look at an example that shows how to define and sort fermions. The tests of this package contains a similar example with majorana fermions, and another test shows how one can collect powers.

## Example: Fermions

Let's see how to define fermions which satisfy
```math
\begin{align*}
\{c_i^\dagger,c_j\}&=\delta_{ij}\\ 
\{c_i,c_j\}&=0\\ 
\{c_i^\dagger,c_j^\dagger\}&=0.
\end{align*}
```

We'll sort them in normal order i.e. all creation operators appear before annihilation operators. First let's define a struct representing a fermion.
```julia
struct Fermion
	label::Int
    dagger::Bool
end
Base.adjoint(x::Fermion) = Fermion(x.label, !x.dagger)
Fermion(k) = Fermion(k, false)
Base.show(io::IO, x::Fermion) = print(io, "c", x.dagger ? "†" : "", "[", x.label, "]")
```
Then, we need to hook up our struct to the package to let it handle the arithmetic. Let's load the package and import the functions we are gonna use.
```julia
using NonCommutativeProducts
import NonCommutativeProducts: Swap, AddTerms, @nc, mul_effect
@nc Fermion
```
Now one can create expressions with these fermions,
```julia
Fermion(1)'*Fermion(1) + 1
```
but they can't be sorted. To sort them in normal order, we define 
```julia
function mul_effect(a::Fermion, b::Fermion)
    # If the fermion is multiplied with itself, we replace it by zero. 
    a.label == b.label && a.dagger == b.dagger && return 0 
    # if a is annihilation and b is creation, we should swap them
    if !a.dagger && b.dagger 
        if a.label != b.label 
            return Swap(-1) # Swaps them and multiplies by -1
        end
        # If their labels are the same, swap and add an extra term
        return AddTerms((Swap(-1), 1)) 
    elseif a.dagger == b.dagger && a.label > b.label
        return Swap(-1) # If b has smaller label than a, swap them and multiply by -1
    else
        return nothing # No effect
    end
end
```
Now we can sort expression involving fermions while respecting the commutation relations.
```julia
Fermion(1)'*Fermion(1) |> sort
#c†[1]*c[1]
Fermion(1)*Fermion(1)' |> sort
#1I - c†[1]*c[1]
```
In order automatically sort them on each multiplication, we can call `enable_autosort!`:
```julia
NonCommutativeProducts.enable_autosort!()
Fermion(1)*Fermion(1)'
prod(Fermion(n) + Fermion(n)' for n in 1:4) 
#=Sum with 16 terms: 
 -c†[1]*c†[2]*c†[4]*c[3] + c†[1]*c[2]*c[3]*c[4] + c†[1]*c†[3]*c†[4]*c[2] + ...=#
```

## Remarks

This package is flexible, but not very efficient. Sorting is done via bubble sort, which is suitable for this use case since it is based on swapping adjacent elements. But it does not scale well with the length of the list, so it won't perform well for products of many elements.

```julia
@time op = prod(Fermion(n) + Fermion(n)' + 1 for n in 1:10)
#=0.264353 seconds (4.34 M allocations: 155.466 MiB, 24.44% gc time, 8.75% compilation time)
Sum with 59049 terms: 
1I-c†[1]*c†[2]*c†[6]*c[5]*c[7]*c[10]+c†[8]*c[2]*c[3]*c[5]*c[6]*c[9]*c[10]+c†[1]*c†[4]*c†[6]*c†[8]*c†[10]*c[2]*c[5]*c[9] + ...=#
```


To cut down on allocations, one can use `NonCommutativeProducts.add!!` which tries to perform addition in place, but widens the type if not possible.
```julia
op = 0
@time for n in 1:100
    op += Fermion(n)'*Fermion(n)
end
 #  0.000251 seconds (2.55 k allocations: 311.127 KiB)

op2 = zero(op)
@time for n in 1:100
    op2 = NonCommutativeProducts.add!!(op2, Fermion(n)'*Fermion(n))
end
#  0.000192 seconds (2.12 k allocations: 106.404 KiB)
op == op2 #true
```
