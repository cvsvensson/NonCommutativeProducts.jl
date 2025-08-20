# NonCommutativeProducts

<!--  [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cvsvensson.github.io/NonCommutativeProducts.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cvsvensson.github.io/NonCommutativeProducts.jl/dev/) -->
[![Build Status](https://github.com/cvsvensson/NonCommutativeProducts.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cvsvensson/NonCommutativeProducts.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cvsvensson/NonCommutativeProducts.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cvsvensson/NonCommutativeProducts.jl)

**NonCommutativeProducts.jl** is a Julia package for sorting non-commuting objects, such as operators in quantum mechanics. Users must specify custom commutation relations and sorting orders.

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
Then, we need to hook up our struct to the package to let it handle the arithmetic. Use either `@nc Fermion` or `@nc_eager Fermion Ordering`. The difference is that `@nc` only sorts things when explicitly prompted, while `@nc_eager` will apply them eagerly according to `Ordering`. Let's do it eagerly here.
```julia
using NonCommutativeProducts
import NonCommutativeProducts: Swap, AddTerms, @nc_eager, mul_effect
struct NormalOrder end
@nc_eager Fermion NormalOrder()
```
Now let's define what happens when we multiply two fermions under normal ordering.
```julia
function mul_effect(a::Fermion, b::Fermion, ::NormalOrder)
    # If the fermion is multiplied with itself, we replace it by zero. 
    a.label == b.label && a.dagger == b.dagger && return 0 
    # if a is annihilation and b is creation, we should swap them
    if !a.dagger && b.dagger 
        if a.label !== b.label 
            return Swap(-1) # Swaps them and multiplies by -1
        end
        # If their labels are the same, swap and add an extra term
        return AddTerms((Swap(-1), 1)) 
    else
        return nothing # No effect
    end
end
```
Now we can multiply and add these fermions. Normal order will be applied directly.
```julia
Fermion(1)'*Fermion(1)
#c†[1]*c[1]
Fermion(1)*Fermion(1)'
#1I - c†[1]*c[1]
prod(Fermion(n) + Fermion(n)' for n in 1:4)
#=Sum with 16 terms: 
- c†[2]*c†[4]*c[1]*c[3] + c†[3]*c†[4]*c[1]*c[2] - c†[2]*c†[3]*c†[4]*c[1] + ...=#
```
Note that with this ordering, some terms might be equivalent to others under further swaps. One could also sort the terms via their label to get a unique representation, which is done in the examples in the tests of this package. 

## Remarks

This package is flexible, but not very efficient. Sorting is done via bubble sort, which is suitable for this use case since it is based on swapping adjacent elements. But it does not scale well with the length of the list, so it won't perform well for products of many elements.

```julia
@time op = prod(Fermion(n) + Fermion(n)' + 1 for n in 1:10)
#=0.327861 seconds (3.05 M allocations: 129.336 MiB, 12.88% compilation time)
Sum with 59048 terms: 
c†[8]*c[2]*c[3]*c[5]*c[6]*c[9]*c[10] + c†[1]*c†[4]*c†[6]*c†[8]*c†[10]*c[2]*c[5]*c[9] - c†[1]*c†[2]*c†[6]*c[5]*c[7]*c[10] + ...=#
```


To cut down on allocations, one can use `NonCommutativeProducts.add!!` which tries to perform addition in place, but widens the type if not possible.
```julia
op = 0
@time for n in 1:100
    op += Fermion(n)'*Fermion(n)
end
 #  0.000192 seconds (1.84 k allocations: 371.699 KiB)

op2 = zero(op)
@time for n in 1:100
    op2 = NonCommutativeProducts.add!!(op2, Fermion(n)'*Fermion(n))
end
#  0.000130 seconds (1.41 k allocations: 83.998 KiB)
op == op2 #true
```
