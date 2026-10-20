# API reference

```@meta
CurrentModule = NonCommutativeProducts
```

No names are exported, so access them as `NonCommutativeProducts.name` or import them explicitly, e.g. `using NonCommutativeProducts: @nc, terms`.

```@docs
NonCommutativeProducts
```

## Setting up types

```@docs
NonCommutativeProducts.@nc
NonCommutativeProducts.@commutative
NonCommutativeProducts.mul_effect
NonCommutativeProducts.Swap
NonCommutativeProducts.AddTerms
```

## Gradings

```@docs
NonCommutativeProducts.grade
NonCommutativeProducts.mul_grade
NonCommutativeProducts.adjoint_grade
NonCommutativeProducts.grading
NonCommutativeProducts.ProductGrade
NonCommutativeProducts.grades
NonCommutativeProducts.grade_part
NonCommutativeProducts.check_grading
```

## Sorting

```@docs
Base.sort(::NonCommutativeProducts.NCMul)
NonCommutativeProducts.enable_autosort!
NonCommutativeProducts.disable_autosort!
NonCommutativeProducts.autosort
NonCommutativeProducts.with_autosort
```

## Expressions

```@docs
NonCommutativeProducts.NCMul
NonCommutativeProducts.NCAdd
NonCommutativeProducts.terms
NonCommutativeProducts.factors
NonCommutativeProducts.coefficient
NonCommutativeProducts.scalar
NonCommutativeProducts.add!!
NonCommutativeProducts.ncmap
NonCommutativeProducts.ncmapreduce
```

## Index

```@index
Pages = ["api.md"]
```
