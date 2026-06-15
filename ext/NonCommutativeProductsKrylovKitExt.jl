module NonCommutativeProductsKrylovKitExt

using NonCommutativeProducts
import KrylovKit

const NC = NonCommutativeProducts

# Treat NC expressions as linear operators acting by left multiplication.
KrylovKit.apply(operator::NC.MulAdd, x) = operator * x
KrylovKit.apply_normal(operator::NC.MulAdd, x) = operator * x
KrylovKit.apply_adjoint(operator::NC.MulAdd, x) = operator' * x

end
