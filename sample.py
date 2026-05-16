"""

@author: Denver-James Logan Marchment

This file contains functions related to the modular representation theory
of the group G = SL_2(F_p).

For the theoretical details, see: https://arxiv.org/abs/2503.07581

"""

import FBmodule
import FGmodule

print()
# Example usage defining and restricting an F[G] module
# For p=7, we restrict the module:
# M(0,0,4,-1) + P_{V_3}
p = 7
M = FGmodule.FGmodule(p,[[0,0,4,-1], [3]])
# We can print the composition factors of M
print("Composition factors of M: " + str(M.get_composition_factors()))
# ... and it's dimension
print("Dimension of M: " + str(M.get_dimension()))
# We restrict M
N = M.res()
# ... and print the resulting decomposition
print("Decomposition of N := Res^G_B(M): " + str(N.get_decomposition()))
print()

# Example usage defining and inducting an F[B] module
# For p=11, we induct the module:
# U_{0,3} + U_{6,4} + U_{9,11}
p = 11
W = FBmodule.FBmodule(p, [[0,3], [6,4], [9,11]])
# We can print the composition factors of W
print("Composition factors of W: " + str(W.get_composition_factors()))
# ... and it's dimension
print("Dimension of W: " + str(W.get_dimension()))
# We induct W
X = W.ind()
# ... and print the resulting decomposition
print("Decomposition of X := Res^G_B(W): " + str(X.get_decomposition()))
# Notice that indeed, the dimension of X is
# (p+1) * dim(W), which is what you expect
# since: |G/B| = p+1
dimW = W.get_dimension()
dimX = X.get_dimension()
if (dimW * (p+1) == dimX):
    print("We have dim(W) * (p+1) = " + str(dimW*(p+1)) + ", and dim(X) = " + str(dimX) + ", which is what you'd expect.")
print()