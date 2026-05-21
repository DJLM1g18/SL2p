# SL(2,p)
A framework for the p-modular representation theory of the finite group G = SL(2,p), based on my [research paper](https://arxiv.org/abs/2503.07581). p is an odd prime, F is an algebraically closed
field of characteristic p > 2, and B is the Borel subgroup of upper triangular matrices of G.
## Files

### `FBmodule.py`.

The file `FBmodule.py` contains the class `FBmodule` which provides the user a mechanism to define F[B] modules. Defining a new F[B] module is as simple as:
```python
from FBmodule import FBmodule
# You can choose whatever odd prime you'd like to consider
p = 11
# The indecomposable modules for F[B] are of the form
# U_{a,b}, where 0 <= a <= p-2, and 1 <= b <= p.
# The parameter a is controlling the socle of the module,
# and the parameter b is controlling its dimension
a_1,b_1, a_2,b_2, a_3,b_3 = 0,1, 3,5, 2,11
modules = [ [a_1,b_1], [a_2,b_2], [a_3,b_3] ]
N = FBmodule(p, modules)
```
Methods include:

- `N.get_decomposition()`. Returns the triple `(p, non_projectives, projectives)`. `p` is the prime currently being considered, `non_projectives` is a dictionary containing the
non-projective indecomposable modules that occur as summands of `N` as keys, and the number of times they occur in the decomposition of `N` as values. `projectives` is a dictionary
containing the projective indecomposable modules that occur as summands of `N` as keys, and the number of times they occur in the decomposition of `N` as values.
- `N.get_dimension()`. Returns the dimension of `N` as a vector space over F.
- `N.get_composition_factors()`. Returns a tuple `[n_0, n_1, n_2, ..., n_{p-2}]` where n_a denotes the number of times the simple F[B] module S_a (see [Proposition 2.1](https://arxiv.org/abs/2503.07581))
occurs as a composition factor of `N`.
- `N.add(module_list)`. Allows the user to update `N` by adding modules to its direct-sum decomposition. The input `module_list` takes the format `module_list = [ [a_1,b_1], ..., [a_n,b_n] ]`.
- `N.remove(module_list)`. Allows the user to update `N` by removing modules from its direct sum decomposition. The input `module_list` takes the format `module_list = [ [a_1,b_1], ..., [a_n,b_n] ]`.
It removes the module `[a,b]` by however many times you include it in `module_list`. If you try to remove a module that is not currently a part of the direct sum decomposition of `N`, nothing
will happen.
- `N.ind()`. Returns an object of class `FGmodule`, which is the module `N` induced to the full group G = SL(2,p).

### `FGmodule.py`

The file `FGmodule.py` contains the class `FGmodule` which provides the user a mechanism to define F[G] modules. Defining a new F[G] module is as simple as:
```python
from FGmodule import FGmodule
# You can choose whatever odd prime you'd like to consider
p = 11
# There are two different parameterisation regimes for the indecomposable modules for G.
# The projective indecomposable modules are of the form P_{V_{t}}, where 1 <= t <= p.
# The non-projective indecomposable modules are of the form M(i,l,s,e), where:
# i is 0 or 1, corresponding to the two non-semisimple blocks of F[G].
# l is an integer satisfying 0 <= l <= (p-3)/2, corresponding to the index of the vertex
# that the walk on the Brauer tree starts on.
# s is an integer satisfying 1 <= s <= p-1-l, corresponding to the length of the walk
# on the Brauer tree.
# e is either +1 or -1, controlling whether the walk starts with a composition factor
# in the top or bottom.

i_1,l_1,s_1,e_1 = 0,0,4,-1
i_2,l_2,s_2,e_2 = 1,3,4,1
t_1,t_2,t_3 = 5,2,11
modules = [ [i_1,l_1,s_1,e_1], [i_2,l_2,s_2,e_2], [t_1], [t_2], [t_3] ]
M = FGmodule(p, modules)
```

Methods include:

- `M.get_decomposition()`. Returns the triple `(p, non_projectives, projectives)`. `p` is the prime currently being considered, `non_projectives` is a dictionary containing the
non-projective indecomposable modules that occur as summands of `M` as keys, and the number of times they occur in the decomposition of `M` as values. `projectives` is a dictionary
containing the projective indecomposable modules that occur as summands of `M` as keys, and the number of times they occur in the decomposition of `M` as values.
- `M.get_dimension()`. Returns the dimension of `N` as a vector space over F.
- `M.get_composition_factors()`. Returns a tuple `[n_1, n_2, n_3, ..., n_p]` where n_t denotes the number of times the simple F[G] module V_t (see [Proposition 2.3](https://arxiv.org/abs/2503.07581))
occurs as a composition factor of `M`.
- `M.add(module_list)`. Allows the user to update `M` by adding modules to its direct-sum decomposition. The input `module_list` takes the same format as in the variable `modules` defined in the
above sample code.
- `M.remove(module_list)`. Allows the user to update `M` by removing modules from its direct sum decomposition. The input `module_list` takes the same format as in the variable `modules` defined
in the above sample code. It removes the module by however many times you include it in `module_list`. If you try to remove a module that is not currently a part of the direct sum decomposition of `M`,
nothing will happen.
- `M.res()`. Returns an object of class `FBmodule`, which is the module `M` restricted to the Borel subgroup of upper triangular matrices B.

### `sample.py`

The file `sample.py` simply contains example usage of the above `FGmodule` and `FBmodule` classes.

### `symbolic-proof.py`

The file `symbolic-proof.py` was used to derive [Proposition 3.16](https://arxiv.org/abs/2503.07581). This result explicitly describes the Green correspondent of a given non-projective indecomposable
F[G] module. Describing this direction of the Green correspondence bijection involved many case distinctions (about 60), and so rather than check all 60 by hand - I found a way to enumerate the
cases, and then put the results into a table!

### `cabt.py`

Let c_{a,b,t} denote how many times the simple F[G] module V_t occurs as a composition factor of the Green correspondent of U_{a,b} (0 <= a <= p-2, 1 <= b <= p-1). Then,
[Corollary 3.18](https://arxiv.org/abs/2503.07581) gives an explicit formula for c_{a,b,t}. On the other hand, the above `FGmodule` and `FBmodule` frameworks allow us
to manually compute the c_{a,b,t} by:
1. Defining `U = FBmodule(p, [ [a,b] ])`.
2. Computing `M = U.ind()`.
3. Obtaining the unique non-projective indecomposable summand of `M` (for example, you can do `non_projective = list(M.get_decomposition()[1].keys())[0]`).
4. Then `cabt = (FGmodule(p [non_projective])).get_composition_factors()[t-1]`.

For a prime odd prime p, this file simply enumerates c_{a,b,t} using both Corollary 3.18 and manual computation, and verifies that indeed, both methods agree.
The code as written is O(p^3): I think I can make it O(p^2), and I may come back to do this. But for now, I have tested primes up to 31, and indeed,
both methods agree, as expected.
