"""

@author: Denver-James Logan Marchment

This file contains functions related to the modular representation theory
of the group G = SL_2(F_p).

For the theoretical details, see: https://arxiv.org/abs/2503.07581

"""

import FGmodule
import FBmodule
from sys import exit

'''

Let c_{a,b,t} denote the multiplicity of V_t
as a composition factors of the Green correspondent
V_{a,b} of U_{a,b}.

Given a prime p, this file verifies that the coefficients
c_{a,b,t} as given in Corollary 3.18 of https://arxiv.org/abs/2503.07581
matches that given by the FGmodule and FBmodule framework I've
written

I've tested primes p in {5,7,11,13,17,19,23,29,31},
and indeed, both methods agree. It should be noted
that testing all combinations has computational
complexity O(p^3), so testing primes much beyond
becomes impactical without maybe some
multithreading

'''

# The two functions from the Corollary
def c_ab(p, a, b, t):
    lsum = 0
    rsum = 0

    # Case: a in [0,1]
    if 0 <= a <= 1:
        # b in [1, (p-1)/2]
        if 1 <= b <= (p - 1) / 2:
            # lsum condition on t
            if 1 <= t <= 2 * b + a - 1:
                lsum += 1
            # rsum condition on t
            if p - a - 2 * b + 1 <= t <= p - 1:
                rsum += 1

        # b in [(p+1)/2, p-1]
        elif (p + 1) / 2 <= b <= p - 1:
            # lsum condition on t
            if 1 <= t <= 2 * (p - b) - a:
                lsum += 1
            # rsum condition on t
            if a + 2 * b - p <= t <= p - 1:
                rsum += 1

    # Case: a in [2, (p-1)/2]
    elif 2 <= a <= (p - 1) / 2:
        # b in [1, (p-a)/2]
        if 1 <= b <= (p - a) / 2:
            if a <= t <= a - 1 + 2 * b:
                lsum += 1
            if p - a - 2 * b + 1 <= t <= p - 1:
                rsum += 1

        # b in [(p-a+1)/2, p-a]
        elif (p - a + 1) / 2 <= b <= p - a:
            if a <= t <= 2 * (p - b) - a:
                lsum += 1
            if a + 2 * b - p <= t <= p - 1:
                rsum += 1

        # b in [p-a, p-(a+1)/2]
        elif p - a <= b <= p - (a + 1) / 2:
            if 2 * (p - b) - a <= t <= a:
                lsum += 1
            if p - a <= t <= p - 1:
                rsum += 1

        # b in [p - a/2, p-1]
        elif p - a / 2 <= b <= p - 1:
            if 2 * (b - p) + a + 1 <= t <= a:
                lsum += 1
            if p - a <= t <= p - 1:
                rsum += 1

    # Case: a in [(p+1)/2, p-2]
    elif (p + 1) / 2 <= a <= p - 2:
        # b in [1, (p-a)/2]
        if 1 <= b <= (p - a) / 2:
            if p - a - 2 * b + 1 <= t <= p - a:
                lsum += 1
            if a <= t <= p - 1:
                rsum += 1

        # b in [(p-a+1)/2, p-a]
        elif (p - a + 1) / 2 <= b <= p - a:
            if a + 2 * b - p <= t <= p - a:
                lsum += 1
            if a <= t <= p - 1:
                rsum += 1

        # b in [p-a, p-(a+1)/2]
        elif p - a <= b <= p - (a + 1) / 2:
            if p - a <= t <= a + 2 * b - p:
                lsum += 1
            if 2 * (p - b) - a <= t <= p - 1:
                rsum += 1

        # b in [p - a/2, p-1]
        elif p - a / 2 <= b <= p - 1:
            if p - a <= t <= 2 * (p - b) - a - 1 + p:
                lsum += 1
            if 2 * (b - p) + a + 1 <= t <= p - 1:
                rsum += 1

    return lsum + rsum

def c_abt(p,a,b,t):
    if (t % 2 == a % 2):
        return 0
    if (1 <= t <= (p-1)/2):
        return c_ab(p,a,b,t)
    if ((p+1)/2 <= t <= p-1):
        return c_ab(p,a,b,p-t)

# We are comparing it to what happens when we manually
# calculate by inducing
def c_abt_manual(p,a,b,t):
    # First, if we are dealing with a projective, return
    # 0, since there projective modules do not have
    # Green correspondents
    if (b == p):
        return 0
    # Define the module
    N = FBmodule.FBmodule(p, [ [a,b] ])
    # Induce
    M = N.ind()
    # Get the non-projective summand of M
    M = FGmodule.FGmodule(p, [list(M.get_decomposition()[1].keys())[0]])
    # And now return the multiplicity of V_t
    # as a composition factor of M
    return M.get_composition_factors()[t-1]

# The prime being considered in this run
p = 11
# We check over all the relevant a,b,t
for t in range(1,p):
    for a in range(0,p-1):
        for b in range(1,p):
            if not (c_abt(p,a,b,t) == c_abt_manual(p,a,b,t)):
                print("When p=" + str(p) + ", a=" + str(a) + ", b=" + str(b) + ", t = " + str(t) + ", the Corollary and explicit computation don't match!")
                print("Manual: " + str(c_abt_manual(p,a,b,t)))
                print("Corollary: " + str(c_abt(p,a,b,t)))
                print("Quitting now.")
                exit()
print("When p=" + str(p) + ", the Corollary and the manual computation match! Success!")