# -*- coding: utf-8 -*-
"""

@author: Denver-James Logan Marchment

This file contains functions related to the modular representation theory
of the group G = SL_2(F_p).

"""

'''
This script symbolically takes an F[G] module M(i,l,s,e)
and gives a formula for the corresponding a value as in U_{a,b}.
'''

from sympy import *
import FGmodule
import FBmodule

p,i,l,s,e,d,a = symbols("p i l s e d a")

'''

iflag: 0 means i=0, 1 means i=1

lflag: 0 means l is zero, 1 means l is odd,
2 means l is even >= 2.

sflag: 0 means s is even >= 2, 1 means s odd.

eflag: -1 means e=-1, +1 means e=1

sumflag: 0 means l+s <= (p-1)/2,
1 means (p-1)/2 < l + s < (p-1),
2 meams l + s = p-1

'''

# This function checks if two tuples agree everywhere except at most
# one spot, and if so returns the spot (returns -1 otherwise)
def onedifference(tuple1,tuple2):
    n = len(tuple1)
    if not (n == len(tuple2)):
        return -1
    differences = 0
    last_difference_spot = -1
    for i in range(0,n):
        if not (tuple1[i] == tuple2[i]):
            differences += 1
            last_difference_spot = i
    if (differences == 1):
        return last_difference_spot 
    return -1

# This modified substitution function ensures that
# an error doesn't occur if the expression we have
# is actually just an integer
def subs(expression,var1,var2):
    if not isinstance(expression,int):
        expression = expression.subs(var1,var2)
    return expression

# This function allows us to find out if the expression
# in l,s is odd or even
def collapse_expression_F2(expression, lflag, sflag):
    if isinstance(expression,int):
        return expression % 2
    expression = expression.subs(p,1)
    expression = expression.subs(l,lflag)
    expression = expression.subs(s,sflag)
    return expression % 2

# Get the left hook of the module
def lefthook(iflag,lflag, eflag):
    if (lflag == 0 and eflag == 1):
        return (iflag,0,1,-1)
    if not(lflag == 0) and (eflag == 1):
        return (iflag,l-1,2,1)
    return (iflag,l,2,-1)

# Get the right hook of the module
def righthook(iflag,lflag, sflag, eflag, sumflag):
    # We first deal with the case that we don't loop around the exceptional vertex
    if (sumflag == 0 and eflag == 1 and sflag == 0):
        # dflag == -1
        return (iflag,l+s-2,2,1)
    if (sumflag == 0 and eflag == 1 and sflag == 1):
        # dflag == 1
        return (iflag,l+s-1,2,-1)
    if (sumflag == 0 and eflag == -1 and sflag == 0):
        # dflag == 1
        return (iflag,l+s-1,2,-1)
    if (sumflag == 0 and eflag == -1 and sflag == 1):
        # dflag == -1
        return (iflag,l+s-2,2,1)
    # If we do loop around the exceptional vertex, we first deal
    # with the case that we don't come all the way back.
    if (sumflag == 1):
        lnew = 2
        if (collapse_expression_F2(p-1-l-s,lflag,sflag) == 1):
            lnew = 1
        enew = eflag
        if (sflag == 0):
            enew = -eflag
        righthook = lefthook(iflag,lnew,enew)
        righthook = list(righthook)
        righthook[1] = subs(righthook[1],l,p-1-l-s)
        righthook = tuple(righthook)
        return righthook
    if (sumflag == 2):
        if (sflag == 0):
            return lefthook(iflag,0,-eflag)
        if (sflag == 1):
            return lefthook(iflag,0,eflag)

# Returns the left and right hooks
def hooks(iflag,lflag, sflag, eflag, sumflag):
    return lefthook(iflag,lflag,eflag),righthook(iflag,lflag,sflag,eflag,sumflag)

# Get's the left and right hooks, and determines which is
# the hook corresponding to a V_{a,p-1} boundary
def determine_hook(iflag, lflag, sflag, eflag,sumflag):
    lefthook,righthook = hooks(iflag,lflag,sflag,eflag,sumflag)
    if(iflag == 0):
        if (lefthook == (0,0,2,-1)):
            return lefthook
        if (collapse_expression_F2(lefthook[1],lflag,sflag) == 0 and not(lefthook[1] == 0)):
            return lefthook
        return righthook
    if(iflag == 1):
        if (lefthook == (1,0,1,-1)):
            return lefthook
        if (collapse_expression_F2(lefthook[1],lflag,sflag) == 1 and not(lefthook[1] == 0)):
            return lefthook
        return righthook


# Given the choice of paramters, we find the corresponding V_{a,p-1} hook,
# and computes the a value
def filltable(iflag, lflag, sflag, eflag,sumflag):
    hook = determine_hook(iflag,lflag,sflag,eflag,sumflag)
    hook = list(hook)
    if (lflag == 0):
        hook[1] = subs(hook[1],l,0)
    a_value = 0
    if (hook == [0,0,2,-1]):
       return tuple(hook),0
    if (hook == [1,0,1,-1]):
        return tuple(hook),1
    if (hook[3] == 1):
        eq = a-2-hook[1]
        a_value = solve(eq,a)[0]
    if (hook[3] == -1):
        eq = p-1-a-hook[1]
        a_value = solve(eq,a)[0]
    return tuple(hook),a_value

'''
Now we wish to enumerate the results and collect similar cases.
'''
cases = {} # dictionary, keys = hooks, values = list containing the cases that apply
a_value_dic = {} # contains the a_value for the associated hook
count = 0
for iflag in range(0,2):
    for lflag in range(0,3):
        for sflag in range(0,2):
            for eflag in range(0,2):
                for sumflag in range(0,3):
                    if (eflag == 0):
                        eflag = -1
                    if not (sumflag == 2 and (lflag + sflag) % 2 == 1): # Filter out the obvious impossible cases
                        hook,a_value = filltable(iflag,lflag,sflag,eflag,sumflag)
                        # If the hook isn't one of these cases that ends up being
                        # dependent on i, we want to keep the general symbol i
                        if not (hook == (0,0,2,-1) or hook == (1,0,1,-1)):
                            hook = list(hook)
                            hook[0] = i
                            hook = tuple(hook)
                        # We add the case just enumerated to the dictionary and
                        # ensure we have recorded the a value
                        if hook in cases:
                            cases[hook].append((iflag,lflag,sflag,eflag,sumflag))
                        else:
                            cases[hook] = [(iflag,lflag,sflag,eflag,sumflag)]
                            a_value_dic[hook] = a_value

# Complete, now we wish to simplify
for hook in list(cases.keys()):
    cases_for_hook = cases[hook]
    for iterations in range(0,1000): # how many simplification iterations to perform
        # Time to perform the simplifications
        for i in range(0,len(cases_for_hook)):
            for j in range(0,len(cases_for_hook)):
                if not (i == j): # make sure we don't do a self comparison
                    try:
                        case1 = cases_for_hook[i]
                        case2 = cases_for_hook[j]
                        differencespot = onedifference(case1,case2)
                        if not (differencespot == -1 or differencespot == 1 or differencespot == 4): # there is only one difference
                            case1 = list(case1)
                            case1[differencespot] = '-'
                            case1 = tuple(case1)
                            cases_for_hook.remove(case2)
                            cases_for_hook[i] = case1
                    except:
                        # Time to just reiterate
                        None
    cases[hook] = cases_for_hook

# Now that we've simplified, we go through and print out table friendly
# data.
for hook in list(cases.keys()):
    cases_for_hook = cases[hook]
    print("Hook: " + str(hook) + ", cases: " + str(cases_for_hook) + ", a: " + str(a_value_dic[hook]))
