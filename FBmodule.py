# -*- coding: utf-8 -*-
"""

@author: Denver-James Logan Marchment

This file contains functions related to the modular representation theory
of the group G = SL_2(F_p).

"""

import FGmodule

class FBmodule:
    
    # Constructor
    def __init__(self, prime, module_list = None, projective_factors = False):
        self.p = prime
        self.non_projectives = {}
        self.projectives = {}
        # If the module to be constructed is based on
        # projective factors
        if (projective_factors):
            # break up the factors into blocks
            principal_proj_factors = []
            nonprincipal_proj_factors = []
            principal_block_summands = []
            nonprincipal_block_summands = []
            for a in range(0,self.p-1):
                if (a % 2 == 0):
                    principal_proj_factors.append(module_list[a])
                elif(a % 2 == 1):
                    nonprincipal_proj_factors.append(module_list[a])
            # Time to do matrix multiplication
            for n in range(1,int((self.p+1)/2)):
                sum1 = 0
                sum2 = 0
                for m in range(1,int((self.p+1)/2)):
                    sum1 += self.__gamma(n,m) * principal_proj_factors[m-1]
                    sum2 += self.__gamma(n,m) * nonprincipal_proj_factors[m-1]
                principal_block_summands.append(int(round(sum1)))
                nonprincipal_block_summands.append(int(round(sum2)))
            # Time to add in the summands
            for j in range(0,int((self.p-1)/2)):
                principal_a = 2 * j
                nonprincipal_a = 2*j+1
                principal_multi = principal_block_summands[j]
                nonprincipal_multi = nonprincipal_block_summands[j]
                self.add([ [principal_a,self.p]  ] * principal_multi)
                self.add([ [nonprincipal_a,self.p]  ] * nonprincipal_multi)
        else:
            # If there is some modules we need to add now,
            # we should add them
            if not(module_list == None):
                for module in module_list:
                    a,b = module[0],module[1]
                    # Check whether it's a projective one
                    # or not
                    if (b <= self.p-1):
                        if ((a,b) in self.non_projectives):
                            self.non_projectives[(a,b)] += 1
                        else:
                            self.non_projectives[(a,b)] = 1
                    elif (b == self.p):
                        if ((a,b) in self.projectives):
                            self.projectives[(a,b)] += 1
                        else:
                            self.projectives[(a,b)] = 1
    
    # This is just the inverse of the Cartan matrix for a block of
    # FB
    def __gamma(self,m,n):
        if (n == m):
            return (self.p-2)/(self.p)
        else:
            return -2/self.p
    
    # Returns the composition factors of a
    # U_{a,b}
    def __get_composition_factors(self,a,b):
        composition_factors = [0] * (self.p-1)
        for j in range(0,b):
            next_factor = (a + 2*j) % (self.p-1)
            composition_factors[next_factor] += 1
        return composition_factors
    
    # Returns the composition factors of a
    # ind^G_B(U_{a,b})
    def __get_composition_factors_ind(self,a,b):
        composition_factors = [0] * (self.p)
        for j in range(0,b):
            next_factor = (a + 2*j) % (self.p-1)
            composition_factors[next_factor] += 1
            composition_factors[self.p - next_factor - 1] += 1
        return composition_factors
    
    # Function which returns the Green correspondent of
    # a U_{a,b} (b <= p-1)
    # All the modules returned by this are already
    # normalised
    def __greencor(self,a,b):
        i = a % 2
        if (0 <= a <= 1):
            if (1 <= b <= (self.p-1)/2):
                return (i,0,2*b + i-1,2*i-1)
            elif ((self.p+1)/2 <= b <= self.p-1):
                return (i,0,2*(self.p-b)-i,2*i-1)
        elif (2 <= a <= (self.p-1)/2):
            if (1 <= b <= (self.p-a)/2):
                return (i,a-1,2*b,1)
            elif((self.p-a+1)/2 <= b <= self.p-a):
                return (i, a - 1, 2*(self.p - a - b) + 1, 1)
            elif(self.p-a <= b <= self.p - (a+1)/2):
                return  (i, 2*(self.p - b) - a - 1, 2*(a + b - self.p) + 1, -1)
            elif (self.p-a/2 <= b <= self.p-1):
                return (i, a + 2*(b - self.p), 2*(self.p - b), 1)
        elif ((self.p+1)/2 <= a <= self.p-2):
            if (1 <= b <= (self.p-a)/2):
                return (i, self.p - a - 2*b, 2*b, -1)
            elif((self.p-a+1)/2 <= b <= self.p-a):
                return (i, a + 2*b - self.p - 1, 2*(self.p - a - b) + 1, 1)
            elif(self.p-a <= b <= self.p - (a+1)/2):
                return (i, self.p - 1 - a, 2*(a + b - self.p) + 1, -1)
            elif (self.p-a/2 <= b <= self.p-1):
                return (i, self.p - 1 - a, 2*(self.p - b), -1)
    
    # Get the current module information
    def get_decomposition(self):
        return self.p,self.non_projectives,self.projectives
    
    # Add more modules
    def add(self,module_list):
        for module in module_list:
            a,b = module[0],module[1]
            # Check whether it's a projective one
            # or not
            if (b <= self.p-1):
                if ((a,b) in self.non_projectives):
                    self.non_projectives[(a,b)] += 1
                else:
                    self.non_projectives[(a,b)] = 1
            elif (b == self.p):
                if ((a,b) in self.projectives):
                    self.projectives[(a,b)] += 1
                else:
                    self.projectives[(a,b)] = 1
    
    # Function to remove modules.
    # Note that, if you try to remove a module which
    # does not exist in the current decomposition, it will
    # simply be ignored.
    def remove(self,module_list):
        for module in module_list:
            a,b = module[0],module[1]
            # Check whether it's a projective one
            # or not
            if (b <= self.p-1):
                if ((a,b) in self.non_projectives):
                    if (self.non_projectives[(a,b)] == 1):
                        self.non_projectives.pop((a,b),None)
                    else:
                        self.non_projectives[(a,b)] -= 1
            elif (b == self.p):
                if ((a,b) in self.projectives):
                    if (self.projectives[(a,b)] == 1):
                        self.projectives.pop((a,b),None)
                    else:
                        self.projectives[(a,b)] -= 1
    
    # Returns the dimension of the module
    def get_dimension(self):
        dimension = 0
        for module in self.non_projectives:
            dimension += self.non_projectives[module] * module[1]
        for module in self.projectives:
            dimension += self.projectives[module] * module[1]
        return dimension
    
    # Returns the composition factors of the module
    def get_composition_factors(self):
        composition_factors = [0] * (self.p-1)
        # First deal with non-projectives
        for module in self.non_projectives:
            multiplicity = self.non_projectives[module]
            a,b = module[0],module[1]
            module_factors = self.__get_composition_factors(a,b)
            for m in range(0,multiplicity):
                composition_factors = [x + y for x,y in zip(composition_factors, module_factors)]
        # Now the projectives
        for module in self.projectives:
            multiplicity = self.projectives[module]
            a,b = module[0],module[1]
            module_factors = self.__get_composition_factors(a,b)
            for m in range(0,multiplicity):
                composition_factors = [x + y for x,y in zip(composition_factors, module_factors)]
        return composition_factors
    
    # Function to induct our B module to a G module
    def ind(self):
        non_projective_modules = [] # This will store all of the non-projective modules to add
        projective_composition_factors = [0] * (self.p)
        # quickly add the composition factors from the projective modules
        for item in self.projectives:
            m = self.projectives[item]
            for k in range(0,m):
                projective_composition_factors = [x+y for x,y in zip(projective_composition_factors,self.__get_composition_factors_ind(item[0],item[1]))]
        for module in self.non_projectives:
            a,b = module[0],module[1]
            (i,l,s,epsilon) = self.__greencor(a,b)
            green_factors = (FGmodule.FGmodule(self.p,[ [i,l,s,epsilon] ])).get_composition_factors()
            ind_factors = self.__get_composition_factors_ind(a, b)
            leftover_factors = [x-y for x,y in zip(ind_factors,green_factors)]
            multiplicity = self.non_projectives[module]
            for m in range(0,multiplicity):
                projective_composition_factors = [x+y for x,y in zip(projective_composition_factors,leftover_factors)]
                non_projective_modules.append([i,l,s,epsilon])
        induced_module = FGmodule.FGmodule(self.p,projective_composition_factors,True)
        induced_module.add(non_projective_modules)
        return induced_module