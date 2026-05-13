# -*- coding: utf-8 -*-
"""

@author: Denver-James Logan Marchment

This file contains functions related to the modular representation theory
of the group G = SL_2(F_p).

"""

import FBmodule

class FGmodule:
    
    # Returns the left boundary module of a
    # M(i,l,s,epsilon)
    def __leftboundary(self,i,l,s,epsilon):
        H_l = (-1,-1,-1,-1)
        if (s == 1):
            if (l == 0):
                H_l = (i,0,1,-1)
            elif (l > 0):
                H_l = (i,l-1,2,1)
        elif (s > 1):
            if (l == 0 and epsilon == 1):
                H_l = (i,0,1,-1)
            elif (l > 0 and epsilon == 1):
                H_l = (i,l-1,2,1)
            elif (epsilon == -1):
                H_l = (i,l,2,-1)
        return H_l
    
    # Returns the right boundary module of a
    # M(i,l,s,epsilon)
    def __rightboundary(self,i,l,s,epsilon):
        H_r = (-1,-1,-1,-1)
        if (s == 1):
            H_r = (i,l,2,-1)
        elif (s > 1):
            epsilon_2 = epsilon
            if (s % 2 == 0):
                epsilon_2 = -epsilon
            if (l + s > (self.p-1)/2):
                H_r = self.__leftboundary(i, self.p-1-l-s, s, epsilon_2)
            elif (l + s <= (self.p-1)/2):
                if (epsilon_2 == -1):
                    H_r = (i,l+s-2,2,1)
                elif (epsilon_2 == 1):
                    H_r = (i,l+s-1,2,-1)
        return H_r
    
    # Returns both boundaries of a
    # M(i,l,s,epsilon)
    def boundaries(self,i,l,s,epsilon):
        H_l = self.__leftboundary(i,l,s,epsilon)
        H_r = self.__rightboundary(i, l, s, epsilon)
        return H_l,H_r
    
    # Returns both boundaries of a
    # M(i,l,s,epsilon)
    def __boundaries(self,i,l,s,epsilon):
        H_l = self.__leftboundary(i,l,s,epsilon)
        H_r = self.__rightboundary(i, l, s, epsilon)
        return H_l,H_r
    
    # Returns the hook that is given by a
    # V_{a,p-1}
    # When we go to check with boundary is of the form V_{a,p-1}, note that because of
    # way my algorithms compute what the boundaries are, we can directly check if they
    # match as given in (insert corollary number here).
    def __hookofinterest(self,H_l,H_r):
        H = (-1,-1,-1,-1) # specific boundary of interest
        i = H_l[0]
        if (i == 0):
            # First see if H_l is the boundary of interest, if not, it's H_r
            if (H_l == (0,0,2,-1)):
                H = H_l
            elif (H_l[2] > 1 and H_l[1] % 2 == 0 and 0 <= H_l[1] <= (self.p-5)/2 and H_l[3] == 1):
                H = H_l
            elif (H_l[2] > 1 and H_l[1] % 2 == 0 and 0 <= H_l[1] <= (self.p-3)/2 and H_l[3] == -1):
                H = H_l
            else :
                H = H_r
        elif (i == 1):
            if (H_l == (1,0,1,-1)):
                H = H_l
            elif (H_l[2] > 1 and H_l[1] % 2 == 1 and 0 <= H_l[1] <= (self.p-5)/2 and H_l[3] == 1):
                H = H_l
            elif (H_l[2] > 1 and H_l[1] % 2 == 1 and 0 <= H_l[1] <= (self.p-3)/2 and H_l[3] == -1):
                H = H_l
            else :
                H = H_r
        return H
    
    # Function which returns a normalized form of
    # M(i,l,s,epsilon).
    # In short, this simply means that if
    # l + s > (p-1)/2, we always ensure the left
    # boundary module is of the form V_{a,p-1}.
    def __normalise(self,i,l,s,epsilon):
        if (l + s <= (self.p-1)/2):
            return (i,l,s,epsilon)
        else:
            H_l,H_r = self.__boundaries(i,l,s,epsilon)
            if (H_l == self.__hookofinterest(H_l, H_r)):
                return (i,l,s,epsilon)
            else:
                # We need to normalise
                epsilon2 = epsilon
                if (s % 2 == 0):
                    epsilon2 = -epsilon
                return (i,self.p-1-l-s,s,epsilon2)
    
    # This function returns the simple module 
    # represented by the j^th edge in the Brauer tree
    # for G
    # Note if i > (p-1)/3, we loop round.
    def __get_edge(self,i,j):
        if (i % 2 == 0):
            if (j <= (self.p-1)/2):
                if(j % 2 == 1):
                    return j
                elif(j % 2 == 0):
                    return self.p-j
            elif (j > (self.p-1)/2):
                return self.__get_edge(i,self.p-j)
        elif (i % 2 == 1):
            if (j <= (self.p-1)/2):
                if(j % 2 == 0):
                    return j
                elif(j % 2 == 1):
                    return self.p-j
            elif (j > (self.p-1)/2):
                return self.__get_edge(i,self.p-j)
    
    # Function which returns the dimension of a M(i,l,s,epsilon)
    def __get_dimension(self,i,l,s):
        dimension = 0
        for j in range(l+1,l+s+1):
            factor = self.__get_edge(i,j)
            dimension += factor
        return dimension
    
    # Function which returns the composition factors of a
    # M(i,l,s,epsilon)
    def __get_composition_factors(self,i,l,s):
        composition_factors = [0] * (self.p)
        for j in range(l+1,l+s+1):
            factor = self.__get_edge(i,j)
            composition_factors[factor-1] += 1
        return composition_factors
    
    # This function returns the entries of the inverse Cartan matrix
    # for G
    def __Gamma(self,n,m):
        if (n <= m):
            return (-1) **(n+m) * (n - (2*n*m)/self.p)
        elif (n > m):
            return (-1) **(n+m) * (m - (2*n*m)/self.p)
    
    # This function returns the dimension of M(i,l,s,epsilon)
    # modulo p
    def __L(self,i,l,s):
        # l,s both odd
        if (l % 2 == 1 and s % 2 == 1):
            return (1-i)*self.p + (-1)**(i+1) * (2*l+s+1)/2
        # l,s both even
        if (l % 2 == 0 and s % 2 == 0):
            return (1-i)*self.p + (-1)**(i+1) * s/2
        # l odd, s even
        if (l % 2 == 1 and s % 2 == 0):
            return i*self.p + (-1)**(i) * s/2
        # l even, s odd
        if (l % 2 == 0 and s % 2 == 1):
            return i*self.p + (-1)**(i) * (2*l+s+1)/2
    
    # This function takes a M(i,l,s,epsilon) and returns the corresponding
    # U_{a,b} under the Green correspondence
    def __greencor(self,i,l,s,epsilon):
        # b is easy to find, just the dimension modulo p
        b = int(self.__L(i,l,s))
        a = -1
        # a is harder, requires us to look at the boundary modules.
        # We will first compute both the boundary modules.
        # Then, we determine which boundary module is 
        # of the form V_{a,p-1}.
        # Finally, we find the variable 'a' from this.
        H_l,H_r = self.__boundaries(i,l,s,epsilon) # get the boundaries
        H = self.__hookofinterest(H_l, H_r) # specific boundary of interest
        # determined the hook, now we compute 'a'.
        if(H == (0,0,2,-1)):
            a = 0
        elif(H[2] == 1):
            a = 1
        elif (H[3] == 1):
            a = H[1] + 2
        elif (H[3] == -1):
            a = self.p - 1 - H[1]
        return (a,b)
    
    # Returns the composition factors of
    # Res^G_B(V_i)
    def __get_composition_factors_res_simple(self,i):
        a = 0
        if (i == 1):
            a = 0
        else:
            a = self.p-i
        res_factors = (FBmodule.FBmodule(self.p, [ [a,i] ] )).get_composition_factors()
        return res_factors
    
    # Returns the composition factors of a
    # res^G_B(M(i,l,s,epsilon))
    def __get_composition_factors_res(self,i,l,s,epsilon):
        res_factors = [0] * (self.p-1)
        # Get the factors as an FGmodule first
        FGfactors = (FGmodule(self.p, [ [i,l,s,epsilon] ] )).get_composition_factors()
        # Now we go through each one
        for j in range(0,self.p):
            j_actual = j+1
            a = -1
            if (j_actual == 1):
                a = 0
            else:
                a = self.p-j_actual
            res_Vj_factors = (FBmodule.FBmodule(self.p, [ [a,j_actual] ] )).get_composition_factors()
            multiplicity = FGfactors[j]
            for m in range(0,multiplicity):
                res_factors = [x+y for x,y in zip(res_factors,res_Vj_factors)]
        return res_factors
    
    # Returns the composition factors of a
    # res^G_B(P_{V_i})
    def __get_composition_factors_res_proj(self,i):
        res_factors = [0] * (self.p-1)
        if (i == 1):
            res_factors = [x+y for x,y in zip(res_factors,self.__get_composition_factors_res_simple(1))]
            res_factors = [x+y for x,y in zip(res_factors,self.__get_composition_factors_res_simple(1))]
            res_factors = [x+y for x,y in zip(res_factors,self.__get_composition_factors_res_simple(self.p-2))]
        if (i == self.p):
            res_factors = [x+y for x,y in zip(res_factors,self.__get_composition_factors_res_simple(self.p))]
        if (1 < i < self.p):
            res_factors = [x+y for x,y in zip(res_factors,self.__get_composition_factors_res_simple(i))]
            res_factors = [x+y for x,y in zip(res_factors,self.__get_composition_factors_res_simple(i))]
            res_factors = [x+y for x,y in zip(res_factors,self.__get_composition_factors_res_simple(self.p+1-i))]
            res_factors = [x+y for x,y in zip(res_factors,self.__get_composition_factors_res_simple(self.p-1-i))]
        return res_factors
    
    # Standard constructor
    def __init__(self, prime, module_list = None, projective_factors = False):
        self.p = prime
        self.non_projectives = {}
        self.projectives = {}
        # If the projective_factors flag is set to true,
        # this means the 'module_list' variable contains a list of the
        # composition factors of a projective FG module,
        # which we compute the decomposition of and set as our module
        if (projective_factors):
            proj_module_list = [0] * (self.p)
            # V_p occurs as many times as a summand as it does as a composition factor
            proj_module_list[self.p-1] = module_list[self.p-1]
            # Gotta now deal with the projectives belonging to the other two blocks
            principal_block_factors = [0] * int((self.p-1)/2)
            nonprincipal_block_factors = [0] * int((self.p-1)/2)
            principal_block_summands = [0] * int((self.p-1)/2)
            nonprincipal_block_summands = [0] * int((self.p-1)/2)
            # We have to start by sorting the other factors into their two blocks
            for j in range(1,int((self.p+1)/2)):
                if (j % 2 == 1):
                    principal_block_factors[j-1] = module_list[j-1]
                    nonprincipal_block_factors[j-1] = module_list[(self.p-j)-1]
                elif (j % 2 == 0):
                    principal_block_factors[j-1] = module_list[(self.p-j)-1]
                    nonprincipal_block_factors[j-1] = module_list[j-1]
            # Time to do matrix multiplication
            for n in range(1,int((self.p+1)/2)):
                sum1 = 0
                sum2 = 0
                for m in range(1,int((self.p+1)/2)):
                    sum1 += self.__Gamma(n,m) * principal_block_factors[m-1]
                    sum2 += self.__Gamma(n,m) * nonprincipal_block_factors[m-1]
                principal_block_summands[n-1] = int(round(sum1))
                nonprincipal_block_summands[n-1] = int(round(sum2))
            # Done, now we have to add it all back in
            for j in range(1,self.p):
                if (j % 2 == 1):
                    # j is the index^th edge of the tree
                    index = j
                    if (j > (self.p-1)/2):
                        index = self.p - j
                    proj_module_list[j-1] = principal_block_summands[index-1]
                elif (j % 2 == 0):
                    # j is the index^th edge of the tree
                    index = j
                    if (j > (self.p-1)/2):
                        index = self.p - j
                    proj_module_list[j-1] = nonprincipal_block_summands[index-1]
            # Now we simply add in all the things
            for j in range(1,self.p+1):
                m = proj_module_list[j-1]
                self.add( [ [j]  ] * m )
        else:
            # If there is some modules we need to add now,
            # we should add them
            if not(module_list == None):
                for module in module_list:
                    # If projective
                    if (len(module) == 1):
                        j = module[0]
                        if (j in self.projectives):
                            self.projectives[j] += 1
                        else:
                            self.projectives[j] = 1
                    # Otherwise, a non-projective module is being added
                    else:
                        module = self.__normalise(module[0],module[1],module[2],module[3])
                        if (module in self.non_projectives):
                            self.non_projectives[module] += 1
                        else:
                            self.non_projectives[module] = 1
    
    # Return the current decomposition
    def get_decomposition(self):
        return self.p,self.non_projectives,self.projectives
    
    # Function to add more modules
    def add(self, modules):
        for module in modules:
            # Check if non-projective or projective
            if ( len(module) == 1 ):
                # is projective
                j = module[0]
                if ( j in self.projectives ):
                    self.projectives[j] += 1
                else:
                    self.projectives[j] = 1
            else:
                # non-projective
                module = self.__normalise(module[0],module[1],module[2],module[3])
                if (module in self.non_projectives):
                    self.non_projectives[module] += 1
                else:
                    self.non_projectives[module] = 1
    
    # Function to remove modules
    # if you try to remove a module that doesn't
    # exist in the decomposition, we will simply ignore it
    def remove(self,modules):
        for module in modules:
            # Check if non-projective or projective
            if ( len(module) == 1 ):
                # is projective
                j = module[0]
                if ( j in self.projectives ):
                    # If it's about to go to 0, remove it
                    # otherwise, just subtract 1
                    if ( self.projectives[j] == 1 ):
                        self.projectives.pop(j, None)
                    else:
                        self.projectives[j] -= 1
            else:
                # non-projective
                module = self.__normalise(module[0],module[1],module[2],module[3])
                if (module in self.non_projectives):
                    # If it's about to go to 0, remove it
                    # otherwise, just subtract 1
                    if ( self.non_projectives[module] == 1 ):
                        self.non_projectives.pop(module, None)
                    else:
                        self.non_projectives[module] -= 1
    
    # Function to return the dimension of the module
    def get_dimension(self):
        dimension = 0
        # First add up the dimensions contributed by the non-projectives
        for module in self.non_projectives:
            multiplicity = self.non_projectives[module]
            dimension += multiplicity * self.__get_dimension(module[0],module[1],module[2])
        # Next, add up the dimensions contributed by the projectives
        for module in self.projectives:
            multiplicity = self.projectives[module]
            if (module == 1 or module == self.p):
                dimension += self.p * multiplicity
            else:
                dimension += self.p * 2 * multiplicity
        return dimension
    
    # Function to return the composition factors
    def get_composition_factors(self):
        composition_factors = [0] * (self.p)
        # Count up the composition factors from the
        # non-projectives first
        for module in self.non_projectives:
            multiplicity = self.non_projectives[module]
            module_factors = self.__get_composition_factors(module[0], module[1], module[2])
            for m in range(0,multiplicity):
                composition_factors = [x+y for x,y in zip(composition_factors,module_factors)]
        # Now do the projectives
        for module in self.projectives:
            multiplicity = self.projectives[module]
            module_factors = [0] * (self.p)
            if (module == 1):
                module_factors[0] = 2
                module_factors[self.p-2-1] = 1
            elif (module == self.p):
                module_factors[self.p-1] = 1
            else:
                module_factors[module-1] += 2
                module_factors[self.p+1-module-1] += 1
                module_factors[self.p-1-module-1] += 1
            for m in range(0,multiplicity):
                composition_factors = [x+y for x,y in zip(composition_factors,module_factors)]
        return composition_factors
    
    # Function to restrict to a B module
    def res(self):
        non_projective_modules = [] # This will store all of the non-projective modules to add
        projective_composition_factors = [0] * (self.p)
        # quickly add the composition factors from the projective modules
        for item in self.projectives:
            m = self.projectives[item]
            for k in range(0,m):
                projective_composition_factors = [x+y for x,y in zip(projective_composition_factors,self.__get_composition_factors_res_proj(item))]
        for module in self.non_projectives:
            i,l,s,epsilon = module[0],module[1],module[2],module[3]
            a,b = self.__greencor(i,l,s,epsilon)
            green_factors = (FBmodule.FBmodule(self.p, [ [a,b] ])).get_composition_factors()
            res_factors = self.__get_composition_factors_res(i, l, s, epsilon)
            leftover_factors = [x-y for x,y in zip(res_factors,green_factors)]
            multiplicity = self.non_projectives[module]
            for m in range(0,multiplicity):
                projective_composition_factors = [x+y for x,y in zip(projective_composition_factors,leftover_factors)]
                non_projective_modules.append([a,b])
        restricted_module = FBmodule.FBmodule(self.p,projective_composition_factors,True)
        restricted_module.add(non_projective_modules)
        return restricted_module