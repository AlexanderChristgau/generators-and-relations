import numpy as np
from math import factorial as fact
from functools import lru_cache
from itertools import product
from tqdm import tqdm
import pickle
from utils import * #Contains relevant functions for A(d) such as norm and cardinality

arr = np.array

'''
This script contains a series of functions used for computing monomials of generators in terms of basis elements in the algebra A(d).
'''

@lru_cache(maxsize=1000000) # use memoization since small values will be called repeatedly.
def Lambda(l: int,a: tuple,b: tuple):
    '''
    Computes normalized coefficient of l-cycle in the product a*b based on recursion.
    '''
    a_arr, b_arr = arr(a),arr(b) # convert to arrays in order to apply standard functions
    
    # base case
    if card(a_arr)==0 or card(b_arr)==0 or l<2:
        return fact(l-1)
    
    # recursion formula
    if norm(a_arr)+norm(b_arr)==l-1:
        if magn(a_arr)<l and magn(b_arr)<=l:
            answer = 0
            for i in range(2,min(l+1,len(b)+2)):
                if b[i-2] > 0:
                    answer += (i-1)*partialterm(l,i,b_arr)*Lambda(l-1,a,tuple(partial(i,b)))
            return l*answer // (l-magn(a_arr))

        elif magn(a_arr)<=l and magn(b_arr)<l:
            answer = 0
            for i in range(2,min(l+1,len(a)+2)):
                if a[i-2] > 0:
                    answer += (i-1)*partialterm(l,i,a_arr)*Lambda(l-1,b,tuple(partial(i,a)))
            return l*answer//(l-magn(b_arr))
    return 0



# non-normalized coefficient of l-cycle in the product a*b
cycle_coeff = lambda l,a,b: Lambda(l,a,b)//fact(l-1)



def gen_cdecomp(a: tuple,c: tuple):
    '''
    Generates all c-decompositions of a as a list of pairs (A,A').
    '''
    answer = []
    a_arr = arr(a)
    nu = np.nonzero(c)[0][0] +2 #smallest non-zero index
    
    # Create a generator for possible values of A:
    intervals = [range(0,min(a[i],nu)*(i<nu)+1) for i in range(0,len(a))]
    cart_prod = product(*intervals)
    
    
    for A in cart_prod:
        A_arr = arr(A)
        if magn(A_arr)<=nu:
            answer.append((A,tuple(a_arr-A_arr)))
        
    return answer



@lru_cache(maxsize=1000000)
def Lambda2(c,a,b):
    '''
    Computes coefficient of $[c]$ in product $[a]\cdot [b]$
    '''
    if norm(a) + norm(b) != norm(c):
        return 0
    
    if sum(a)==0 and sum(b)==0 and sum(c) == 0:
        return 1

    nu = np.nonzero(c)[0][0] +2
    a_decomps = gen_cdecomp(a,c)
    b_decomps = gen_cdecomp(b,c)
    
    answer = 0
    for A1,A2 in a_decomps:
        for B1,B2 in b_decomps:
            if norm(arr(A1))+norm(arr(B1)) == nu-1:
                answer += cycle_coeff(nu,A1,B1) * Lambda2(tuple(c - np.eye(1,len(c),nu-2,dtype=int)[0]), A2,B2)
    return answer



def multiply_basis(a: tuple,b: tuple):
    '''
    Given two basis elements a and b, compute the product a*b as a dictionary.
    '''
    product = {}
    
    for cycle_type in generate_cycletypes(norm(a)+norm(b)):
        if coeff := Lambda2(cycle_type,a,b):
            product[cycle_type] = coeff
    
    return product



def multiply_generator(generator: tuple, g: dict, d: int=100):
    product = {}
    
    for basis1, coeff1 in g.items():
        for basis2, coeff2 in multiply_basis(generator,basis1).items():
            if magn(basis2) > d:
                pass
            elif basis2 in product.keys():
                product[basis2] += coeff1*coeff2
            else:
                product[basis2] = coeff1*coeff2
    return product


def generate_monomials(d=12,use_cycles=False, max_norm=100):
    # Whether to use cycle generators or transpositions as generators.
    if use_cycles:
        generators = [(c,(0,)*i + (1,)) for i,c in enumerate(alphabet[0:d//2])] #cycle generators
    else:
        generators = [(c,(i+1,)) for i,c in enumerate(alphabet[0:d//2])] #(a_2,0,0,...)
    
    # The list of monomials is initialized as a list of the generators
    monomials = [(i,{j:1}) for i,j in generators]
    uniquenames = set([c for c in alphabet[0:d//2]])
    
    # Repeatedly multiply each generator on list of monomials and add to list of monomials
    for name1, generator in tqdm(generators):
        for name2,g in monomials:
            if norm(list(g.keys())[0]) + norm(generator) <= min(d,max_norm):
                newname = "".join(sorted(name1+name2))
                if newname not in uniquenames:
                    uniquenames.add(newname)
                    prod = multiply_generator(generator, g,d)
                    if prod != {}:
                        monomials.append((newname, prod))

    return monomials

def generate_monomials_sorted(d=12,use_cycles=False,max_norm=100):
    monomials = generate_monomials(d,use_cycles,max_norm)
    sorted_monoms = {i:[] for i in range(1,len(monomials)+1)}
    for name, dic in monomials:
        n = get_norm(name)
        sorted_monoms[n] += [(name, dic)]
    return sorted_monoms


if __name__ == "__main__":
    monomials = generate_monomials_sorted(d = 30, use_cycles=True, max_norm=16)
    with open('monomials_cylces.pkl', 'wb') as file:
        pickle.dump(monomials, file)