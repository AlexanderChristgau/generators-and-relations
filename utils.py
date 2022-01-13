import numpy as np

alphabet = "abcdefghijklmnopq" 

#Define norm, cardinality and magnitude of a tuple
def norm(a):
    return ((1+np.arange(len(a)))*a).sum()

def card(a):
    return a.sum()

def magn(a):
    return ((2+np.arange(len(a)))*a).sum()



#Define the partial opeartor used in the recursion formula
def partial(i,a):
    if i==2:
        return a - np.eye(1,len(a),0,dtype=int)[0]
    return a + np.eye(1,len(a),i-3,dtype=int)[0] - np.eye(1,len(a),i-2,dtype=int)[0]


def partialterm(l,i,b):
    if i==2:
        return l-magn(b)+1
    return b[i-3]+1


#Define functions for generating all cycle types
def accel_asc(n):
    '''
    Outputs a generator of all partitions of norm n.
    See also https://arxiv.org/abs/0909.2331
    '''
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]


def to_ctype(l):
    answer = np.zeros(max(l),dtype=int)
    for i in l:
        answer[i-1] += 1
    return tuple(answer)


def generate_cycletypes(d):
    return (to_ctype(partition) for partition in accel_asc(d))



# Format monomial "aa" into "a^2" etc.
def format_monomname(string):
    counts = {c:0 for c in alphabet}
    for char in string:
        counts[char] +=1
    s = ""
    for char in alphabet:
        if counts[char] == 1:
            s += char
        elif counts[char]>1:
            s += char + "^{" + str(counts[char]) + "}"
    return s



# Computes norm from label in order to sort monomials
D1 = {c:(i+1) for i,c in enumerate(alphabet)}
def get_norm(monomname):
    return sum(D1[char] for char in monomname)