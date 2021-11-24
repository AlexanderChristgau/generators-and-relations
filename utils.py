import numpy as np
from math import factorial as fact


alphabet = "abcdefghijklmnopq" 


def norm(a):
    return ((1+np.arange(len(a)))*a).sum()

def card(a):
    return a.sum()

def magn(a):
    return ((2+np.arange(len(a)))*a).sum()



def partial(i,a):
    if i==2:
        return a - np.eye(1,len(a),0,dtype=int)[0]
    return a + np.eye(1,len(a),i-3,dtype=int)[0] - np.eye(1,len(a),i-2,dtype=int)[0]



def partialterm(l,i,b):
    if i==2:
        return l-magn(b)+1
    return b[i-3]+1


def accel_asc(n):
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



# Turns monomial "aa" into "a^2" etc.
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
def get_norm(monomname):
    answer = 0
    D = {c:(i+1) for i,c in enumerate(alphabet)}
    for char in monomname:
        answer += D[char]
    return answer
