from recursion_formula import *
import sympy as sp

# Define a function which eliminates elements with too large support and finds relations
def relations(d,sorted_monoms,use_cycles = True,n_max=None):
    '''
    Returns the relations between monomials of generators for each norm
    '''
    
    # Create a string of allowed symbols in monomial names
    if use_cycles:
        allowed_generators = alphabet[0:d//2] + "^1234567890" #Generators that vanish in lower dimension d.
    else:
        allowed_generators = alphabet[0:d//2] + "^1234567890" #Generators that vanish in lower dimension d.
    
    if n_max is None:
        n_max = d+1
    
    output = []
    for n in range(2,n_max):
        key_names = set([])
        monoms = []

        # Get the relevant monomials and basis elements in A(d) of norm n.
        for name, dic in sorted_monoms[n]:
            ## Consider only monomials of _non-zero_ generators:
            monom_non_triv = True
            for char in name:
                if char not in allowed_generators:
                    monom_non_triv = False
            if monom_non_triv:
                monoms.append((name,dic))
                # Take all summands with valid support
                for element in dic.keys():
                    if magn(element) <= d: #Sort out elements with support equal to d.
                        key_names.add(element)
        key_names = sorted(list(key_names))
        monoms = sorted(monoms)
        
        # Create a matrix with coefficients from the monomial equations
        A = np.zeros((len(key_names),len(monoms)), dtype = int)

        for i, key in enumerate(key_names):
            for j, equation in enumerate(monoms):
                if key in (equation[1]).keys():
                    A[i,j] = (equation[1])[key]

        output.append(([m[0] for m in monoms], sp.Matrix(A).nullspace()))
    return output



def relations_to_latex(rel,names,n,add_commas=False):
    lines = []
    lines.append("Norm "+str(n+2))
    lines.append(r"\[ \begin{bmatrix}")
    s = ""
    for x in names:
        s += format_monomname(x) + " & "
    lines.append(s[:-2] + r"\\")
    for row in rel:
        rel = ""
        common_denom = sp.ilcm(*tuple(x.q for x in row))
        for x in common_denom * row:
            rel += str(x) + add_commas*"," + " & "
        lines.append(rel[:-2] + r"\\")
    lines.append(r"\end{bmatrix} \]")
    
    return lines


def write_relations_to_file(d,sorted_monoms, use_cycles = False, n_max=None, add_commas = False):
    L = relations(d,sorted_monoms, use_cycles, n_max)
    lines = []
    for n, (names, null_space) in enumerate(L):
        if len(null_space) and n+2!=d:
            lines += relations_to_latex(null_space,names,n,add_commas)
        
        elif n+2==d:
            lines.append("monomials in norm " + str(n+2))
            lines.append(r"\[")
            lines.append("".join(format_monomname(x)+", " for x in names)[:-2])
            lines.append(r"\]")
    return lines



def compute_minimal_relations(d,sorted_monoms,extra_relations=True,
                                use_cycles=True,return_relations=False):
    
    n_max = 3*d//2 if extra_relations else d+1
    sufficient_relations = relations(d,sorted_monoms,use_cycles,n_max)[(d-1)//2:]
    
    name_tables = []
    lifted_relations = []
    for names,A in sufficient_relations:
        #sp.pprint(A)
        name_tables.append({name:i for i,name in enumerate(names)})
        lifted_relations.append(sp.Matrix(len(names),0,[]))
    
    for i,char in enumerate(alphabet):
        for j, (names_lower, rels) in enumerate(sufficient_relations):
            if i+j+1 < len(sufficient_relations):
                names_higher = name_tables[i+j+1]
                conversion_matrix = sp.zeros(len(names_higher),len(names_lower))

                for x,name in enumerate(names_lower):
                    new_name = "".join(sorted(name+char))
                    if new_name in names_higher:
                        conversion_matrix[names_higher[new_name],x] = 1
                
                rels_lifted = []
                for rel in rels:
                    rels_lifted.append(conversion_matrix*rel)

                lifted_relations[i+j+1] = lifted_relations[i+j+1].row_join(sp.Matrix([rels_lifted]))
    
def compute_minimal_relations(d,sorted_monoms,extra_relations=True,use_cycles=True,return_relations=False):
    n_max = 3*d//2 if extra_relations else d+1
    sufficient_relations = relations(d,sorted_monoms,use_cycles,n_max)[(d-1)//2:]
    
    name_tables = []
    lifted_relations = []
    for names,A in sufficient_relations:
        #sp.pprint(A)
        name_tables.append({name:i for i,name in enumerate(names)})
        lifted_relations.append(sp.Matrix(len(names),0,[]))
    
    for i,char in enumerate(alphabet):
        for j, (names_lower, rels) in enumerate(sufficient_relations):
            if i+j+1 < len(sufficient_relations):
                names_higher = name_tables[i+j+1]
                conversion_matrix = sp.zeros(len(names_higher),len(names_lower))

                for x,name in enumerate(names_lower):
                    new_name = "".join(sorted(name+char))
                    if new_name in names_higher:
                        conversion_matrix[names_higher[new_name],x] = 1
                
                rels_lifted = []
                for rel in rels:
                    rels_lifted.append(conversion_matrix*rel)

                lifted_relations[i+j+1] = lifted_relations[i+j+1].row_join(sp.Matrix([rels_lifted]))
    
    L = []
    for (names,A),(n,B) in zip(sufficient_relations,enumerate(lifted_relations)):
        if return_relations:
            
            if B.shape[1]>0:
                AB = sp.Matrix(sp.BlockMatrix((B,sp.BlockMatrix(A))))
                pivots = AB.rref()[1]
                min_rels = [AB[:,pivot] for pivot in pivots if pivot>=B.shape[1]]
                L += relations_to_latex(min_rels,names,(d-1)//2 + n)
            else:
                L += relations_to_latex(A,names,(d-1)//2 + n)
        else:
            L.append(len(A)-B.rank())

    return L
