from sympy import *
from functools import reduce

X = [X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13] = (
    symbols('X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12, X13')
)

def gen_four_body_ham(gen_stabs, sigmas):
    # H = - J (Σ syndrome * stab) - h (Σ sigma_i)
    syndromes = symbols(' '.join([f's{i}' for i in range(1, len(gen_stabs) + 1)]))
    J, h = symbols('J h')
    H1 = -J * reduce(lambda x, y: x + y, [syn * stab for syn, stab in zip(syndromes, gen_stabs)])
    H2 = -h * reduce(lambda x, y: x + y, sigmas)
    return H1 + H2

def find_degree_of_expr(expr):
    # Fully expand to get the terms in the expression
    terms = expand(expr).args

    max_terms_multiplied = 0
    max_term_example = None

    for term in terms:
        # check if the term is a multiplication
        if isinstance(term, Mul):
            # get the number of factors in the multiplication
            num_factors = len(term.args)
            # update the maximum number of terms multiplied
            if num_factors > max_terms_multiplied:
                max_terms_multiplied = num_factors
                max_term_example = term

    print(max_terms_multiplied, max_term_example)

stab_gen = [X1*X4*X6, X2*X4*X5*X7, X3*X5*X8, X6*X9*X11, X7*X9*X10*X12, X8*X10*X13]

H = gen_four_body_ham(stab_gen, X)
print(H)
find_degree_of_expr(H)

