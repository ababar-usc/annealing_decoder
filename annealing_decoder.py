from sympy import *
from functools import reduce
import re

X = [X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13] = (
    symbols('X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12, X13')
)

Z = [Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z12, Z13] = (
    symbols('Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9 Z10 Z11 Z12, Z13')
)

def gen_four_body_ham(gen_stabs, sigmas):
    # H = - J (Σ syndrome * stab) - h (Σ sigma_i)
    syndromes = symbols(' '.join([f's{i}' for i in range(1, len(gen_stabs) + 1)]), commutative=False)
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

def find_necessary_substitutions_commutative(terms_needing_substitution):
    subs = set()
    for term in terms_needing_substitution:
        sub1 = (term[0], term[1])
        reversed_sub1 = (term[1], term[0])
        # sub for first two terms if none of them are in subs
        if len(term) == 3:
            sub2 = (term[1], term[2])
            reversed_sub2 = (term[2], term[1])
            if not any(sub in subs for sub in [sub1, sub2, reversed_sub1, reversed_sub2]):
                subs.add(sub1)
        # add both sub1 and sub2 for 4 terms
        if len(term) == 4:
            sub2 = (term[2], term[3])
            reversed_sub2 = (term[3], term[2])
            if not reversed_sub1 in subs:
                subs.add(sub1)
            if not reversed_sub2 in subs:
                subs.add(sub2)
    return subs

def find_necessary_substitutions(terms_needing_substitution):
    subs = set()
    for term in terms_needing_substitution:
        sub1 = (term[0], term[1])
        # sub for first two terms if none of them are in subs
        if len(term) == 3:
            sub2 = (term[1], term[2])
            if not any(sub in subs for sub in [sub1, sub2]):
                subs.add(sub1)
        # add both sub1 and sub2 for 4 terms
        if len(term) == 4:
            sub2 = (term[2], term[3])
            subs.update([sub1, sub2])
    return subs

def check_if_hamiltonian_2_body(H):
    for term in expand(H).args:
        relevant_symbols = [symbol for symbol in term.free_symbols if re.match(r'z|x', str(symbol))]
        if len(relevant_symbols) > 2:
            print(relevant_symbols, term)
            return False
    return True

def generate_penalty_hamiltonian(subs, z_terms):
    penalty_term = 0
    for (sub, z_term) in zip(subs, z_terms):
        penalty_term += sub[0] * sub[1] - 2 * z_term * (sub[0] + sub[1]) + 3 * z_term
    alpha = symbols('a')
    return alpha * penalty_term

def convert_to_qubo(H):
    # create binary variables for each ising spin
    binary_vars = symbols(' '.join([f'x{i}' for i in range(1, len(Z)+1)]), commutative=False)
    # replace each ising spin with binary variable
    binary_vars_subs = map(lambda bin_var: 1-2*bin_var, binary_vars)
    H_after_binary_vars = H.subs(dict(zip(Z, binary_vars_subs)))

    # replace with z_i
    terms_needing_substitution = []
    for term in expand(H_after_binary_vars).args:
        x_symbols_in_term = [symbol for symbol in term.free_symbols if 'x' in str(symbol)]
        x_symbols_strings_ordered = [symbol for symbol in str(term).split('*') if 'x' in symbol]
        x_symbols_in_term_ordered = sorted(x_symbols_in_term, key=lambda x: x_symbols_strings_ordered.index(str(x)))
        if len(x_symbols_in_term_ordered) > 2:
            # print(term.args, term.free_symbols)
            terms_needing_substitution.append(x_symbols_in_term_ordered)

    subs = find_necessary_substitutions(terms_needing_substitution)
    z_terms = symbols(' '.join([f'z{i}' for i in range(1, len(subs)+1)]), commutative=False)
    z_term_subs = list(map(lambda x_pair: x_pair[0] * x_pair[1], subs))
    H_after_z_subs = expand(H_after_binary_vars).subs(dict(zip(z_term_subs, z_terms)))
    # print(check_if_hamiltonian_2_body(H_after_z_subs))

    # create penalty hamiltonian
    penalty_ham = generate_penalty_hamiltonian(subs, z_terms)
    # return H + H penalty
    qubo_H = H_after_z_subs + expand(penalty_ham)

    # replace with y_i
    y_vars = symbols(' '.join([f'y{i}' for i in range(1, len(binary_vars) + len(z_terms) + 1)]), commutative=False)
    qubo_H = qubo_H.subs(dict(zip(binary_vars + z_terms, y_vars)))
    return qubo_H

def convert_to_dwave_params(qubo_H, J, h, alpha, syndromes):
    dwave_params = {}
    for term in expand(qubo_H).args:
        y_symbols_in_term = [symbol for symbol in term.free_symbols if 'y' in str(symbol)]

        if len(y_symbols_in_term) > 2:
            return ValueError("Qubo ham with more than 2 y terms: error in conversion process.")
        # TODO: what do we do with this?
        if len(y_symbols_in_term) == 0:
            print(term)
            continue
        if len(y_symbols_in_term) == 2:
            key = tuple([str(y_symbols_in_term[0]), str(y_symbols_in_term[1])])
        if len(y_symbols_in_term) == 1:
            key = tuple([str(y_symbols_in_term[0]), str(y_symbols_in_term[0])])\
        # make substitutions
        value = term.subs(dict(zip(y_symbols_in_term, [1]*len(y_symbols_in_term))))
        # substitute for J, h, alpha
        value = value.subs({symbols('J'): J, symbols('h'): h, symbols('a'): alpha})
        # substitute for syndromes
        # NB: the line below doesn't work for some reason
        # value = value.subs({symbols(f's{i+1}'): syndromes[i] for i in range(len(syndromes))})
        syndrome_symbols_in_term = [symbol for symbol in value.free_symbols if 's' in str(symbol)]
        syndromes_indexes = [int(re.search(r'\d+', str(s)).group()) for s in syndrome_symbols_in_term]
        syndromes_in_term = [syndromes[i-1] for i in syndromes_indexes]
        value = value.subs(dict(zip(syndrome_symbols_in_term, syndromes_in_term)))
        dwave_params[key] = value
    return dwave_params

stab_gen_x = [X1*X4*X6, X2*X4*X5*X7, X3*X5*X8, X6*X9*X11, X7*X9*X10*X12, X8*X10*X13]
stab_gen_z = [Z1*Z4*Z6, Z2*Z4*Z5*Z7, Z3*Z5*Z8, Z6*Z9*Z11, Z7*Z9*Z10*Z12, Z8*Z10*Z13]

H_x = gen_four_body_ham(stab_gen_x, X)
H_z = gen_four_body_ham(stab_gen_z, Z)

qubo_H = convert_to_qubo(H_z)
# for term in qubo_H.args: print(term)
dwave_params = convert_to_dwave_params(qubo_H, 10, 5, 1, [1, 1, 1, 1, 1, 1])
for key, val in dwave_params.items():
    print(key, val)
