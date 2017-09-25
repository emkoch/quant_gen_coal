from collections import Counter
import sys
import itertools
import sympy
from geneMGF import *

def make_subscr(branch_comb):
    branch_strs = ['.'.join([str(bb) for bb in bc]) for bc in branch_comb]
    return sympy.symbols('t_' + 'x'.join(branch_strs))

def make_subscr_ex(branch_comb):
    indivs = list(indiv for branch in branch_comb for indiv in branch)
    unique_indivs = list(set(indivs))
    ind_mults = Counter(indivs)
    ordering = sorted(unique_indivs, key=ind_mults.get, reverse=True)
    alt_branch_comb = [sorted([ordering.index(bb) for bb in bc]) for bc in branch_comb]
    branch_strs = ['.'.join([str(bb) for bb in bc]) for bc in alt_branch_comb]
    return sympy.symbols('t_' + 'x'.join(branch_strs))

class mgfApprox(object):
    def __init__(self, nIndiv, mOrd):
        self.n_indivs = nIndiv
        self.moment_order = mOrd
        self.approx = self.create_mgf_expr()
        self.mgf_type = 'taylor'

    def create_mgf_expr(self):
        indivs = range(self.n_indivs)
        branches = []
        for branch_size in range(1, self.n_indivs):
            branches += itertools.combinations(indivs, branch_size)
        mgf_base = sympy.Integer(1)
        num_loci = sympy.symbols("L", integer=True)
        for mom_order in range(1, self.moment_order + 1):
            m_max = self.moment_order - mom_order + 1
            branch_combos = itertools.combinations_with_replacement(branches, mom_order)
            for branch_combo in branch_combos:
                branch_mult = list(Counter(branch_combo).values())
                term_factor = sympy.Integer(1)
                for multi in branch_mult:
                    term_factor *= sympy.factorial(multi)
                bc_term = (make_subscr(branch_combo)*
                           sympy.symbols('theta')**mom_order*term_factor**-1)
                for branch in branch_combo:
                    branch_sum = sympy.Integer(0)
                    for indiv in branch:
                        branch_sum += sympy.symbols('k_' + str(indiv))
                    branch_term = sympy.Integer(0)
                    for mut_order in range(1, m_max + 1):
                        branch_term += (sympy.symbols('m_' + str(mut_order))/
                                        sympy.factorial(mut_order)*branch_sum**mut_order)
                    bc_term *= branch_term
                mgf_base += bc_term
        return mgf_base**num_loci

## True full can only happen if we're willing to specify the exact for
## of both the gene mgf and the mutational mgf, this isn't implemented for
## the time being
class mgfApproxFull(mgfApprox):
    def __init__(self, nIndiv, mOrd, gMGF):
        self.gene_mgf = gMGF
        mgfApprox.__init__(self, nIndiv, mOrd)
        self.mgf_type = 'full'

    def create_mgf_expr(self):
        return (self.gene_mgf.gene_mgf_to_trait(self.moment_order)**
                sympy.symbols("L", integer=True))

class mgfApproxLMR(mgfApprox):
    def __init__(self, nIndiv, mOrd):
        mgfApprox.__init__(self, nIndiv, mOrd)
        self.mgf_type = 'lmr'

    def create_mgf_expr(self):
        indivs = range(self.n_indivs)
        branches = []
        for branch_size in range(1, self.n_indivs):
            branches += itertools.combinations(indivs, branch_size)
        mgf_base = sympy.Integer(1)
        num_loci = sympy.symbols("L", integer=True)
        for branch in branches:
            branch_sum = sympy.Integer(0)
            for indiv in branch:
                branch_sum += sympy.symbols('k_' + str(indiv))
            mut_term = sympy.Integer(0)
            for mut_order in range(1, self.moment_order + 1):
                mut_term += (sympy.symbols('m_' + str(mut_order))/sympy.factorial(mut_order)*
                             branch_sum**mut_order)
            mgf_base += make_subscr([branch])*sympy.symbols('theta')*mut_term
        return mgf_base**num_loci

class mgfApproxExchange(mgfApprox):
    def __init__(self, nIndiv, mOrd):
        mgfApprox.__init__(self, nIndiv, mOrd)
        self.mgf_type = "exchangeable"

    def create_mgf_expr(self):
        indivs = range(self.n_indivs)
        branches = []
        for branch_size in range(1, self.n_indivs):
            branches += itertools.combinations(indivs, branch_size)
        mgf_base = sympy.Integer(1)
        num_loci = sympy.symbols("L", integer=True)
        for mut_order in range(1, self.moment_order + 1):
            m_max = self.moment_order - mut_order + 1
            branch_combos = itertools.combinations_with_replacement(branches, mut_order)
            for branch_combo in branch_combos:
                branch_mult = list(Counter(branch_combo).values())
                term_factor = sympy.Integer(1)
                for multi in branch_mult:
                    term_factor *= sympy.factorial(multi)
                bc_term = (make_subscr_ex(branch_combo)*
                           sympy.symbols('theta')**mut_order*term_factor**-1)
                for branch in branch_combo:
                    branch_sum = sympy.Integer(0)
                    for indiv in branch:
                        branch_sum += sympy.symbols('k_' + str(indiv))
                    branch_term = sympy.Integer(0)
                    for mut_order in range(1, m_max + 1):
                        branch_term += (sympy.symbols('m_' + str(mut_order))/
                                        sympy.factorial(mut_order)*branch_sum**mut_order)
                    bc_term *= branch_term
                mgf_base += bc_term
        return mgf_base**num_loci

class traitMGF(object):
    def __init__(self, nIndiv):
        self.num_indiv = nIndiv
        self.mgf = {}
        self.mgf['taylor'] = {}
        self.mgf['lmr'] = {}
        self.mgf['full'] = {}
        self.mgf['exchangeable'] = {}
        self.moments = {}
        self.moments['taylor'] = {}
        self.moments['lmr'] = {}
        self.moments['full'] = {}
        self.moments['exchangeable'] = {}

    def make_mgf(self, order, approx_type = 'taylor', gene_mgf = None):
        if order in self.mgf[approx_type].keys():
            print("mgf already created!")
        else:
            if approx_type == 'taylor':
                self.mgf[approx_type][order] = mgfApprox(self.num_indiv, order)
            elif approx_type == 'lmr':
                self.mgf[approx_type][order] = mgfApproxLMR(self.num_indiv, order)
            elif approx_type == 'full':
                self.mgf[approx_type][order] = mgfApproxFull(self.num_indiv, order, gene_mgf)
            else:
                self.mgf[approx_type][order] = mgfApproxExchange(self.num_indiv, order)

    def calc_moment(self, pows, approx_type = 'taylor', gene_mgf = None):
        """ Calculate a moment specifying an approximation type.
        Arguments:
        pows        -- list of powers (ints) for trait values in each individual
        approx_type -- what approximation to use for calculating the moment
        gene_mgf    -- geneMGF object necessary if using 'full' approximation
        """
        mom_hash = '.'.join([str(power) for power in pows])
        if mom_hash in self.moments[approx_type].keys():
            return self.moments[approx_type][mom_hash]
        print('making mgf approx...')
        self.make_mgf(sum(pows), approx_type, gene_mgf)
        d_mgf = self.mgf[approx_type][sum(pows)].approx
        dummies_nonzero = sympy.symbols(['k_' + str(ii) for ii in range(self.num_indiv)
                                         if pows[ii] > 0])
        dummies_zero = sympy.symbols(['k_' + str(ii) for ii in range(self.num_indiv)
                                      if pows[ii] == 0])
        dummies = sympy.symbols(['k_' + str(ii) for ii in range(self.num_indiv)])
        print('simplifying mgf...')
        d_mgf = d_mgf.subs([(kk, 0) for kk in dummies_zero])
        for indiv in range(self.num_indiv):
            if pows[indiv] > 0:
                for power in range(pows[indiv]):
                    print('taking derivative ' + str(dummies[indiv]) + ' order ' + str(power + 1))
                    d_mgf = sympy.diff(d_mgf, dummies[indiv], 1)
        result = d_mgf.subs([(kk, 0) for kk in dummies_nonzero])
        self.moments[approx_type][mom_hash] = result
        return result

    def mom_lmr(self, pows, approx_type = 'taylor', gene_mgf = None):
        mom_hash = '.'.join([str(power) for power in pows])
        if mom_hash not in self.moments[approx_type].keys():
            self.calc_moment(pows, approx_type, gene_mgf)
        #if approx_type == 'lmr':
        #    return self.moments[approx_type][mom_hash]
        mom_full = self.moments[approx_type][mom_hash]
        num_loci = sympy.symbols('L', integer=True)
        theta = sympy.symbols('theta')
        result = sympy.expand(mom_full).subs(num_loci*theta, sympy.symbols('M'))
        result = result.subs(theta, 0).subs(sympy.symbols('M'), theta*num_loci)
        return result
