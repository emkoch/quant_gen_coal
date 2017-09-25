from copy import deepcopy
import geneMGF
import sympy
import itertools

def make_partition_term(deme_part, coal_rates, migr_rates):
    result = sympy.Integer(0)
    num_demes = len(coal_rates)
    for deme_idx, deme in enumerate(deme_part):
        n_deme_lineages = len(deme)
        result += sympy.binomial(n_deme_lineages, 2)*coal_rates[deme_idx]
        for other_deme_idx in range(num_demes):
            if deme_idx != other_deme_idx:
                result += n_deme_lineages*migr_rates[deme_idx, other_deme_idx]
        for lineage in deme:
            result -= sympy.symbols('s_' + '.'.join([str(tip) for tip in sorted(lineage)]))
    return result 

def extract_lineages(deme_part):
    """ Extract the individual lineages from a deme partitioning """
    return [lin for deme in deme_part for lin in deme]
            
def phi_swap_list(deme_part_1, deme_part_2):
    """ Create a list for swapping mgf variables. """
    num_demes = len(deme_part_1)
    lineages_1 = extract_lineages(deme_part_1)
    indv_partitions_1 = list(geneMGF.partition(lineages_1))
    indv_partitions_1.remove([[individual for lineage in lineages_1
                                 for individual in lineage]])
    deme_partitions_1 = [deme_part for indv_partition in indv_partitions_1
                           for deme_part in geneMGF.deme_partition(indv_partition, num_demes)]
    
    lineages_2 = extract_lineages(deme_part_2)
    indv_partitions_2 = list(geneMGF.partition(lineages_2))
    indv_partitions_2.remove([[individual for lineage in lineages_2
                                 for individual in lineage]])
    deme_partitions_2 = [deme_part for indv_partition in indv_partitions_2
                           for deme_part in geneMGF.deme_partition(indv_partition, num_demes)]
    result = []
    for ii in range(len(deme_partitions_1)):
        result.append((geneMGF.deme_part_to_symbol(deme_partitions_1[ii]),
                       geneMGF.deme_part_to_symbol(deme_partitions_2[ii])))
    return result

def dummy_swap_list(deme_part_1, deme_part_2):
    """ Create a list for swapping dummy variables in mgf solutions
    Arguments:
    deme_part_1 -- initial deme partition
                   ex: [[[1], [2]], [[3]]]
    deme_part_2 -- deme partition to swap with
                   ex: [[['a'], ['b']], [['c']]]
    Result: A list of lists of length two specifying the swapping
    """
    num_demes = len(deme_part_1)
    lineages_1, lineages_2 = extract_lineages(deme_part_1), extract_lineages(deme_part_2)
    symbols_1, symbols_2 = [], []
    for num_lins in range(1, len(lineages_1)):
        for lins_1 in itertools.combinations(lineages_1, num_lins):
            symbols_1.append("s_" + ".".join([str(lin) for lin in 
                                              sorted(list(itertools.chain.from_iterable(lins_1)))]))
        for lins_2 in itertools.combinations(lineages_2, num_lins):
            symbols_2.append("s_" + ".".join([str(lin) for lin in 
                                              sorted(list(itertools.chain.from_iterable(lins_2)))]))
    result = []
    for idx_sym in range(len(symbols_1)):
        result.append((symbols_1[idx_sym], symbols_2[idx_sym]))
    return result
            
def deme_part_after_coal(deme_part, coal_pair):
    result = deepcopy(deme_part)
    for deme_idx, deme in enumerate(result):
        if coal_pair[0] in deme and coal_pair[1] in deme:
            result[deme_idx].remove(coal_pair[0])
            result[deme_idx].remove(coal_pair[1])
            result[deme_idx].append(sorted(coal_pair[0] + coal_pair[1]))
    return result

def make_sol_key(symbol):
    """ Make a solution key from a given symbol. """
    counts = [str(deme_portion.count("(")) for 
              deme_portion in str(symbol).split("_")[1].split("].[")]
    return ".".join(counts)
    
class geneMGF_struct(geneMGF.geneMGF):
    def __init__(self, lineages, **kwargs):
        geneMGF.geneMGF.__init__(self, lineages, **kwargs)
        self.mgf_type = 'structured'
        
    def create_mgf_expr(self, coal_rates, migration_rates, initial_config):
        """ Generates the mgf, only should be called by __init__
        Arguments:
        coal_rates      -- a list giving the coalescent rate in each deme
                           ex: [1, 1]
        migration_rates -- a matrix giving the migration rates from each deme to each other deme
                           ex: sympy.Matrix([[0, 1], [1, 0]])
        initial_config  -- a list giving the starting configuration of the lineages
                           ex: [ [[1], [2]], [[3]] ]
        """
        
        num_demes = len(coal_rates)
        
        ## SET UP HELPER FUNCTIONS TO SOLVE THE SYSTEM RECURSIVELY
        def update_sol_dict(sol_dict, partition):
            # To update the solution dictionary one has to solve the system of equations for that level
            num_lineages = len(extract_lineages(partition))
            dummy_lineages = [[str(ii)] for ii in range(num_lineages)]
            system_partitions, system_solution = solve_markov(dummy_lineages, sol_dict)
            for ii in range(len(system_partitions)):
                sol_key = make_sol_key(geneMGF.deme_part_to_symbol(system_partitions[ii]))
                if sol_key not in sol_dict.keys():
                    sol_dict[sol_key] = [system_partitions[ii], 
                                         system_solution[geneMGF.deme_part_to_symbol(system_partitions[ii])]]
        
        def solve_markov(lineages, sol_dict):
            # Set up the system of equations
            Eqns = []
            deme_partitions = list(geneMGF.deme_partition(lineages, num_demes))
            deme_part_symbols = []
            for deme_partition in deme_partitions:
                deme_part_symbol = geneMGF.deme_part_to_symbol(deme_partition)
                deme_part_symbols.append(deme_part_symbol)
                deme_part_term = make_partition_term(deme_partition, coal_rates, migration_rates)
                deme_part_eqn = -deme_part_term*deme_part_symbol
                # All of the migration terms (which don't require known solutions)
                for deme_idx, deme in enumerate(deme_partition):
                    for lin in deme:
                        for deme_to_idx in range(num_demes):
                            deme_part_eqn += (migration_rates[deme_idx, deme_to_idx]*
                                              geneMGF.deme_part_symbol_after_migr(deme_partition, lin,
                                                                      deme_idx, deme_to_idx))
                # All of the coalescent terms (which do require known solutions)
                for deme_idx, deme in enumerate(deme_partition):
                    for pair in itertools.combinations(deme, 2):
                        symbol_after_coal = geneMGF.deme_part_symbol_after_coal(deme_partition, pair)
                        partition_after_coal = deme_part_after_coal(deme_partition, pair)
                        sol_key = make_sol_key(symbol_after_coal)
                        # If we already have a solution for this configuration after coalescence
                        if sol_key in sol_dict.keys():
                            solution = sol_dict[sol_key]
                            deme_part_eqn += coal_rates[deme_idx]*solution[1].subs(dummy_swap_list(solution[0], 
                                                                                                   partition_after_coal))
                        # If we don't already have a solution for this configuration
                        else:
                            update_sol_dict(sol_dict, partition_after_coal)
                            solution = sol_dict[sol_key]
                            deme_part_eqn += coal_rates[deme_idx]*solution[1].subs(dummy_swap_list(solution[0], 
                                                                                                   partition_after_coal))
                Eqns.append(deme_part_eqn)
            # Solve the system of equations
            print("Trying to solve a system of " + str(len(Eqns)) + " equations")
            # for Eqn in Eqns:
                # print(Eqn)
            solution = sympy.solve(Eqns, deme_part_symbols, simplify=False)
            return deme_partitions, solution
        
        solution_set = dict()
        for ii in range(num_demes):
            partition = []
            for jj in range(num_demes):
                if jj == ii:
                    partition.append(['1'])
                else: 
                    partition.append([])
            solution_set[make_sol_key(geneMGF.deme_part_to_symbol(partition))] = [partition, sympy.Integer(1)]
        # print(solution_set)
        
        partition, solution = solve_markov(self.lineages, solution_set)
        return solution#[geneMGF.deme_part_to_symbol(initial_config)]