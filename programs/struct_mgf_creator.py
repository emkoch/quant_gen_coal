import geneMGF
import geneMGF_struct
import pickle
import mgf_math
import sympy
import time

eta = sympy.symbols("eta")
coal_rates = [eta, eta]
lineages = [[1], [2], [3], [4]]
initial_config = [[[1], [2], [3], [4]], []]
migration_rates = sympy.Matrix([[0, "m"], ["m", 0]])

start = time.time()
mgf_four_lin = geneMGF_struct.geneMGF_struct(lineages = lineages, coal_rates = coal_rates, 
                                         migration_rates = migration_rates, initial_config = initial_config)
print(time.time() - start)

with open("mgf_four_lin.pyc", "wb") as fout:
    pickle.dump(mgf_four_lin, fout)