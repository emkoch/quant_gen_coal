import pickle
import sympy
import mgf_math
import geneMGF
import geneMGF_faststruct

eta = sympy.symbols("eta")
coal_rates = [eta, eta]
lineages = [["0"], ["1"], ["2"], ["3"]]
initial_config = [[["0"], ["1"], ["2"], ["3"]], []]
migration_rates = sympy.Matrix([[0, "m"], ["m", 0]])
gene_mgf_deme_4 = geneMGF_faststruct.geneMGF_faststruct(lineages=lineages,
                                                        coal_rates=coal_rates,
                                                        migration_rates=migration_rates,
                                                        initial_config=initial_config)

with open("gene_mgf_deme_4.pyc", "wb") as fout:
    pickle.dump(gene_mgf_deme_4, fout)
