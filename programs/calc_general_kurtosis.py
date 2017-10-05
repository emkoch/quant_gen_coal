import pickle
import sympy
import mgf_math

# Generate the fourth order mgf approximation
test_mgf_4 = mgf_math.traitMGF(nIndiv=4)
test_mgf_4.make_mgf(order=4, approx_type="exchangeable")

# Save at this point
with open("test_mgf_4.pyc", "wb") as fout:
    pickle.dump(test_mgf_4, fout)

# Calculate moments necessary for kurtosis
mom_x1x1x1x1_4 = test_mgf_4.calc_moment(pows=[4,0,0,0],
                                        approx_type="exchangeable").simplify().expand()
with open("mom_x1x1x1x1_4.pyc", "wb") as fout:
    pickle.dump(mom_x1x1x1x1_4, fout)

mom_x1x1x1x2_4 = test_mgf_4.calc_moment(pows=[3,1,0,0],
                                        approx_type="exchangeable").simplify().expand()
with open("mom_x1x1x1x2_4.pyc", "wb") as fout:
    pickle.dump(mom_x1x1x1x2_4, fout)

mom_x1x1x2x3_4 = test_mgf_4.calc_moment(pows=[2,1,1,0],
                                        approx_type="exchangeable").simplify().expand()
with open("mom_x1x1x2x3_4.pyc", "wb") as fout:
    pickle.dump(mom_x1x1x2x3_4, fout)

mom_x1x2x3x4_4 = test_mgf_4.calc_moment(pows=[1,1,1,1],
                                        approx_type="exchangeable").simplify().expand()
with open("mom_x1x2x3x4_4.pyc", "wb") as fout:
    pickle.dump(mom_x1x2x3x4_4, fout)

# Save the overall mgf again
with open("test_mgf_4.pyc", "wb") as fout:
    pickle.dump(test_mgf_4, fout)
