# distribution modules
import os
import sys

# local user modules
cmd_folder = os.path.dirname(os.path.abspath('../Modules/hx.py'))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

import te_pair
reload(te_pair)

te_pair = te_pair.TE_Pair()

te_pair.Ntype.material = 'MgSi'
te_pair.Ptype.material = 'HMS'

te_pair.T_c_conv = 300.
te_pair.T_h_conv = 800.

te_pair.U_cold = 800000.
te_pair.U_hot = 20000.

te_pair.solve_te_pair_for_real()

# print "\nLoad resistance is ", te_pair.R_load
# print "\nNtype T distribution is \n", te_pair.Ntype.T_x
# print "\nPtype T distribution is \n", te_pair.Ptype.T_x

# print "\nte_pair q_h is ", te_pair.q_h
# print "\nte_pair q_c is ", te_pair.q_c
# print "\nPower output for te_pair is ", te_pair.P
