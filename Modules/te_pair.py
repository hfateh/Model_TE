import numpy as np
from scipy.optimize import fsolve

import leg
reload(leg)

class TE_Pair(object):
    """ TE_pair
    """

    def __init__(self):

        """ """
        self.R_load_total = 1.0
        self.pairs = 128.0
        self.length = 1.e-3
        self.leg_area_ratio = 1.
        self.fill_fraction = .3
        self.Vs = 1.64/256. 
        self.R_internal = 1./256.
        self.Vs = 1.64/256.
        
        self.Ptype = leg.Leg()
        self.Ntype = leg.Leg()

        self.Ptype.material = 'HMS'
        self.Ntype.material = 'MgSi'

        self.nodes = 50
        self.t_array = np.linspace(0., 5., 20)
        self.T_h_conv = 500.
        self.T_c_conv = 300.
        self.U_hot = 500.
        self.U_cold = 500.

    def set_R_load(self):
        """ """
        self.R_load = (
            self.R_load_total / self.pairs
            )

    def set_I(self):

        """ """
        self.I = (
            self.Vs/(self.R_load + self.R_internal)
            )
        self.Ntype.I = -self.I
        self.Ptype.I = self.I

    def set_leg_areas(self):
        """ """
        self.Ntype.area = self.Ptype.area * self.leg_area_ratio
        self.area_void = (
            ((1.- self.fill_fraction) * (self.Ntype.area +
        self.Ptype.area)) / (self.fill_fraction)
            )
        self.area = (
            self.Ntype.area + self.Ptype.area + self.area_void
            )

    def set_constants(self):

        """ """
        self.set_R_load()
        self.set_leg_areas()
        self.set_I()
        self.Ntype.length = self.length
        self.Ptype.length = self.length
        self.Ptype.nodes = self.nodes
        self.Ntype.nodes = self.nodes
        self.Ntype.set_constants()
        self.Ptype.set_constants()
        self.Ntype.Vs = -self.Vs
        self.Ptype.Vs = self.Vs
        self.Ntype.R_internal = self.R_internal
        self.Ptype.R_internal = self.R_internal
        self.Ptype.t_array = self.t_array
        self.Ntype.t_array = self.t_array



        # self.Ptype.U_hot = self.U_hot
        # self.Ntype.U_hot = self.U_hot
        # self.Ptype.U_cold = self.U_cold
        # self.Ntype.U_cold = self.U_cold
        # self.Ptype.T_h_conv = self.T_h_conv
        # self.Ntype.T_h_conv = self.T_h_conv
        # self.Ptype.T_c_conv = self.T_c_conv
        # self.Ntype.T_c_conv = self.T_c_conv

    def set_q_guess(self):

        """ """
        self.Ntype.set_q_guess()
        self.Ptype.set_q_guess()

    def solve_te_pair_once(self):

        """ """
        self.Ntype.solve_leg_once(self.Ntype.q_h)
        self.Ptype.solve_leg_once(self.Ptype.q_h)
        self.T_c = self.Ntype.T_c

        self.q_h = (
            (self.Ptype.q_h * self.Ptype.area + self.Ntype.q_h *
             self.Ntype.area) / self.area * 0.001
            )

        self.q_c = (
            (self.Ptype.q_c * self.Ptype.area + self.Ntype.q_c *
             self.Ntype.area) / self.area * 0.001
            )

        self.Vs = -self.Ntype.Vs + self.Ptype.Vs

        self.R_internal = ( 
            self.Ntype.R_internal + self.Ptype.R_internal
            )

    def get_error(self, knob_arr):

        """ """
        self.Ntype.q_h = knob_arr[0]
        self.Ptype.q_h = knob_arr[1]
        self.Ntype.T_h = knob_arr[2]
        self.Ptype.T_h = knob_arr[3]
        self.I = knob_arr[4]

        # Following two lines are very important
        self.Ptype.I = self.I
        self.Ntype.I = -self.I

        #self.Ptype.T_h = self.T_h
        #self.Ntype.T_h = self.T_h

        self.solve_te_pair_once()

        self.Ntype.q_h_conv = (
            self.U_hot * (self.T_h_conv - self.Ntype.T_h)
            )
        self.Ptype.q_h_conv = (
            self.U_hot * (self.T_h_conv - self.Ptype.T_h)
            )
        self.Ntype.q_c_conv = (
            self.U_hot * (self.Ntype.T_c - self.T_c_conv)
            )
        self.Ptype.q_c_conv = (
            self.U_cold * (self.Ptype.T_c - self.T_c_conv)
            )

        #self.q_c_conv = self.U_cold * (self.T_c - self.T_c_conv)
        #self.q_h_conv = self.U_hot * (self.T_h_conv - self.T_h)

        #T_c_error = self.Ntype.T_c - self.Ptype.T_c
        #T_c_error = self.Ntype.T_c - self.Ptype.T_c

        Nq_c_error = self.Ntype.q_c - self.Ntype.q_c_conv
        Pq_c_error = self.Ptype.q_c - self.Ptype.q_c_conv
        Nq_h_error = self.Ntype.q_h_conv - self.Ntype.q_h
        Pq_h_error = self.Ptype.q_h_conv - self.Ptype.q_h

        self.I_correct = (
            self.Vs / (self.R_load + self.R_internal)
            )

        #q_h_error = self.q_h - self.q_h_conv
        I_error = self.I_correct - self.I

        self.error = (
            np.array([Nq_c_error, Pq_c_error, Nq_h_error, Pq_h_error, I_error]).flatten()
            )

        return self.error

    def solve_te_pair(self):

        """ """
        self.set_constants()
        self.Ptype.T_h = self.T_h_conv 
        self.Ntype.T_h = self.T_h_conv
        self.Ptype.T_c = self.T_c_conv
        self.Ntype.T_c = self.T_c_conv
        self.set_q_guess()

        knob_arr0 = (
            np.array([self.Ntype.q_h_guess, self.Ptype.q_h_guess,
            self.Ntype.T_h_conv, self.Ptype.T_h_conv, self.I])
            )

        self.fsolve_output = fsolve(self.get_error, x0=knob_arr0)
        self.V = self.I * self.R_load
        self.P = self.I * self.V
        self.P_total = self.P * self.pairs

    def solve_te_pair_transient_once(self):
        """ """
        self.Ntype.solve_leg_transient_once()
        self.Ptype.solve_leg_transient_once()
        
        self.Ntype.T_h_xt = self.Ntype.Txt[:,0]
        self.Ptype.T_h_xt = self.Ptype.Txt[:,0]
        self.Ntype.T_c_xt = self.Ntype.Txt[:,-1]
        self.Ptype.T_c_xt = self.Ptype.Txt[:,-1]
        
        self.Vs_transient = (
            - self.Ntype.Vs_transient + self.Ptype.Vs_transient
            )
        
        self.R_internal_transient = (
            self.Ntype.R_internal_transient +
            self.Ptype.R_internal_transient
            )

        self.I_transient = (
            self.Vs_transient / (self.R_load +
            self.R_internal_transient)
            )
        
        self.Power_transient = (
            self.I_transient * self.R_load
            )
        
        # Setting BCs for a new run of te_pair
        # for a new run of transient solution
        # self.Ntype.T_x = self.Ntype.Txt[-1,:]
        # self.Ptype.T_x = self.Ptype.Txt[-1,:]
        # self.Ntype.q_x = self.Ntype.qxt[-1,:]
        # self.Ptype.q_x = self.Ptype.qxt[-1,:]
        # self.Ntype.Vs_x = self.Ntype.Vsxt[-1,:]
        # self.Ptype.Vs_x = self.Ptype.Vsxt[-1,:]
        # self.Ntype.R_x = self.Ntype.Rxt[-1,:]
        # self.Ptype.R_x = self.Ptype.Rxt[-1,:]


    # def get_error_transient(self, error_transient):
    #     """ """
        
    #     self.Ntype.I_transient = error_transient[:,0]
    #     self.Ptype.I_transient = error_transient[:.1]
    #     self.Ntype.T_h_xt = error_transient[:,2]
    #     self.Ptype.T_h_xt = error_transient[:,3]

    #     self.solve_te_pair_transient_once()

    #     self.Ntype.T_c_xt = self.Ntype.Txt[:,-1]
    #     self.Ptype.T_c_xt = self.Ptype.Txt[:,-1]

    #     self.Ntype.T_h_xt = self.Ntype.Txt[:,0]
    #     self.Ptype.T_h_xt = self.Ptype.Txt[:,0]

    #     T_c_xt_error = (
    #         self.Ntype.T_c_xt - self.Ptype.T_c_xt
    #         )

    #     T_h_xt_error = (
    #         self.Ntype.T_h_xt - self.Ptype.T_h_xt
    #         )
        
    #     I_error_transient = (
    #         self.Ntype.I_transient - self.Ptype.I_transient
    #         )
        
    #     self.transient_te_error = (
    #         np.array([T_c_xt_error,T_h_xt_error,I_error_transient]).flatten()
    #         )

    #     return self.transient_te_error

    # def solve_te_pair_transient(self):
    #     """ """
        
    #     self.guess_transient = (np.array([]))
        
    #     self.fsolve_transient_output = (
    #         fsolve(self.get_error_transient, x0=self.guess_transient)
    #         )






# Bunch of errors which I need to organize and find a way to equate
# with respect to correct values














































































#=================================================
# This is a perfectly working code
# Separated on 14/02/2013
# After separating, I tried to solve te_pair my way, separate n and p
# type heat transfer rates
#=================================================

# import numpy as np
# from scipy.optimize import fsolve

# import leg
# reload(leg)

# class TE_Pair(object):
#     """ TE_pair
#     """

#     def __init__(self):

#         """ """
#         self.R_load_total = 1.0
#         self.pairs = 128.0
#         self.length = 1.e-3
#         self.leg_area_ratio = 1.
#         self.fill_fraction = .3
#         self.Vs = 1.64/256. 
#         self.R_internal = 1./256.
#         self.Vs = 1.64/256.
        
#         self.Ptype = leg.Leg()
#         self.Ntype = leg.Leg()

#         self.Ptype.material = 'HMS'
#         self.Ntype.material = 'MgSi'

#         self.nodes = 50
#         self.t_array = np.linspace(0., 5., 20)
#         self.T_h_conv = 500.
#         self.T_c_conv = 300.
#         self.U_hot = 500.
#         self.U_cold = 500.

#     def set_R_load(self):
#         """ """
#         self.R_load = (
#             self.R_load_total / self.pairs
#             )

#     def set_I(self):

#         """ """
#         self.I = (
#             self.Vs/(self.R_load + self.R_internal)
#             )
#         self.Ntype.I = -self.I
#         self.Ptype.I = self.I

#     def set_leg_areas(self):
#         """ """
#         self.Ntype.area = self.Ptype.area * self.leg_area_ratio
#         self.area_void = (
#             ((1.- self.fill_fraction) * (self.Ntype.area +
#         self.Ptype.area)) / (self.fill_fraction)
#             )
#         self.area = (
#             self.Ntype.area + self.Ptype.area + self.area_void
#             )

#     def set_constants(self):

#         """ """
#         self.set_R_load()
#         self.set_leg_areas()
#         self.set_I()
#         self.Ntype.length = self.length
#         self.Ptype.length = self.length
#         self.Ptype.nodes = self.nodes
#         self.Ntype.nodes = self.nodes
#         self.Ntype.set_constants()
#         self.Ptype.set_constants()
#         self.Ntype.Vs = -self.Vs
#         self.Ptype.Vs = self.Vs
#         self.Ntype.R_internal = self.R_internal
#         self.Ptype.R_internal = self.R_internal
#         self.Ptype.t_array = self.t_array
#         self.Ntype.t_array = self.t_array



#         # self.Ptype.U_hot = self.U_hot
#         # self.Ntype.U_hot = self.U_hot
#         # self.Ptype.U_cold = self.U_cold
#         # self.Ntype.U_cold = self.U_cold
#         # self.Ptype.T_h_conv = self.T_h_conv
#         # self.Ntype.T_h_conv = self.T_h_conv
#         # self.Ptype.T_c_conv = self.T_c_conv
#         # self.Ntype.T_c_conv = self.T_c_conv

#     def set_q_guess(self):

#         """ """
#         self.Ntype.set_q_guess()
#         self.Ptype.set_q_guess()

#     def solve_te_pair_once(self):

#         """ """
#         self.Ntype.solve_leg_once(self.Ntype.q_h)
#         self.Ptype.solve_leg_once(self.Ptype.q_h)
#         self.T_c = self.Ntype.T_c

#         self.q_h = (
#             (self.Ptype.q_h * self.Ptype.area + self.Ntype.q_h *
#              self.Ntype.area) / self.area * 0.001
#             )

#         self.q_c = (
#             (self.Ptype.q_c * self.Ptype.area + self.Ntype.q_c *
#              self.Ntype.area) / self.area * 0.001
#             )

#         self.Vs = -self.Ntype.Vs + self.Ptype.Vs

#         self.R_internal = ( 
#             self.Ntype.R_internal + self.Ptype.R_internal
#             )

#     def get_error(self, knob_arr):

#         """ """
#         self.Ntype.q_h = knob_arr[0]
#         self.Ptype.q_h = knob_arr[1]
#         self.T_h = knob_arr[2]
#         self.I = knob_arr[3]

#         # Following two lines are very important
#         self.Ptype.I = self.I
#         self.Ntype.I = -self.I

#         self.Ptype.T_h = self.T_h
#         self.Ntype.T_h = self.T_h

#         self.solve_te_pair_once()

#         self.q_c_conv = self.U_cold * (self.T_c - self.T_c_conv)
#         self.q_h_conv = self.U_hot * (self.T_h_conv - self.T_h)
        
#         self.I_correct = (
#             self.Vs / (self.R_load + self.R_internal)
#             )

#         T_c_error = self.Ntype.T_c - self.Ptype.T_c
#         # T_c_error = 0
#         q_c_error = self.q_c - self.q_c_conv
#         #q_c_error = 0
#         q_h_error = self.q_h - self.q_h_conv
#         #q_h_error = 0
#         I_error = self.I_correct - self.I

#         self.error = (
#             np.array([T_c_error, q_c_error, q_h_error, I_error]).flatten()
#             )
#         return self.error

#     def solve_te_pair(self):

#         """ """
#         self.set_constants()
#         self.Ptype.T_h = self.T_h_conv 
#         self.Ntype.T_h = self.T_h_conv
#         self.Ptype.T_c = self.T_c_conv
#         self.Ntype.T_c = self.T_c_conv
#         self.set_q_guess()

#         knob_arr0 = (
#             np.array([self.Ntype.q_h_guess, self.Ptype.q_h_guess,
#         self.T_h_conv, self.I])
#             )

#         self.fsolve_output = fsolve(self.get_error, x0=knob_arr0)
#         self.V = self.I * self.R_load
#         self.P = self.I * self.V
#         self.P_total = self.P * self.pairs

#     def solve_te_pair_transient_once(self):
#         """ """
#         self.Ntype.solve_leg_transient_once()
#         self.Ptype.solve_leg_transient_once()

#         # CALCULATE EACH VARIABLE FOR THE PAIR
        
#         self.Ntype.T_h_xt = self.Ntype.Txt[:,0]
#         self.Ptype.T_h_xt = self.Ptype.Txt[:,0]
#         self.Ntype.T_c_xt = self.Ntype.Txt[:,-1]
#         self.Ptype.T_c_xt = self.Ptype.Txt[:,-1]
        
#         self.Vs_transient = (
#             - self.Ntype.Vs_transient + self.Ptype.Vs_transient
#             )
        
#         self.R_internal_transient = (
#             self.Ntype.R_internal_transient +
#             self.Ptype.R_internal_transient
#             )

#         self.I_transient = (
#             self.Vs_transient / (self.R_load +
#             self.R_internal_transient)
#             )
        
#         self.Power_transient = (
#             self.I_transient * self.R_load
#             )
        
#         # Setting BCs for a new run of te_pair
#         # for a new run of transient solution
#         # self.Ntype.T_x = self.Ntype.Txt[-1,:]
#         # self.Ptype.T_x = self.Ptype.Txt[-1,:]
#         # self.Ntype.q_x = self.Ntype.qxt[-1,:]
#         # self.Ptype.q_x = self.Ptype.qxt[-1,:]
#         # self.Ntype.Vs_x = self.Ntype.Vsxt[-1,:]
#         # self.Ptype.Vs_x = self.Ptype.Vsxt[-1,:]
#         # self.Ntype.R_x = self.Ntype.Rxt[-1,:]
#         # self.Ptype.R_x = self.Ptype.Rxt[-1,:]







































































# there is a working te_pair at the bottom
# scroll down
































































# import numpy as np
# from scipy.optimize import fsolve

# import leg
# reload(leg)

# class TE_Pair(object):
#     """ TE_pair
#     """

#     def __init__(self):

#         """ """
#         self.R_load_total = 1.0
#         self.pairs = 128.0
#         # self.R_load = 1.0/128.0
#         self.length = 1.e-3
#         self.leg_area_ratio = 0.8
#         self.fill_fraction = 0.3
#         self.Vs = 1.64/256. 
#         self.R_internal = 1./256.
#         self.Vs = 1.64/256.
        
#         self.Ptype = leg.Leg()
#         self.Ntype = leg.Leg()

#         self.Ptype.material = 'HMS'
#         self.Ntype.material = 'MgSi'

#         self.nodes = 10
#         # self.set_constants() was initially here

#     def set_R_load(self):
#         """ """
#         self.R_load = (
#             self.R_load_total / self.pairs
#             )

#     def set_I(self):

#         """ """
#         self.I = (
#             self.Vs/(self.R_load + self.R_internal)
#             )
#         # self.set_leg_areas() #was initially here
#         self.Ntype.I = -self.I
#         self.Ptype.I = self.I

#     def set_leg_areas(self):
#         """ """
#         self.Ntype.area = self.Ptype.area * self.leg_area_ratio
#         self.area_void = (
#             ((1.- self.fill_fraction) * (self.Ntype.area +
#         self.Ptype.area)) / (self.fill_fraction)
#             )
#         self.area = (
#             self.Ntype.area + self.Ptype.area + self.area_void
#             )

#     def set_constants(self):

#         """ """
#         self.set_R_load()
#         self.set_leg_areas()
#         self.set_I()
#         self.Ntype.length = self.length
#         self.Ptype.length = self.length
#         self.Ptype.nodes = self.nodes
#         self.Ntype.nodes = self.nodes
#         self.Ntype.set_constants()
#         self.Ptype.set_constants()
#         self.Ntype.Vs = -self.Vs
#         self.Ptype.Vs = self.Vs
#         self.Ntype.R_internal = self.R_internal
#         self.Ptype.R_internal = self.R_internal

#     def set_q_guess(self):

#         """ """
#         self.Ntype.set_q_guess()
#         self.Ptype.set_q_guess()

#     # I dont know why this is here, its unnecessary, delete this
#     # def set_TEproperties(self, T_props):

#     #     """ """
#     #     self.Ntype.set_TEproperties(T_props)
#     #     self.Ptype.set_TEproperties(T_props)

#     def solve_te_pair_once(self):

#         """ """
#         self.Ntype.solve_leg_once(self.Ntype.q_h)
#         self.Ptype.solve_leg_once(self.Ptype.q_h)
#         self.T_c = self.Ntype.T_c

#         self.q_h = (
#             (self.Ptype.q_h * self.Ptype.area + self.Ntype.q_h *
#              self.Ntype.area) / self.area * 0.001
#             )

#         self.q_c = (
#             (self.Ptype.q_c * self.Ptype.area + self.Ntype.q_c *
#              self.Ntype.area) / self.area * 0.001
#             )

#         self.Vs = -self.Ntype.Vs + self.Ptype.Vs

#         self.R_internal = ( 
#             self.Ntype.R_internal + self.Ptype.R_internal
#             )

#     def get_error(self, knob_arr):

#         """ """
#         self.Ntype.q_h = knob_arr[0]
#         self.Ptype.q_h = knob_arr[1]
#         self.T_h = knob_arr[2]
#         self.I = knob_arr[3]

#         self.Ptype.T_h = self.T_h
#         self.Ntype.T_h = self.T_h

#         self.solve_te_pair_once()

#         self.q_c_conv = self.U_cold * (self.T_c - self.T_c_conv)
#         self.q_h_conv = self.U_hot * (self.T_h_conv - self.T_h)
#         #print "I in this run is ", self.I
        
#         self.I_correct = (
#             self.Vs / (self.R_load + self.R_internal)
#             )

#         T_c_error = self.Ntype.T_c - self.Ptype.T_c
#         q_c_error = self.q_c - self.q_c_conv
#         q_h_error = self.q_h - self.q_h_conv
#         I_error = self.I_correct - self.I
#         # print "\n"
#         # print "\n"
#         # print "\nI_correct is ", self.I_correct
#         # print "I_guess is ", self.I
#         # print "I_error is ", I_error
#         # print "\n"
#         # print "\n"

#         self.error = (
#             np.array([T_c_error, q_c_error, q_h_error, I_error]).flatten()
#             )
#         return self.error

#     def solve_te_pair(self):

#         """ """
#         self.set_constants()
#         self.Ptype.T_h = self.T_h_conv 
#         self.Ntype.T_h = self.T_h_conv
#         self.Ptype.T_c = self.T_c_conv
#         self.Ntype.T_c = self.T_c_conv
#         # self.set_I()  # was initially here
#         self.set_q_guess()

#         knob_arr0 = (
#             np.array([self.Ntype.q_h_guess, self.Ptype.q_h_guess,
#         self.T_h_conv, self.I])
#             )

#         # self.Ptype.T_c_goal = None
#         # self.Ntype.T_c_goal = None

#         self.fsolve_output = fsolve(self.get_error, x0=knob_arr0)
#         self.V = self.I * self.R_load
#         self.P = self.I * self.V
#         self.P_total = self.P * self.pairs





