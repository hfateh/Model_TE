import numpy as np
from scipy.optimize import fsolve
import time

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
        self.R_internal = 1./256.
        self.Vs = 1.64/256.
        
        self.Ptype = leg.Leg()
        self.Ntype = leg.Leg()

        # All other void area, Ntype area, total te_pair area are
        # calculated based on Ptype area
        self.Ptype.area = (3.0e-3) ** 2 

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
	leg_area_ratio = self.leg_area_ratio
	fill_fraction = self.fill_fraction
        self.Ntype.area = self.Ptype.area * leg_area_ratio
        self.area_void = (
            ((1.- fill_fraction) * (self.Ntype.area +
        self.Ptype.area)) / (fill_fraction)
            )
        self.area = (
            self.Ntype.area + self.Ptype.area + self.area_void
            )

    def set_t_array(self):
        """only used when transient solutions are stacked"""
        self.Ntype.t_array = self.t_array
        self.Ptype.t_array = self.t_array

    def set_constants(self):

        """ """
        self.set_leg_areas()
        self.set_R_load()
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

        self.Ptype.U_hot = self.U_hot
        self.Ntype.U_hot = self.U_hot
        self.Ptype.U_cold = self.U_cold
        self.Ntype.U_cold = self.U_cold
        self.Ptype.T_h_conv = self.T_h_conv
        self.Ntype.T_h_conv = self.T_h_conv
        self.Ptype.T_c_conv = self.T_c_conv
        self.Ntype.T_c_conv = self.T_c_conv

    def set_q_guess(self):

        """ """
        self.Ntype.set_q_guess()
        self.Ptype.set_q_guess()

    def solve_te_pair_once(self):

        """ """
        self.Ntype.solve_leg_once(self.Ntype.q_h)
        self.Ptype.solve_leg_once(self.Ptype.q_h)
        # self.T_c = self.Ntype.T_c

        # self.q_h = (
        #     (self.Ptype.q_h * self.Ptype.area + self.Ntype.q_h *
        #      self.Ntype.area) / self.area * 0.001
        #     )

        # self.q_c = (
        #     (self.Ptype.q_c * self.Ptype.area + self.Ntype.q_c *
        #      self.Ntype.area) / self.area * 0.001
        #     )

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

        Nq_c_error = self.Ntype.q_c - self.Ntype.q_c_conv
        Pq_c_error = self.Ptype.q_c - self.Ptype.q_c_conv
        Nq_h_error = self.Ntype.q_h_conv - self.Ntype.q_h
        Pq_h_error = self.Ptype.q_h_conv - self.Ptype.q_h

        self.I_correct = (
            self.Vs / (self.R_load + self.R_internal)
            )

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
        self.P_flux = self.P / self.area
        self.P_total = self.P * self.pairs
        # self.P_flux = (self.P_total / self.area)

    def solve_te_pair_transient_once(self):
        """ """
        self.Ntype.T_h_conv = self.T_h_conv
        self.Ptype.T_h_conv = self.T_h_conv
        self.Ntype.T_c_conv = self.T_c_conv
        self.Ptype.T_c_conv = self.T_c_conv

        self.Ntype.solve_leg_transient_once()
        self.Ptype.solve_leg_transient_once()
        
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

        self.Vs_transient_total = (
            self.Vs_transient * self.pairs
            )

        # this bottom code is just a check
        self.Power_transient_total1 = (
            self.I_transient * self.R_load_total
            )

        self.Power_transient_total = (
            self.Power_transient * self.pairs
            )

        self.q_c_transient = (
            (self.Ntype.qxt[:,-1] * self.Ntype.area +
            self.Ptype.qxt[:,-1] * self.Ptype.area) / self.area
            )

        self.q_h_transient = (
            (self.Ntype.qxt[:,0] * self.Ntype.area +
            self.Ptype.qxt[:,0] * self.Ptype.area) / self.area
            )
        
        # set BCs for another run of te_transient_inst
        self.Ntype.T_x = self.Ntype.Txt[-1,:]
        self.Ptype.T_x = self.Ptype.Txt[-1,:]
        self.Ntype.q_x = self.Ntype.qxt[-1,:]
        self.Ptype.q_x = self.Ptype.qxt[-1,:]
        self.Ntype.Vs_x = self.Ntype.Vsxt[-1,:]
        self.Ptype.Vs_x = self.Ptype.Vsxt[-1,:]
        self.Ntype.R_x = self.Ntype.Rxt[-1,:]
        self.Ptype.R_x = self.Ptype.Rxt[-1,:]


# The results of this transient model are transient
# Nodal temperature distribution
# Nodal Seebeck voltage distribution
# Nodal heat flux distribution
# Nodal internal resistance distribution
# Internal resistance
# Current
# Power
# Hot side and cold side heat flux distribution

# I need to  concatenate these for various runs
    

    # def solve_te_pair_transient(self):
    #     """ """

    #     try:
    #         self.NTxt_array

    #     except AttributeError:
    #         self.NT_xt_array = self.Ntype.T_xt

    #     else:
    #         self.NT_xt_array = np.concatenate((self.Ntype.T_xt,            

    #     try: 
    #         self.T_xt

    #     except AttributeError:
    #         self.odeint_output = odeint(
    #             self.get_dTx_dt, y0=self.y0, t=self.t_array,
    #             full_output=1 
    #             )
    #         self.T_xt = self.odeint_output[0]

    #     # doesnt really get here
    #     # else:
    #     #     self.y0 = self.T_xt[-1,:]
    #     #     self.odeint_output = odeint(
    #     #         self.get_dTx_dt, y0=self.y0, t=self.t_array,
    #     #         full_output=1 
    #     #         )
    #     #     self.T_xt = np.concatenate((self.T_xt, self.odeint_output[0]))




















#====================================================
#====================================================
#Optimization works
#====================================================
#OPTIMIZATION
#====================================================
#calling get_min_par should be able to vary the variables that wea
#are trying to optimize.

    def optimize(self):
        """ """
        self.optimize1()
        # self.optimize2()

    def set_opt_constants(self):
        """ """
        self.set_leg_areas()
        self.set_R_load()

    def get_minpar1(self, apar):

        """Returns inverse of power flux.

        Methods:

        self.set_leg_areas

        Used by method self.optimize

        self.length = apar[0]
        self.fill_fraction = apar[1]
        self.I = apar[2]
        self.leg_area_ratio = apar[3]

        Use with scipy.optimize.fmin to find optimal set of input
        parameters.

        This method uses power flux rather than power because for
        optimal power, leg height approaches zero and void area
        approaches infinity.  This trivial result is not useful."""

        self.opt_iter = self.opt_iter + 1
        if self.opt_iter % 15 == 0:
            print "\noptimizaton iteration ", self.opt_iter
            print "fill_fraction ", self.fill_fraction
            print "R_load_total ", self.R_load_total
            print "leg_area_ratio ", self.leg_area_ratio
            print "length ", self.length
        #     print "leg length =", self.length, "m"
        #     print "fill fraction =", self.fill_fraction * 100., "%"
        #     print "current =", self.I, "A"
        #     print "area ratio =", self.leg_area_ratio
        #     print "power flux (kW/m^2)", self.P_flux
        # # apar = np.array(apar)

        self.R_load_total = apar[0]
        self.leg_area_ratio = apar[1]
        self.length = apar[2]
        #self.fill_fraction = apar[3]
        #self.R_load_total = apar[0]
        #self.fill_fraction = apar[1]
        #self.I = apar[2]
        #self.leg_area_ratio = apar[1]
        #self.fill_fraction = apar[0]

        # reset surrogate variables
        # self.set_constants()
        self.set_opt_constants()

        self.solve_te_pair()

        if (apar <= 0.).any():
            minpar = np.abs(self.P_flux) ** 3. + 100
            print "Encountered impossible value."

        else:
            minpar = - self.P_flux

        return minpar

    def optimize1(self):

        """ Optimizes R_load_total, leg_area_ratio, and length """

        time.clock()

        # dummy function that might be used with minimization
        def fprime():
            return 1

        self.opt_iter = 0

        # self.x0 = (
        #     np.array([self.leg_area_ratio,
        #     self.length, self.fill_fraction])
        #     )

        self.x0 = (
            np.array([self.R_load_total, self.leg_area_ratio, self.length])
            )

        from scipy.optimize import fmin

        self.xmin = fmin(self.get_minpar1, self.x0)

        t1 = time.clock()

        print '\n'

        print "Optimized parameters:"
        # print "leg length =", self.length, "m"
        # print "fill fraction =", self.fill_fraction * 100., "%"
        # print "current =", self.I, "A"
        # print "area ratio =", self.leg_area_ratio

        # print "\npower:", self.P * 1000., 'W'
        # print "power flux:", self.P_total, "kW/m^2"
        # print "Optimum length is", self.length

        print "R_load_total ", self.R_load_total        
        print "leg_area_ratio ", self.leg_area_ratio
        print "Length ", self.length
        # print "Optimum fill_fraction is ", self.fill_fraction
        print """Elapsed time solving xmin1 =""", t1

     

    # def get_minpar2(self, apar):

    #     """Returns inverse of power flux.

    #     Methods:

    #     self.set_leg_areas

    #     Used by method self.optimize

    #     self.length = apar[0]
    #     self.fill_fraction = apar[1]
    #     self.I = apar[2]
    #     self.leg_area_ratio = apar[3]

    #     Use with scipy.optimize.fmin to find optimal set of input
    #     parameters.

    #     This method uses power flux rather than power because for
    #     optimal power, leg height approaches zero and void area
    #     approaches infinity.  This trivial result is not useful."""

    #     self.opt_iter = self.opt_iter + 1
    #     if self.opt_iter % 15 == 0:
    #         print "\noptimizaton iteration ", self.opt_iter
    #         print "fill_fraction ", self.fill_fraction

    #     self.fill_fraction = apar[0]

    #     # reset surrogate variables
    #     # self.set_constants()
    #     self.set_opt_constants()

    #     self.solve_te_pair()

    #     if (apar >= 1.).any():
    #         minpar = np.abs(self.P_flux) ** 3. + 100
    #         print "Encountered impossible value."

    #     else:
    #         minpar = - self.P_flux

    #     return minpar

    # def optimize2(self):

    #     """ Optimizes R_load_total, leg_area_ratio, and length """

    #     time.clock()

    #     # dummy function that might be used with minimization
    #     def fprime():
    #         return 1

    #     self.opt_iter = 0

    #     self.x0 = (
    #         np.array([self.fill_fraction])
    #         )

    #     from scipy.optimize import fmin

    #     self.xmin = fmin(self.get_minpar2, self.x0)

    #     t1 = time.clock()

    #     print '\n'

    #     print "Optimized parameters:"
        
    #     print "R_load_total ", self.R_load_total        
    #     print "leg_area_ratio ", self.leg_area_ratio
    #     print "Length ", self.length
    #     print "fill_fraction ", self.fill_fraction
    #     print """Elapsed time solving xmin1 =""", t1



















































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





