import numpy as np
from scipy.optimize import fsolve
import leg
reload(leg)

class TE_Pair(object):
    """ TE_pair
    """

    def __init__(self):

        """ """
        self.R_load = 1.0/256.0
        self.length = 1.e-3
        self.leg_area_ratio = 0.8
        self.fill_fraction = 0.3
        self.Vs = 1.64/256.         # initial guess for Voc
        self.R_internal = 1./256.   # initial guess for R_internal

        self.Ptype = leg.Leg()
        self.Ntype = leg.Leg()

        self.Ptype.material = 'HMS'
        self.Ntype.material = 'MgSi'

        self.nodes = 10
        self.set_J()
        self.set_constants()

    def set_J(self):

        """ """
        self.J = self.Vs / (self.R_load + self.R_internal)

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
        self.set_leg_areas()
        self.Ntype.length = self.length
        self.Ptype.length = self.length
        self.Ptype.nodes = self.nodes
        self.Ntype.nodes = self.nodes
        self.Ntype.set_constants()
        self.Ptype.set_constants()
        self.Ntype.R_internal = self.R_internal
        self.Ptype.R_internal = self.R_internal
        self.Ntype.J = - self.J
        self.Ptype.J = self.J

    def set_q_guess(self):

        """ """
        self.Ntype.set_q_guess()
        self.Ptype.set_q_guess()

    # I dont know why this is here, its unnecessary, delete this
    # def set_TEproperties(self, T_props):

    #     """ """
    #     self.Ntype.set_TEproperties(T_props)
    #     self.Ptype.set_TEproperties(T_props)

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

    def get_error(self, knob_arr):

        """ """
        self.Ntype.q_h = knob_arr[0]
        self.Ptype.q_h = knob_arr[1]
        self.T_h = knob_arr[2]
        self.J = knob_arr[3]

        self.Ptype.T_h = self.T_h
        self.Ntype.T_h = self.T_h

        self.solve_te_pair_once()

        self.q_c_conv = self.U_cold * (self.T_c - self.T_c_conv)
        self.q_h_conv = self.U_hot * (self.T_h_conv - self.T_h)

        self.J_correct = (
            self.Vs / (self.R_load + self.R_internal)
            )

        T_c_error = self.Ntype.T_c - self.Ptype.T_c
        q_c_error = self.q_c - self.q_c_conv
        q_h_error = self.q_h - self.q_h_conv
        J_error = self.J_correct - self.J

        # print "Error in T_c is", T_c_error
        # print "Error in q_c is", q_c_error
        # print "Error in q_h is", q_h_error
        # print "Error in J is", J_error

        self.error = (
            np.array([T_c_error, q_c_error, q_h_error, J_error]).flatten()
            )
        return self.error

    def solve_te_pair(self):

        """ """

        self.Ptype.T_h = self.T_h_conv 
        self.Ntype.T_h = self.T_h_conv
        self.Ptype.T_c = self.T_c_conv
        self.Ntype.T_c = self.T_c_conv
        self.set_q_guess()

        knob_arr0 = (
            np.array([self.Ntype.q_h_guess, self.Ptype.q_h_guess,
        self.T_h_conv, self.J])
            )

        self.Ptype.T_c_goal = None
        self.Ntype.T_c_goal = None

        self.fsolve_output = fsolve(self.get_error, x0=knob_arr0)

        self.P = (self.Ntype.P + self.Ptype.P) * 0.001
        self.P_flux = self.P / self.area
        self.Vs = -self.Ntype.Vs + self.Ptype.Vs
        self.V = self.J * self.R_load / self.area
        self.R_internal = ( 
            self.Ntype.R_internal + self.Ptype.R_internal
            )


    # Now we don't need this following technique since I have already
    # included the J_error inside fsolve


    # def get_J_error(self, J):
    #     """Return the error in actual and guessed J value
    #     """
    #     # print "New Guess for J is", self.J
    #     # print "Internal resistance is ", self.R_internal
        
    #     self.solve_te_pair()
        
    #     self.J_correct = (
    #         self.Vs / (self.R_load + self.R_internal)
    #         )
        
    #     self.J_error = self.J_correct - self.J
    #     print "Seebeck voltage is ", self.Vs
    #     print "The error in J is ", self.J_error
    #     return self.J_error

    # def solve_te_pair_for_real(self):
    #     """ """
    #     self.fsolve_output0 = fsolve(self.get_J_error, x0= self.J)

