import types
import numpy as np
from scipy.integrate import odeint
from numpy.testing import assert_approx_equal
from scipy.optimize import fsolve

import mat_prop
reload (mat_prop)


class Leg(object):

    """ Class for individual TE leg.

    Methods:
    


    """

    def __init__(self):

        """ sets constants and binds methods

        """
        # Initial guess
        self.I_guess = 1.5
        self.Vs = 1.64/256.0
        self.R_internal = 1.0/256

        self.material = 'HMS'
        self.T_h_conv = 750.
        self.T_c_conv = 360.
        self.U_hot = 54.e2
        self.U_cold = 253.e3
        self.R_load = 1.0
        self.nodes = 10
        self.length = 1.5e-3
        self.area = (3.0e-3) ** 2 
        self.set_constants()
        
        self.import_raw_property_data = ( 
            types.MethodType(mat_prop.import_raw_property_data, self)
            )
        self.set_properties_v_temp = (
            types.MethodType(mat_prop.set_properties_v_temp, self)
            )
        self.set_TEproperties = (
            types.MethodType(mat_prop.set_TEproperties, self)
            )

    def set_constants(self):

        self.x = np.linspace(0., self.length, self.nodes)

    def set_ZT(self):

        self.ZT = (
            self.alpha ** 2. * self.T_props / (self.k * self. rho)
            )

    def set_J(self):

        self.J = self.Vs / (self.R_load + self.R_internal)

    def solve_leg(self):
        
        self.T_h = self.T_h_conv
        self.T_c = self.T_c_conv

        self.fsolve_output = fsolve(self.get_error, x0=self.T_h-1.)

    def get_error(self, T_h):

        self.T_h = T_h[0]
        self.q_h_conv = (
            self.U_hot * (self.T_h_conv - self.T_h)
            )
        self.q_h = self.q_h_conv
        self.solve_leg_once(self.q_h)
        self.q_c_conv = self.U_cold * (self.T_c - self.T_c_conv)
        self.q_c_error = self.q_c - self.q_c_conv
        
        return self.q_c_error
    
    def solve_leg_once(self, q_h):

        self.q_h = q_h
        self.y0 = np.array([self.T_h, self.q_h, 0, 0])
        self.y = odeint(self.get_dTq_dx, y0=self.y0, t=self.x)
        
        self.T_x = self.y[:,0]
        self.q_x = self.y[:,1]
        self.Vs_x = self.y[:,2]
        self.R_int_x = self.y[:,3]
        
        self.T_c = self.T_x[-1]
        self.q_c = self.q_x[-1]
        
        self.Vs = self.Vs_x[0] - self.Vs_x[-1]
        self.V = self.R_load * (self.J * self.area)
        self.R_internal = self.R_int_x[-1]

        # Multiply q_h and q_c by area to get rid of this error

        # self.P = self.R_load * (self.J / self.area) ** 2
        self.P = self.V * (self.J * self.area)

        # # Sanity check.  q_h - q_c should be nearly equal but not
        # # exactly equal to P.  It is not exact because of spatial
        # # asymmetry in electrical resistivity along the leg.  I
        # # imported assert_approx_equal in the front matter to make
        # # this print an error if there is too much disagreement.

        # self.P_from_heat = (self.q_h - self.q_c) * self.area

        # sig_figs = 3
        # try:
        #     assert_approx_equal(self.P, self.P_from_heat, sig_figs)
        # except AssertionError:
        #     print "\nPower from q_h - q_c and I**2 * R disagree."
        #     print "Consider reducing sig_figs under solve_leg_once"
        #     print "in leg.py if you think this is an error."

    def get_dTq_dx(self, Tq, x):

        #print "Tq array inside dTq_dx is \n", Tq
        T = Tq[0]
        q = Tq[1]
        
        self.set_TEproperties(T)
        self.set_ZT()
        self.I = self.J * self.area

        dT_dx = (
            (self.J * self.alpha * T - q) / self.k
            )
        dq_dx = (
            (self.rho * self.J ** 2. * (1+self.ZT)) - (self.J *
            self.alpha * q / self.k)
            )
        dVs_dx = (self.alpha * dT_dx)
        dR_dx = (self.rho / self.area)

        return dT_dx, dq_dx, dVs_dx, dR_dx
    
    def set_q_guess(self):

        """Sets guess for q_c to be used by iterative solutions.
        """

        self.T_props = 0.5 * (self.T_h + self.T_c)
        self.set_TEproperties(T_props=self.T_props)
        delta_T = self.T_h - self.T_c

        self.q_c = - (
            self.alpha * self.T_c * self.J - delta_T / self.length *
        self.k - self.J ** 2 * self.length * self.rho
            )

        self.q_h = - (
            self.alpha * self.T_h * self.J - delta_T / self.length *
            self.k + self.J ** 2. * self.length * self.rho / 2.
            )

        self.q_c_guess = self.q_c
        self.q_h_guess = self.q_h
        self.q_guess = self.q_h

    def solve_leg_for_real(self):

        self.fsolve_output0 = fsolve(self.get_error_J, x0= self.J)
        
    def get_error_J(self, J):
        
        print "J in this step is ", self.J
        self.J = J
        self.solve_leg()

        # this is what J needs to be equal to 
        self.J_correct = (
            self.Vs / (self.R_load + self.R_internal)
            )

        # the difference between the correct J and calculated J
        self.J_error = self.J_correct - self.J

        return self.J_error

 
