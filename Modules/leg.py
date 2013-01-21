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
        self.C = 1.e7
        self.t_array = np.linspace(0., 5., 10)



        self.Vs = 1.64/256.0
        self.R_internal = 1.0/256

        self.material = 'HMS'
        self.T_h_conv = 443.
        self.T_c_conv = 323.
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

    def set_ZT(self):

        self.ZT = (
            self.alpha ** 2. * self.T_props / (self.k * self. rho)
            )

    def set_I(self):

        self.I = self.Vs / (self.R_load + self.R_internal)

    def set_constants(self):

        self.set_I()
        self.x = np.linspace(0., self.length, self.nodes)

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
        self.V = self.R_load * self.I
        self.R_internal = self.R_int_x[-1]

        # Multiply q_h and q_c by area to get rid of this error

        # self.P = self.R_load * self.I ** 2
        


        # self.P = self.I * self.V
        # print "\nPower output is ", self.P
        # # Sanity check.  q_h - q_c should be nearly equal but not
        # # exactly equal to P.  It is not exact because of spatial
        # # asymmetry in electrical resistivity along the leg.  I
        # # imported assert_approx_equal in the front matter to make
        # # this print an error if there is too much disagreement.

        # # at the beginning, they do not match so we get an error, as
        # # the iteration continues, the error disappears

        # self.P_from_heat = (self.q_h - self.q_c) * self.area
        # print "P from heat is ", self.P_from_heat

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
        J = self.I / self.area

        dT_dx = (
            (J * self.alpha * T - q) / self.k
            )
        dq_dx = (
            (self.rho * J ** 2. * (1+self.ZT)) - (J * self.alpha * q /
            self.k)
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
        J = self.I / self.area

        self.q_c = - (
            self.alpha * self.T_c * J - delta_T / self.length *
        self.k - J ** 2 * self.length * self.rho
            )

        self.q_h = - (
            self.alpha * self.T_h * J - delta_T / self.length *
            self.k + J ** 2. * self.length * self.rho / 2.
            )

        self.q_c_guess = self.q_c
        self.q_h_guess = self.q_h
        self.q_guess = self.q_h


    # ===========================================
    # everything above this is steady state 
    # everything below this is for transient only
    # ===========================================


    def solve_leg_transient(self):

        """Solves leg based on array of transient BC's."""

        self.delta_x = self.x[1] - self.x[0]

        self.y0 = self.T_x

        try: 
            self.T_xt

        except AttributeError:
            self.odeint_output = odeint(
                self.get_dTx_dt, y0=self.y0, t=self.t_array,
                full_output=1 
                )
            self.T_xt = self.odeint_output[0]

        else:
            self.y0 = self.T_xt[-1,:]
            self.odeint_output = odeint(
                self.get_dTx_dt, y0=self.y0, t=self.t_array,
                full_output=1 
                )
            self.T_xt = np.concatenate((self.T_xt, self.odeint_output[0]))

    def get_dTx_dt(self, T, t):

        """Returns derivative of array of T wrt time.
        """
        
        J = self.I/self.area
        self.dT_dt = np.zeros(T.size)
        self.q0 = np.zeros(T.size)
        self.dq_dx_ss = np.zeros(T.size)
        self.dq_dx = np.zeros(T.size)
        self.dT_dx = np.zeros(T.size)

        self.dT_dx[1:-1] = 0.5 * (T[2:] - T[:-2]) / self.delta_x  
        self.dT_dx[0] = (T[1] - T[0]) / self.delta_x
        self.dT_dx[-1] = (T[-1] - T[-2]) / self.delta_x

        for i in range(self.nodes):

            T_props = T[i]  # i for central differencing
            self.set_TEproperties(T_props)
            self.set_ZT()

            self.q0[i] = (
                J * T[i] * self.alpha - self.k * self.dT_dx[i]
                ) 

            self.dq_dx_ss[i] = (
                (self.rho * J ** 2. * (1. + self.ZT)) - J *
                self.alpha * self.q0[i] / self.k
                )

        # hot side BC, q_h
        self.q0[0] = self.U_hot * (self.T_h_conv - T[0]) 

        # cold side BC, q_c 
        self.q0[-1] = self.U_cold * (T[-1] - self.T_c_conv)

        self.dq_dx[1:-1] = (
            (self.q0[2:] - self.q0[:-2]) / (2. * self.delta_x)
            )
        self.dq_dx[0] = (
            (self.q0[1] - self.q0[0]) / self.delta_x
            )
        self.dq_dx[-1] = (
            (self.q0[-1] - self.q0[-2]) / self.delta_x
            )

        for i in range(self.nodes):

            T_props = T[i]  # i for central differencing
            self.set_TEproperties(T_props)
            self.set_ZT()

            self.dT_dt[i] = (
                1. / self.C * (-self.dq_dx[i] + self.dq_dx_ss[i])
                )

        return self.dT_dt




















































    def solve_leg_for_real(self):

        self.fsolve_output0 = fsolve(self.get_error_I, x0= self.I)
        print "Final I is ", self.I
        
    def get_error_I(self, I):
        
        self.I = I
        self.solve_leg()
        print "I in this step is ", self.I

        # this is what I needs to be equal to 
        self.I_correct = self.Vs / (self.R_load + self.R_internal)

        # the difference between the correct J and calculated J
        self.I_error = self.I_correct - self.I

        print "Error in I is ", self.I_error
        return self.I_error

 
