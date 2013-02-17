# working code 7
# Trying to solve leg with correct current this time 
#===========================================
#===========================================
#===========================================
#===========================================
#===========================================

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
        # self.C = 1.e7
        # self.C = 10. * 1.e7

        self.t_array = np.linspace(0., 5., 10)

        self.Vs = 1.64/256.0
        self.R_internal = 1.0/256

        self.material = 'HMS'
        self.T_h_conv = 500.
        self.T_c_conv = 300.
        self.U_hot = 500.
        self.U_cold = 500.
        self.R_load = 1.0/256.0
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
        self.guess_array = np.array([self.T_h, self.I])

        self.fsolve_output = fsolve(self.get_error, x0=self.guess_array)

    def get_error(self, guess):

        self.T_h = guess[0]
        self.I = guess[1]
        self.q_h_conv = (
            self.U_hot * (self.T_h_conv - self.T_h)
            )
        self.q_h = self.q_h_conv
        self.solve_leg_once(self.q_h)
        self.q_c_conv = self.U_cold * (self.T_c - self.T_c_conv)
        self.q_c_error = self.q_c - self.q_c_conv
        self.I_error = self.I_correct - self.I

        self.error = (
            np.array([self.q_c_error, self.I_error])
            )
        return self.error
    
    def solve_leg_once(self, q_h):

        self.q_h = q_h
        self.y0 = np.array([self.T_h, self.q_h, 0, 0])
        self.y = odeint(self.get_dTq_dx, y0=self.y0, t=self.x)
        
        self.T_x = self.y[:,0]
        self.q_x = self.y[:,1]
        self.Vs_x = self.y[:,2]
        self.R_x = self.y[:,3]
        
        self.T_c = self.T_x[-1]
        self.q_c = self.q_x[-1]
        
        self.Vs = self.Vs_x[0] - self.Vs_x[-1]
        self.V = self.R_load * self.I
        self.R_internal = self.R_x[-1]
        self.I_correct = (
            self.Vs / (self.R_load + self.R_internal)
            )
        #print "\nI in this run is", self.I
        #print "I_correct in this run is", self.I_correct

        # print "\nI is ", self.I
        # print "\nq_x is found to be \n", self.q_x

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
# ===========================================
# everything above this is steady state 
# everything below this is for transient only
# ===========================================
# ===========================================


    def get_dTx_dt(self, TqVsR, t):

        """Returns derivative of array of T wrt time.
        """
        # 3 for 3 terms - T_x, q_x, Vs_x, and R_x
        TqVsR.shape = (4, self.nodes)

        T = TqVsR[0,:]
        q0 = TqVsR[1,:]
        # Vs_x = TqVsR[2,:]
        # R_x = TqVsR[3,:]
        
        # ====================================
        # Need to make a 1D array of current
        # ====================================
        #self.I_transient = np.zeros(1,self.t_array.size)
        #self.I_transient[:,0] = self.I
        #J = self.I_transient / self.area
        J = self.I / self.area
        
        dT_dx = np.zeros(T.size)
        dq_dx_ss = np.zeros(T.size)
        dq_dx = np.zeros(T.size)
        dq_dt = np.zeros(T.size)
        dT_dt = np.zeros(T.size)
        dR_dt = np.zeros(T.size)
        dVs_dt = np.zeros(T.size)
        
        dT_dx[1:-1] = 0.5 * (T[2:] - T[:-2]) / self.delta_x  
        dT_dx[0] = (T[1] - T[0]) / self.delta_x
        dT_dx[-1] = (T[-1] - T[-2]) / self.delta_x

        for i in range(self.nodes):

            T_props = T[i]  # i for central differencing
            self.set_TEproperties(T_props)
            self.set_ZT()

            q0[i] = (
                J * T[i] * self.alpha - self.k * dT_dx[i]
                ) 
            # dq_dx_ss based on old q0
            dq_dx_ss[i] = (
                (self.rho * J ** 2. * (1. + self.ZT)) - J *
                self.alpha * q0[i] / self.k
                )

        #print "\nq0 hot is ", q0[0]
        q0[0] = self.U_hot * (self.T_h_conv - T[0]) 
        #print "\nNow q0 hot is ", q0[0], "\n"
        q0[-1] = self.U_cold * (T[-1] - self.T_c_conv)

        dq_dx[1:-1] = (
            (q0[2:] - q0[:-2]) / (2. * self.delta_x)
            )
        dq_dx[0] = (
            (q0[1] - q0[0]) / self.delta_x
            )
        dq_dx[-1] = (
            (q0[-1] - q0[-2]) / self.delta_x
            )

        for i in range(self.nodes):

            # need a delta_t here that changes every loop  so that log
            # scale can be used
            T_props = T[i]  # i for central differencing
            self.set_TEproperties(T_props)
            self.set_ZT()


            dT_dt[i] = (
                1. / self.C * (-dq_dx[i] + dq_dx_ss[i])
                )
            
            dVs_dt[i] = self.alpha * dT_dt[i]

#============================
#CHECK THIS FORMULA FOR dR_dt
#============================
            dR_dt[i] = (
                self.rho * self.delta_x / self.area * self.delta_t
                )

        self.return_array = (
            np.array([dT_dt, dq_dt, dVs_dt, dR_dt]).flatten()
            )

        return self.return_array

    def solve_leg_transient_once(self):

        """Solves leg based on array of transient BC's."""

        self.delta_x = self.x[1] - self.x[0]
        self.delta_t = self.t_array[1] - self.t_array[0]
        self.y0 = np.array([self.T_x, self.q_x, self.Vs_x, self.R_x]).flatten()

        try: 
            self.T_xt

        except AttributeError:
            self.odeint_output = odeint(
                self.get_dTx_dt, y0=self.y0, t=self.t_array,
                full_output=1 
                )
            self.T_xt = self.odeint_output[0]

        # doesnt really get here
        # else:
        #     self.y0 = self.T_xt[-1,:]
        #     self.odeint_output = odeint(
        #         self.get_dTx_dt, y0=self.y0, t=self.t_array,
        #         full_output=1 
        #         )
        #     self.T_xt = np.concatenate((self.T_xt, self.odeint_output[0]))
            
        print "\nDid get through all the calculations without error \n"

        self.Txt = self.T_xt[:, :self.nodes]
        # Don't need this following line anymore
        # self.qxt = self.T_xt[:, self.nodes:2*self.nodes]
        self.Vsxt = self.T_xt[:, 2*self.nodes:3*self.nodes]
        self.Rxt = self.T_xt[:, 3*self.nodes:]

        self.R_internal_transient = self.Rxt[:,-1]        
        self.Vs_transient = self.Vsxt[:,0] - self.Vsxt[:,-1]

        self.I_transient = (
            self.Vs_transient/(self.R_load +
                               self.R_internal_transient)
            )

        dT_dx = np.zeros([self.t_array.size,self.nodes])
        dT_dx[:,0] = (
            (self.Txt[:,1] - self.Txt[:,0]) / self.delta_x
            )
        dT_dx[:,1:-1] = (
            0.5 * (self.Txt[:,2:] - self.Txt[:,:-2]) / self.delta_x
            )
        dT_dx[:,-1] = (
            (self.Txt[:,-1] - self.Txt[:,-2]) / self.delta_x
            )
        
        self.qxt = np.zeros([self.t_array.size, self.nodes])
        for i in range(self.t_array.size):
            T_props = self.Txt[i,0]  # i for central differencing
            self.set_TEproperties(T_props)            
            J = (self.I_transient[i]/self.area)
            self.qxt[i,0] = (
                J * self.Txt[i,0] * self.alpha - self.k * dT_dx[i,0]
                )
            self.qxt[i,-1] = (
                J * self.Txt[i,-1] * self.alpha - self.k * dT_dx[i,-1]
                )
            self.qxt[i,1:-1] = (
                J * self.Txt[i,1:-1] * self.alpha - self.k * dT_dx[i,1:-1]
                )

































































































































# # working code 6
# # Code runs 
# #===========================================
# #===========================================
# #===========================================
# #===========================================
# #===========================================





# import types
# import numpy as np
# from scipy.integrate import odeint
# from numpy.testing import assert_approx_equal
# from scipy.optimize import fsolve

# import mat_prop
# reload (mat_prop)


# class Leg(object):

#     """ Class for individual TE leg.

#     Methods:
#     """

#     def __init__(self):

#         """ sets constants and binds methods

#         """
# #        self.C = 1.e7
#         self.C = 10. * 1.e7

#         self.t_array = np.linspace(0., 5., 10)

#         self.Vs = 1.64/256.0
#         self.R_internal = 1.0/256

#         self.material = 'HMS'
#         self.T_h_conv = 500.
#         self.T_c_conv = 300.
#         self.U_hot = 500.
#         self.U_cold = 500.
#         self.R_load = 1.0/256.0
#         self.nodes = 10
#         self.length = 1.5e-3
#         self.area = (3.0e-3) ** 2 
#         self.set_constants()
        
#         self.import_raw_property_data = ( 
#             types.MethodType(mat_prop.import_raw_property_data, self)
#             )
#         self.set_properties_v_temp = (
#             types.MethodType(mat_prop.set_properties_v_temp, self)
#             )
#         self.set_TEproperties = (
#             types.MethodType(mat_prop.set_TEproperties, self)
#             )

#     def set_ZT(self):

#         self.ZT = (
#             self.alpha ** 2. * self.T_props / (self.k * self. rho)
#             )

#     def set_I(self):

#         self.I = self.Vs / (self.R_load + self.R_internal)

#     def set_constants(self):

#         self.set_I()
#         self.x = np.linspace(0., self.length, self.nodes)

#     def solve_leg(self):
        
#         self.T_h = self.T_h_conv
#         self.T_c = self.T_c_conv

#         self.fsolve_output = fsolve(self.get_error, x0=self.T_h-1.)

#     def get_error(self, T_h):

#         self.T_h = T_h[0]
#         self.q_h_conv = (
#             self.U_hot * (self.T_h_conv - self.T_h)
#             )
#         self.q_h = self.q_h_conv
#         self.solve_leg_once(self.q_h)
#         self.q_c_conv = self.U_cold * (self.T_c - self.T_c_conv)
#         self.q_c_error = self.q_c - self.q_c_conv
        
#         return self.q_c_error
    
#     def solve_leg_once(self, q_h):

#         self.q_h = q_h
#         self.y0 = np.array([self.T_h, self.q_h, 0, 0])
#         self.y = odeint(self.get_dTq_dx, y0=self.y0, t=self.x)
        
#         self.T_x = self.y[:,0]
#         self.q_x = self.y[:,1]
#         self.Vs_x = self.y[:,2]
#         self.R_x = self.y[:,3]
        
#         self.T_c = self.T_x[-1]
#         self.q_c = self.q_x[-1]
        
#         self.Vs = self.Vs_x[0] - self.Vs_x[-1]
#         self.V = self.R_load * self.I
#         self.R_internal = self.R_x[-1]

#         # print "\nI is ", self.I
#         # print "\nq_x is found to be \n", self.q_x

#         # Multiply q_h and q_c by area to get rid of this error

#         # self.P = self.R_load * self.I ** 2
        


#         # self.P = self.I * self.V
#         # print "\nPower output is ", self.P
#         # # Sanity check.  q_h - q_c should be nearly equal but not
#         # # exactly equal to P.  It is not exact because of spatial
#         # # asymmetry in electrical resistivity along the leg.  I
#         # # imported assert_approx_equal in the front matter to make
#         # # this print an error if there is too much disagreement.

#         # # at the beginning, they do not match so we get an error, as
#         # # the iteration continues, the error disappears

#         # self.P_from_heat = (self.q_h - self.q_c) * self.area
#         # print "P from heat is ", self.P_from_heat

#         # sig_figs = 3
#         # try:
#         #     assert_approx_equal(self.P, self.P_from_heat, sig_figs)
#         # except AssertionError:
#         #     print "\nPower from q_h - q_c and I**2 * R disagree."
#         #     print "Consider reducing sig_figs under solve_leg_once"
#         #     print "in leg.py if you think this is an error."

#     def get_dTq_dx(self, Tq, x):

#         T = Tq[0]
#         q = Tq[1]

#         self.set_TEproperties(T)
#         self.set_ZT()
#         J = self.I / self.area

#         dT_dx = (
#             (J * self.alpha * T - q) / self.k
#             )

#         dq_dx = (
#             (self.rho * J ** 2. * (1+self.ZT)) - (J * self.alpha * q /
#             self.k)
#             )
#         dVs_dx = (self.alpha * dT_dx)
#         dR_dx = (self.rho / self.area)

#         return dT_dx, dq_dx, dVs_dx, dR_dx
    
#     def set_q_guess(self):

#         """Sets guess for q_c to be used by iterative solutions.
#         """

#         self.T_props = 0.5 * (self.T_h + self.T_c)
#         self.set_TEproperties(T_props=self.T_props)
#         delta_T = self.T_h - self.T_c
#         J = self.I / self.area

#         self.q_c = - (
#             self.alpha * self.T_c * J - delta_T / self.length *
#             self.k - J ** 2 * self.length * self.rho
#             )

#         self.q_h = - (
#             self.alpha * self.T_h * J - delta_T / self.length *
#             self.k + J ** 2. * self.length * self.rho / 2.
#             )

#         self.q_c_guess = self.q_c
#         self.q_h_guess = self.q_h
#         self.q_guess = self.q_h

# # ===========================================
# # ===========================================
# # everything above this is steady state 
# # everything below this is for transient only
# # ===========================================
# # ===========================================


#     def get_dTx_dt(self, TqVsR, t):

#         """Returns derivative of array of T wrt time.
#         """
#         # 3 for 3 terms - T_x, q_x, Vs_x, and R_x
#         TqVsR.shape = (4, self.nodes)

#         T = TqVsR[0,:]
#         # q0 = TqVsR[1,:]
#         # Vs_x = TqVsR[2,:]
#         # R_x = TqVsR[3,:]
        
#         # ====================================
#         # Need to make a 1D array of current
#         # ====================================
#         #self.I_transient = np.zeros(1,self.t_array.size)
#         #self.I_transient[:,0] = self.I
#         #J = self.I_transient / self.area
#         J = self.I / self.area
        
#         dT_dx = np.zeros(T.size)
#         q0 = np.zeros(T.size)
#         dq_dx_ss = np.zeros(T.size)
#         dq_dx = np.zeros(T.size)
#         dq_dt = np.zeros(T.size)
#         dT_dt = np.zeros(T.size)
#         dR_dt = np.zeros(T.size)
#         dVs_dt = np.zeros(T.size)
        
#         dT_dx[1:-1] = 0.5 * (T[2:] - T[:-2]) / self.delta_x  
#         dT_dx[0] = (T[1] - T[0]) / self.delta_x
#         dT_dx[-1] = (T[-1] - T[-2]) / self.delta_x

#         for i in range(self.nodes):

#             T_props = T[i]  # i for central differencing
#             self.set_TEproperties(T_props)
#             self.set_ZT()

#             q0[i] = (
#                 J * T[i] * self.alpha - self.k * dT_dx[i]
#                 ) 
#             # dq_dx_ss based on old q0
#             dq_dx_ss[i] = (
#                 (self.rho * J ** 2. * (1. + self.ZT)) - J *
#                 self.alpha * q0[i] / self.k
#                 )

#         # update q0
#         # hot side BC, q_h
#         q0[0] = self.U_hot * (self.T_h_conv - T[0]) 

#         # cold side BC, q_c 
#         q0[-1] = self.U_cold * (T[-1] - self.T_c_conv)

#         # this is dq_dx
#         dq_dx[1:-1] = (
#             (q0[2:] - q0[:-2]) / (2. * self.delta_x)
#             )
#         dq_dx[0] = (
#             (q0[1] - q0[0]) / self.delta_x
#             )
#         dq_dx[-1] = (
#             (q0[-1] - q0[-2]) / self.delta_x
#             )

#         for i in range(self.nodes):

#             # need a delta_t here that changes every loop  so that log
#             # scale can be used
#             T_props = T[i]  # i for central differencing
#             self.set_TEproperties(T_props)
#             self.set_ZT()


#             dT_dt[i] = (
#                 1. / self.C * (-dq_dx[i] + dq_dx_ss[i])
#                 )
            
#             dVs_dt[i] = self.alpha * dT_dt[i]

#             dR_dt[i] = (
#                 self.rho * self.delta_x / (self.area * self.delta_t)
#                 )

#         self.return_array = (
#             np.array([dT_dt, dq_dt, dVs_dt, dR_dt]).flatten()
#             )

#         return self.return_array

#     def solve_leg_transient_once(self):

#         """Solves leg based on array of transient BC's."""

#         self.delta_x = self.x[1] - self.x[0]
#         self.delta_t = self.t_array[1] - self.t_array[0]
#         self.y0 = np.array([self.T_x, self.q_x, self.Vs_x, self.R_x]).flatten()

#         try: 
#             self.T_xt
#         # basically this command is being run, not the bottom one
#         except AttributeError:
#             self.odeint_output = odeint(
#                 self.get_dTx_dt, y0=self.y0, t=self.t_array,
#                 full_output=1 
#                 )
#             self.T_xt = self.odeint_output[0]

#         # doesnt really get here
#         # else:
#         #     self.y0 = self.T_xt[-1,:]
#         #     self.odeint_output = odeint(
#         #         self.get_dTx_dt, y0=self.y0, t=self.t_array,
#         #         full_output=1 
#         #         )
#         #     self.T_xt = np.concatenate((self.T_xt, self.odeint_output[0]))
            
#         print "\nDid get through all the calculations without error \n"

#         self.Txt = self.T_xt[:, :self.nodes]
#         self.qxt = self.T_xt[:, self.nodes:2*self.nodes]
#         self.Vsxt = self.T_xt[:, 2*self.nodes:3*self.nodes]
#         self.Rxt = self.T_xt[:, 3*self.nodes:]

#         self.R_internal_transient = self.Rxt[:,-1]        
#         self.Vs_transient = self.Vsxt[:,0] - self.Vsxt[:,-1]

#         self.I_transient = (
#             self.Vs_transient/(self.R_load +
#                                self.R_internal_transient)
#             )

#         self.q_h_xt = self.qxt[:,0]
#         self.q_c_xt = self.qxt[:,-1]

#         dT_dx = np.zeros([self.t_array.size,self.nodes])
#         dT_dx[:,0] = (
#             (self.Txt[:,1] - self.Txt[:,0]) / self.delta_x
#             )
#         dT_dx[:,1:-1] = (
#             0.5 * (self.Txt[:,2:] - self.Txt[:,:-2]) / self.delta_x
#             )
#         dT_dx[:,-1] = (
#             (self.Txt[:,-1] - self.Txt[:,-2]) / self.delta_x
#             )
        
#         self.q = np.zeros([self.t_array.size, self.nodes])
#         for i in range(self.t_array.size):
#             T_props = self.Txt[i,0]  # i for central differencing
#             self.set_TEproperties(T_props)            
#             J = self.I_transient[i]
#             self.q[i,0] = (
#                 J * self.Txt[i,0] * self.alpha - self.k * dT_dx[i,0]
#                 )
#             self.q[i,-1] = (
#                 J * self.Txt[i,-1] * self.alpha - self.k * dT_dx[i,-1]
#                 )
#             self.q[i,1:-1] = (
#                 J * self.Txt[i,1:-1] * self.alpha - self.k * dT_dx[i,1:-1]
#                 )

#         # self.q_h_transient = self.q[:,0]
#         # self.q_c_transient = self.q[:,-1]

#         # self.q_h_transient_correct = (
#         #     self.U_hot * (self.T_h_conv - self.Txt[:,0])
#         #     )

#         # self.q_h_transient_correct = (
#         #     self.U_cold * (self.Txt[:,-1] - self.T_c_conv)
#         #     )
        
#         # # self.T_x = self.Txt[self.t_array.size,:]
#         # # self.q_x = self.q[self.t_array.size,:]
#         # # self.Vs_x = self.Vsxt[self.t_array.size,:]
#         # # self.R_x = self.Rxt[self.t_array.size,:]









# # ======================================
# #        self.dT_dt[0,:] = 
# # ======================================        

# # LAYOUT OF THE PROCEDURE

# # solve_leg_transient()

# # solve_leg_transient_once()
# # fsolve(get_leg_transient_error, guess = )

# # get_leg_transient_error()
# # gets the guess
# # use the first line of guess as T
# # error is the error in hot side convection heat flux based on the new
# # temperature distribution that was calculated


















































































































# # working code 5
# # Code runs 
# #===========================================
# #===========================================
# #===========================================
# # This is able to calculate Vs, R, and wrong q since the BC is forced
# # way ahead, from now on, I am not enforcing a q0[0] and q0[-1] while
# # solving, I will solve it and then force the BCs to get a solution
# # for transient
# #===========================================
# #===========================================





# import types
# import numpy as np
# from scipy.integrate import odeint
# from numpy.testing import assert_approx_equal
# from scipy.optimize import fsolve

# import mat_prop
# reload (mat_prop)


# class Leg(object):

#     """ Class for individual TE leg.

#     Methods:
#     """

#     def __init__(self):

#         """ sets constants and binds methods

#         """
#         self.C = 1.e7
#         self.t_array = np.linspace(0., 5., 10)

#         self.Vs = 1.64/256.0
#         self.R_internal = 1.0/256

#         self.material = 'HMS'
#         self.T_h_conv = 500.
#         self.T_c_conv = 300.
#         self.U_hot = 54.e2
#         self.U_cold = 253.e3
#         self.R_load = 1.0/256.0
#         self.nodes = 10
#         self.length = 1.5e-3
#         self.area = (3.0e-3) ** 2 
#         self.set_constants()
        
#         self.import_raw_property_data = ( 
#             types.MethodType(mat_prop.import_raw_property_data, self)
#             )
#         self.set_properties_v_temp = (
#             types.MethodType(mat_prop.set_properties_v_temp, self)
#             )
#         self.set_TEproperties = (
#             types.MethodType(mat_prop.set_TEproperties, self)
#             )

#     def set_ZT(self):

#         self.ZT = (
#             self.alpha ** 2. * self.T_props / (self.k * self. rho)
#             )

#     def set_I(self):

#         self.I = self.Vs / (self.R_load + self.R_internal)

#     def set_constants(self):

#         self.set_I()
#         self.x = np.linspace(0., self.length, self.nodes)

#     def solve_leg(self):
        
#         self.T_h = self.T_h_conv
#         self.T_c = self.T_c_conv

#         self.fsolve_output = fsolve(self.get_error, x0=self.T_h-1.)

#     def get_error(self, T_h):

#         self.T_h = T_h[0]
#         self.q_h_conv = (
#             self.U_hot * (self.T_h_conv - self.T_h)
#             )
#         self.q_h = self.q_h_conv
#         self.solve_leg_once(self.q_h)
#         self.q_c_conv = self.U_cold * (self.T_c - self.T_c_conv)
#         self.q_c_error = self.q_c - self.q_c_conv
        
#         return self.q_c_error
    
#     def solve_leg_once(self, q_h):

#         self.q_h = q_h
#         self.y0 = np.array([self.T_h, self.q_h, 0, 0])
#         self.y = odeint(self.get_dTq_dx, y0=self.y0, t=self.x)
        
#         self.T_x = self.y[:,0]
#         self.q_x = self.y[:,1]
#         self.Vs_x = self.y[:,2]
#         self.R_x = self.y[:,3]
        
#         self.T_c = self.T_x[-1]
#         self.q_c = self.q_x[-1]
        
#         self.Vs = self.Vs_x[0] - self.Vs_x[-1]
#         self.V = self.R_load * self.I
#         self.R_internal = self.R_x[-1]

#         # print "\nI is ", self.I
#         # print "\nq_x is found to be \n", self.q_x

#         # Multiply q_h and q_c by area to get rid of this error

#         # self.P = self.R_load * self.I ** 2
        


#         # self.P = self.I * self.V
#         # print "\nPower output is ", self.P
#         # # Sanity check.  q_h - q_c should be nearly equal but not
#         # # exactly equal to P.  It is not exact because of spatial
#         # # asymmetry in electrical resistivity along the leg.  I
#         # # imported assert_approx_equal in the front matter to make
#         # # this print an error if there is too much disagreement.

#         # # at the beginning, they do not match so we get an error, as
#         # # the iteration continues, the error disappears

#         # self.P_from_heat = (self.q_h - self.q_c) * self.area
#         # print "P from heat is ", self.P_from_heat

#         # sig_figs = 3
#         # try:
#         #     assert_approx_equal(self.P, self.P_from_heat, sig_figs)
#         # except AssertionError:
#         #     print "\nPower from q_h - q_c and I**2 * R disagree."
#         #     print "Consider reducing sig_figs under solve_leg_once"
#         #     print "in leg.py if you think this is an error."

#     def get_dTq_dx(self, Tq, x):

#         T = Tq[0]
#         q = Tq[1]

#         self.set_TEproperties(T)
#         self.set_ZT()
#         J = self.I / self.area

#         dT_dx = (
#             (J * self.alpha * T - q) / self.k
#             )

#         dq_dx = (
#             (self.rho * J ** 2. * (1+self.ZT)) - (J * self.alpha * q /
#             self.k)
#             )
#         dVs_dx = (self.alpha * dT_dx)
#         dR_dx = (self.rho / self.area)

#         return dT_dx, dq_dx, dVs_dx, dR_dx
    
#     def set_q_guess(self):

#         """Sets guess for q_c to be used by iterative solutions.
#         """

#         self.T_props = 0.5 * (self.T_h + self.T_c)
#         self.set_TEproperties(T_props=self.T_props)
#         delta_T = self.T_h - self.T_c
#         J = self.I / self.area

#         self.q_c = - (
#             self.alpha * self.T_c * J - delta_T / self.length *
#             self.k - J ** 2 * self.length * self.rho
#             )

#         self.q_h = - (
#             self.alpha * self.T_h * J - delta_T / self.length *
#             self.k + J ** 2. * self.length * self.rho / 2.
#             )

#         self.q_c_guess = self.q_c
#         self.q_h_guess = self.q_h
#         self.q_guess = self.q_h

# # ===========================================
# # ===========================================
# # everything above this is steady state 
# # everything below this is for transient only
# # ===========================================
# # ===========================================


#     def get_dTx_dt(self, TqVsR, t):

#         """Returns derivative of array of T wrt time.
#         """
#         # 3 for 3 terms - T_x, q_x, Vs_x, and R_x
#         TqVsR.shape = (4, self.nodes)

#         T = TqVsR[0,:]
#         # Vs_x = TqVsR[2,:]
#         # R_x = TqVsR[3,:]
        
#         # ====================================
#         # Need to make a 1D array of current
#         # ====================================
#         #self.I_transient = np.zeros(1,self.t_array.size)
#         #self.I_transient[:,0] = self.I
#         #J = self.I_transient / self.area
#         J = self.I / self.area
        
#         dT_dx = np.zeros(T.size)
#         q0 = np.zeros(T.size)
#         dq_dx_ss = np.zeros(T.size)
#         dq_dx = np.zeros(T.size)
#         dq_dt = np.zeros(T.size)
#         dT_dt = np.zeros(T.size)
#         dR_dt = np.zeros(T.size)
#         dVs_dt = np.zeros(T.size)
        
#         dT_dx[1:-1] = 0.5 * (T[2:] - T[:-2]) / self.delta_x  
#         dT_dx[0] = (T[1] - T[0]) / self.delta_x
#         dT_dx[-1] = (T[-1] - T[-2]) / self.delta_x

#         for i in range(self.nodes):

#             T_props = T[i]  # i for central differencing
#             self.set_TEproperties(T_props)
#             self.set_ZT()

#             q0[i] = (
#                 J * T[i] * self.alpha - self.k * dT_dx[i]
#                 ) 
#             # dq_dx_ss based on old q0
#             dq_dx_ss[i] = (
#                 (self.rho * J ** 2. * (1. + self.ZT)) - J *
#                 self.alpha * q0[i] / self.k
#                 )

#         # update q0
#         # hot side BC, q_h
#         q0[0] = self.U_hot * (self.T_h_conv - T[0]) 

#         # cold side BC, q_c 
#         q0[-1] = self.U_cold * (T[-1] - self.T_c_conv)

#         # this is dq_dx
#         dq_dx[1:-1] = (
#             (q0[2:] - q0[:-2]) / (2. * self.delta_x)
#             )
#         dq_dx[0] = (
#             (q0[1] - q0[0]) / self.delta_x
#             )
#         dq_dx[-1] = (
#             (q0[-1] - q0[-2]) / self.delta_x
#             )

#         for i in range(self.nodes):

#             # need a delta_t here that changes every loop  so that log
#             # scale can be used
#             T_props = T[i]  # i for central differencing
#             self.set_TEproperties(T_props)
#             self.set_ZT()

#             # This dq_dt equation is good in node direction but it
#             # does not account for the time wise direction
#             # have to work with 2D arrays

#             dq_dt[i] = (
#                 ((self.I * self.alpha * dT_dt[i]) / self.area) -
#                 (self.k/self.delta_x) * (dT_dt[i] - dT_dt[i-1])
#                 )            

#             dT_dt[i] = (
#                 1. / self.C * (-dq_dx[i] + dq_dx_ss[i])
#                 )
            
#             dVs_dt[i] = self.alpha * dT_dt[i]

#             dR_dt[i] = (
#                 self.rho * self.delta_x / self.area * self.delta_t
#                 )

#         # dq_dt[0] = (
#         #     ((self.I * self.alpha * dT_dt[0]) / self.area) -
#         #     (self.k/self.delta_x) * (dT_dt[1] - dT_dt[0])
#         #     )            

#         # dq_dt[1:-1] = (
#         #     ((self.I * self.alpha * dT_dt[1:-1]) / self.area) -
#         #     (self.k/(2*self.delta_x)) * (dT_dt[2:] - dT_dt[:-2])
#         #     )            
        
#         # dq_dt[0] = (
#         #     ((self.I * self.alpha * dT_dt[0]) / self.area) -
#         #     (self.k/self.delta_x) * (dT_dt[1] - dT_dt[0])
#         #     )            

#         # dq_dt[-1] = (
#         #     ((self.I * self.alpha * dT_dt[-1]) / self.area) -
#         #     (self.k/self.delta_x) * (dT_dt[-2] - dT_dt[-1])
#         #     )            

#         self.return_array = (
#             np.array([dT_dt, dq_dt, dVs_dt, dR_dt]).flatten()
#             )

#         return self.return_array

#     def solve_leg_transient_once(self):

#         """Solves leg based on array of transient BC's."""

#         self.delta_x = self.x[1] - self.x[0]
#         self.delta_t = self.t_array[1] - self.t_array[0]
#         self.y0 = np.array([self.T_x, self.q_x, self.Vs_x, self.R_x]).flatten()

#         try: 
#             self.T_xt
#         # basically this command is being run, not the bottom one
#         except AttributeError:
#             self.odeint_output = odeint(
#                 self.get_dTx_dt, y0=self.y0, t=self.t_array,
#                 full_output=1 
#                 )
#             self.T_xt = self.odeint_output[0]

#         # doesnt really get here
#         else:
#             self.y0 = self.T_xt[-1,:]
#             self.odeint_output = odeint(
#                 self.get_dTx_dt, y0=self.y0, t=self.t_array,
#                 full_output=1 
#                 )
#             self.T_xt = np.concatenate((self.T_xt, self.odeint_output[0]))
            
#         print "\nDid get through all the calculations without error \n"

#         self.Txt = self.T_xt[:, :self.nodes]
#         self.qxt = self.T_xt[:, self.nodes:2*self.nodes]
#         self.Vsxt = self.T_xt[:, 2*self.nodes:3*self.nodes]
#         self.Rxt = self.T_xt[:, 3*self.nodes:]

#         self.R_internal_transient = self.Rxt[:,-1]        
#         self.Vs_transient = self.Vsxt[:,0] - self.Vsxt[:,-1]

#         self.I_transient = (
#             self.Vs_transient/(self.R_load +
#                                self.R_internal_transient)
#             )

#         self.q_h_xt = self.qxt[:,0]
#         self.q_c_xt = self.qxt[:,-1]

#         dT_dx = np.zeros([self.t_array.size,self.nodes])
#         dT_dx[:,0] = (
#             (self.Txt[:,1] - self.Txt[:,0]) / self.delta_x
#             )
#         dT_dx[:,1:-1] = (
#             0.5 * (self.Txt[:,2:] - self.Txt[:,:-2]) / self.delta_x
#             )
#         dT_dx[:,-1] = (
#             (self.Txt[:,-1] - self.Txt[:,-2]) / self.delta_x
#             )
        
#         self.q = np.zeros([self.t_array.size, self.nodes])
#         for i in range(self.t_array.size):
#             T_props = self.Txt[i,0]  # i for central differencing
#             self.set_TEproperties(T_props)            
#             J = self.I_transient[i]
#             self.q[i,0] = (
#                 J * self.Txt[i,0] * self.alpha - self.k * dT_dx[i,0]
#                 )
#             self.q[i,-1] = (
#                 J * self.Txt[i,-1] * self.alpha - self.k * dT_dx[i,-1]
#                 )
#             self.q[i,1:-1] = (
#                 J * self.Txt[i,1:-1] * self.alpha - self.k * dT_dx[i,1:-1]
#                 )

#             # self.q[i,0] = (
#             #     J * self.Txt[i,0] * self.alpha
#             #     )
#             # self.q[i,-1] = (
#             #     J * self.Txt[i,-1] * self.alpha
#             #     )
#             # self.q[i,1:-1] = (
#             #     J * self.Txt[i,1:-1] * self.alpha
#             #     )


#         # This following formula is wrong
#         # self.dq_dt[:,-1] = (
#         #     (self.Txt[2:] - self.Txt[:-2]) / (2. * self.delta_t)
#         #     )
#         # self.dq_dt[:,0] = (
#         #     (self.Txt[1] - self.Txt[0]) / self.delta_t
#         #     )
#         # self.dq_dt[:,1:-1] = (
#         #     (self.Txt[-1] - self.Txt[-2]) / self.delta_t
#         #     )

# # ======================================
# #        self.dT_dt[0,:] = 
# # ======================================        

# # LAYOUT OF THE PROCEDURE

# # solve_leg_transient()

# # solve_leg_transient_once()
# # fsolve(get_leg_transient_error, guess = )

# # get_leg_transient_error()
# # gets the guess
# # use the first line of guess as T
# # error is the error in hot side convection heat flux based on the new
# # temperature distribution that was calculated



































# working code 4
# Code runs 
#===========================================
#===========================================
#===========================================
# This is able to calculate Vs and R for transient
# but not the error in current yet or power
#===========================================
#===========================================



    # # def get_dTx_dt(self, T, t):
    # def get_dTx_dt(self, TqVsR, t):

    #     """Returns derivative of array of T wrt time.
    #     """
    #     # 3 for 3 terms - T_x, Vs_x, and R_x
    #     TqVsR.shape = (3, self.nodes)

    #     T = TqVsR[0,:]
    #     # Vs_x = TqVsR[2,:]
    #     # R_x = TqVsR[3,:]
        
    #     # ====================================
    #     # Need to make a 1D array of current
    #     # ====================================
        
    #     J = self.I/self.area

    #     dT_dx = np.zeros(T.size)
    #     q0 = np.zeros(T.size)
    #     dq_dx_ss = np.zeros(T.size)
    #     dq_dx = np.zeros(T.size)
    #     dT_dt = np.zeros(T.size)
    #     dR_dt = np.zeros(T.size)
    #     dVs_dt = np.zeros(T.size)
        
    #     # This is dT_dx
    #     dT_dx[1:-1] = 0.5 * (T[2:] - T[:-2]) / self.delta_x  
    #     dT_dx[0] = (T[1] - T[0]) / self.delta_x
    #     dT_dx[-1] = (T[-1] - T[-2]) / self.delta_x

    #     for i in range(self.nodes):

    #         T_props = T[i]  # i for central differencing
    #         self.set_TEproperties(T_props)
    #         self.set_ZT()

    #         q0[i] = (
    #             J * T[i] * self.alpha - self.k * dT_dx[i]
    #             ) 
    #         # dq_dx_ss based on old q0
    #         dq_dx_ss[i] = (
    #             (self.rho * J ** 2. * (1. + self.ZT)) - J *
    #             self.alpha * q0[i] / self.k
    #             )
    #     # update q0
    #     # hot side BC, q_h
    #     q0[0] = self.U_hot * (self.T_h_conv - T[0]) 

    #     # cold side BC, q_c 
    #     q0[-1] = self.U_cold * (T[-1] - self.T_c_conv)

    #     # this is dq_dx
    #     dq_dx[1:-1] = (
    #         (q0[2:] - q0[:-2]) / (2. * self.delta_x)
    #         )
    #     dq_dx[0] = (
    #         (q0[1] - q0[0]) / self.delta_x
    #         )
    #     dq_dx[-1] = (
    #         (q0[-1] - q0[-2]) / self.delta_x
    #         )

    #     for i in range(self.nodes):

    #         T_props = T[i]  # i for central differencing
    #         self.set_TEproperties(T_props)
    #         self.set_ZT()

    #         # this is dT_dt
    #         dT_dt[i] = (
    #             1. / self.C * (-dq_dx[i] + dq_dx_ss[i])
    #             )

    #         dVs_dt[i] = self.alpha * dT_dt[i]

    #         dR_dt[i] = (
    #             self.rho * self.delta_x / self.area * self.delta_t
    #             )

    #     self.return_array = (
    #         np.array([dT_dt, dVs_dt, dR_dt]).flatten()
    #         )
    #     # print "\nreturn array is \n", self.return_array
    #     return self.return_array

    #     # return dT_dt

    # def solve_leg_transient_once(self):

    #     """Solves leg based on array of transient BC's."""

    #     self.delta_x = self.x[1] - self.x[0]
    #     self.delta_t = self.t_array[1] - self.t_array[0]
    #     self.y0 = np.array([self.T_x, self.Vs_x, self.R_x]).flatten()

    #     try: 
    #         self.T_xt
    #     # basically this command is being run, not the bottom one
    #     except AttributeError:
    #         self.odeint_output = odeint(
    #             self.get_dTx_dt, y0=self.y0, t=self.t_array,
    #             full_output=1 
    #             )
    #         self.T_xt = self.odeint_output[0]

    #     # doesnt really get here
    #     else:
    #         self.y0 = self.T_xt[-1,:]
    #         self.odeint_output = odeint(
    #             self.get_dTx_dt, y0=self.y0, t=self.t_array,
    #             full_output=1 
    #             )
    #         self.T_xt = np.concatenate((self.T_xt, self.odeint_output[0]))
            
    #     print "\nDid get through all the calculations without error \n"

    #     self.Txt = self.T_xt[:, :self.nodes]
    #     self.Vsxt = self.T_xt[:, self.nodes:2*self.nodes]
    #     self.Rxt = self.T_xt[:, 2*self.nodes:]

    #     self.R_internal_transient = self.Rxt[:,-1]        
    #     self.Vs_transient = self.Vsxt[:,0] - self.Vsxt[:,-1]
    #     # self.Vs_power = self.Vs_transient * self.I_transient

        
    # # def get_error_transient(self, guess_array):
    # #     """ get transient I error """

    # def solve_transient_leg(self):
    #     """ transient solution for one leg """

    #     self.I_transient = np.zeros(self.t_array.size)
    #     self.I_transient[:] = self.I
    #     self.fsolve_output_transient = fsolve(self.get_error_transient, x0=self.I_transient)

    # def get_error_transient(self, guess_array):
    #     """get transient error in current"""

    #     self.solve_leg_transient_once()

    #     self.I_correct_transinet = (
    #         self.Vs_transient / (self.R_load + self.R_internal_transient)
    #         )
        
    #     self.I_error_transient = (
    #         self.I_transient - self.I_correct_transinet
    #         )

    #     return self.I_error_transient

        








































# working code 3
# Code runs but the process is not right 
#===========================================
#===========================================
#===========================================
# The following code is correct, if any screw up happens, copy and
# paste this in place of the codes given above. At least the following
# ones work. 
#===========================================
#===========================================






    # # def get_dTx_dt(self, T, t):
    # def get_dTx_dt(self, TqVsR, t):

    #     """Returns derivative of array of T wrt time.
    #     """
    #     # ===========================================
    #     # VERY IMPORTANT - all of TqVsR doesn't need to be used
    #     # ===========================================
    #     # 4 for 4 terms - T_x, q_x, Vs_x, and R_x
    #     TqVsR.shape = (4, self.nodes)

    #     T = TqVsR[0,:]
    #     # following are unnecessary
    #     # q_x = TqVsR[1,:]
    #     # Vs_x = TqVsR[2,:]
    #     # T_x = TqVsR[3,:]

    #     J = self.I/self.area

    #     dT_dx = np.zeros(T.size)
    #     q0 = np.zeros(T.size)
    #     dq_dx_ss = np.zeros(T.size)
    #     dq_dx = np.zeros(T.size)
    #     dT_dt = np.zeros(T.size)
    #     dVs_dx = np.zeros(T.size)
    #     dR_dx = np.zeros(T.size)
        
    #     # This is dT_dx
    #     dT_dx[1:-1] = 0.5 * (T[2:] - T[:-2]) / self.delta_x  
    #     dT_dx[0] = (T[1] - T[0]) / self.delta_x
    #     dT_dx[-1] = (T[-1] - T[-2]) / self.delta_x

    #     for i in range(self.nodes):

    #         T_props = T[i]  # i for central differencing
    #         self.set_TEproperties(T_props)
    #         self.set_ZT()

    #         dVs_dx[i] = self.alpha * dT_dx[i]

    #         dR_dx[i] = self.rho / self.area

    #         q0[i] = (
    #             J * T[i] * self.alpha - self.k * dT_dx[i]
    #             ) 
    #         # dq_dx_ss based on old q0
    #         dq_dx_ss[i] = (
    #             (self.rho * J ** 2. * (1. + self.ZT)) - J *
    #             self.alpha * q0[i] / self.k
    #             )
    #     # update q0
    #     # hot side BC, q_h
    #     q0[0] = self.U_hot * (self.T_h_conv - T[0]) 

    #     # cold side BC, q_c 
    #     q0[-1] = self.U_cold * (T[-1] - self.T_c_conv)

    #     # this is dq_dx
    #     dq_dx[1:-1] = (
    #         (q0[2:] - q0[:-2]) / (2. * self.delta_x)
    #         )
    #     dq_dx[0] = (
    #         (q0[1] - q0[0]) / self.delta_x
    #         )
    #     dq_dx[-1] = (
    #         (q0[-1] - q0[-2]) / self.delta_x
    #         )

    #     for i in range(self.nodes):

    #         T_props = T[i]  # i for central differencing
    #         self.set_TEproperties(T_props)
    #         self.set_ZT()

    #         # this is dT_dt
    #         dT_dt[i] = (
    #             1. / self.C * (-dq_dx[i] + dq_dx_ss[i])
    #             )

    #     # ==================================================
    #     # ALL OF THIS IS WRONG SINCE I NEED TO RETURN t DERIVATES
    #     # AND NOT x DERIVATES
    #     # ==================================================

    #     self.return_array = (
    #         np.array([dT_dt, dq_dx, dVs_dx, dR_dx]).flatten()
    #         )
    #     return self.return_array        
    #     # return dT_dt

    

    # def solve_leg_transient_once(self):

    #     """Solves leg based on array of transient BC's."""

    #     self.delta_x = self.x[1] - self.x[0]
    #     self.y0 = np.array([self.T_x, self.q_x, self.Vs_x, self.R_x]).flatten()
    #     print "\nAfter flatten, y0 reads like \n", self.y0
    #     print "\nsize and shape  of y0 is \n"
    #     print self.y0.size
    #     print self.y0.shape

    #     # self.y0 = self.T_x

    #     try: 
    #         self.T_xt
    #     # basically this command is being run, not the bottom one
    #     except AttributeError:
    #         self.odeint_output = odeint(
    #             self.get_dTx_dt, y0=self.y0, t=self.t_array,
    #             full_output=1 
    #             )
    #         self.T_xt = self.odeint_output[0]            

    #     # doesnt really get here
    #     else:
    #         self.y0 = self.T_xt[-1,:]
    #         self.odeint_output = odeint(
    #             self.get_dTx_dt, y0=self.y0, t=self.t_array,
    #             full_output=1 
    #             )
    #         self.T_xt = np.concatenate((self.T_xt, self.odeint_output[0]))
            
    #     print "\nDid get through all the calculations without error \n"
        
    #     print "\nT_xt array is \n", self.T_xt
        
    #     # ===================================
    #     # This one seems to work
    #     self.Txt = self.T_xt[:, :self.nodes]
    #     self.qxt = self.T_xt[:, self.nodes:2*self.nodes]
    #     self.Vsxt = self.T_xt[:, 2*self.nodes:3*self.nodes]
    #     self.Rxt = self.T_xt[:, 3*self.nodes:]

    #     print "\nT's are \n", self.Txt
    #     print "\nq's are \n", self.qxt
    #     print "\nVs's are \n", self.Vsxt
    #     print "\nR's are \n", self.Rxt












# working code 2 
#===========================================
#===========================================
#===========================================
# The following code is correct, if any screw up happens, copy and
# paste this in place of the codes given above. At least the following
# ones work. 
#===========================================
#===========================================




    # def get_dTx_dt(self, T, t):
    # # def get_dTx_dt(self, TqVsR, t):

    #     """Returns derivative of array of T wrt time.
    #     """
    #     # =============================================
    #     # make sure that a 2D array T is being returned
    #     # 2D array flattened into a very long 1D array
    #     # need to divide this T by number of varables
    #     # =============================================

    #     # T is basically the same as y0
    #     # 4 for 4 terms - T_x, q_x, Vs_x, and R_x
    #     # TqVsR.shape = (4, self.nodes)
    #     # print "\nafter reshape, TqVsR looks like \n", TqVsR
        
    #     # ===========================================
    #     # VERY IMPORTANT - all of TqVsR doesn't need to be used
    #     # ===========================================

    #     # T = TqVsR[0,:]
    #     # q_x = TqVsR[1,:]
    #     # Vs_x = TqVsR[2,:]
    #     # T_x = TqVsR[3,:]

    #     J = self.I/self.area

    #     dT_dx = np.zeros(T.size)
    #     q0 = np.zeros(T.size)
    #     dq_dx_ss = np.zeros(T.size)
    #     dq_dx = np.zeros(T.size)
    #     dT_dt = np.zeros(T.size)

    #     # dR_dx = np.zeros(T.size)
    #     # dVs_dx = np.zeros(T.size)

    #     # dT_dx = np.zeros(u)
    #     # q0 = np.zeros(u)
    #     # dq_dx_ss = np.zeros(u)
    #     # dq_dx = np.zeros(u)
    #     # dT_dt = np.zeros(u)

    #     # dR_dx = np.zeros(u)
    #     # dVs_dx = np.zeros(u)
        
    #     # This is dT_dx
    #     dT_dx[1:-1] = 0.5 * (T[2:] - T[:-2]) / self.delta_x  
    #     dT_dx[0] = (T[1] - T[0]) / self.delta_x
    #     dT_dx[-1] = (T[-1] - T[-2]) / self.delta_x

    #     for i in range(self.nodes):

    #         T_props = T[i]  # i for central differencing
    #         self.set_TEproperties(T_props)
    #         self.set_ZT()

    #         # dR_dx[i] = self.rho / self.area

    #         # dVs_dx[i] = self.alpha * dT_dx[i]

    #         q0[i] = (
    #             J * T[i] * self.alpha - self.k * dT_dx[i]
    #             ) 
    #         # dq_dx_ss based on old q0
    #         dq_dx_ss[i] = (
    #             (self.rho * J ** 2. * (1. + self.ZT)) - J *
    #             self.alpha * q0[i] / self.k
    #             )
    #     # update q0
    #     # hot side BC, q_h
    #     q0[0] = self.U_hot * (self.T_h_conv - T[0]) 

    #     # cold side BC, q_c 
    #     q0[-1] = self.U_cold * (T[-1] - self.T_c_conv)

    #     # this is dq_dx
    #     dq_dx[1:-1] = (
    #         (q0[2:] - q0[:-2]) / (2. * self.delta_x)
    #         )
    #     dq_dx[0] = (
    #         (q0[1] - q0[0]) / self.delta_x
    #         )
    #     dq_dx[-1] = (
    #         (q0[-1] - q0[-2]) / self.delta_x
    #         )

    #     for i in range(self.nodes):

    #         T_props = T[i]  # i for central differencing
    #         self.set_TEproperties(T_props)
    #         self.set_ZT()

    #         # this is dT_dt
    #         dT_dt[i] = (
    #             1. / self.C * (-dq_dx[i] + dq_dx_ss[i])
    #             )

    #     # self.return_array = (
    #     #     np.array([dT_dt, dq_dx, dVs_dx, dR_dx]).flatten()
    #     #     )
        
    #     # return self.return_array

    #     return dT_dt

    # def solve_leg_transient_once(self):

    #     """Solves leg based on array of transient BC's."""

    #     self.delta_x = self.x[1] - self.x[0]
    #     self.y0 = np.array([self.T_x, self.q_x, self.Vs_x, self.R_x]).flatten()
    #     print "\nAfter flatten, y0 reads like \n", self.y0

    #     self.y0 = self.T_x

    #     try: 
    #         self.T_xt
    #     # basically this command is being run, not the bottom one
    #     except AttributeError:
    #         self.odeint_output = odeint(
    #             self.get_dTx_dt, y0=self.y0, t=self.t_array,
    #             full_output=1 
    #             )
    #         self.T_xt = self.odeint_output[0]

    #     # doesnt really get here
    #     else:
    #         self.y0 = self.T_xt[-1,:]
    #         self.odeint_output = odeint(
    #             self.get_dTx_dt, y0=self.y0, t=self.t_array,
    #             full_output=1 
    #             )
    #         self.T_xt = np.concatenate((self.T_xt, self.odeint_output[0]))
            
    #     print "\nDid get through all the calculations without error \n"
        
    
    # # def get_transient_error(self):
    # #     """ """
    # #     # find all the errors based on one run

    # #     # can also do something like following 
    # #     self.error = (
    # #         np.array([T_c_error, q_c_error, q_h_error, I_error]).flatten()
    # #         )
    # #     return self.error


    # #     return self.transient_error
















# working code 1 
#===========================================
#===========================================
#===========================================
# The following code is correct, if any screw up happens, copy and
# paste this in place of the codes given above. At least the following
# ones work. 
#===========================================
#===========================================

    # def get_dTx_dt(self, T, t):

    #     """Returns derivative of array of T wrt time.
    #     """
    #     # =============================================
    #     # make sure that a 2D array T is being returned
    #     # 2D array flattened into a very long 1D array
    #     # need to divide this T by number of varables
    #     # =============================================

    #     J = self.I/self.area
    #     dT_dt = np.zeros(T.size)
    #     dR_dx = np.zeros(T.size)
    #     dVs_dx = np.zeros(T.size)
    #     q0 = np.zeros(T.size)
    #     dq_dx_ss = np.zeros(T.size)
    #     dq_dx = np.zeros(T.size)
    #     dT_dx = np.zeros(T.size)

    #     dT_dx[1:-1] = 0.5 * (T[2:] - T[:-2]) / self.delta_x  
    #     dT_dx[0] = (T[1] - T[0]) / self.delta_x
    #     dT_dx[-1] = (T[-1] - T[-2]) / self.delta_x

    #     for i in range(self.nodes):

    #         T_props = T[i]  # i for central differencing
    #         self.set_TEproperties(T_props)
    #         self.set_ZT()

    #         dR_dx[i] = self.rho * self.area

    #         dVs_dx[i] = self.alpha * dT_dx[i]

    #         q0[i] = (
    #             J * T[i] * self.alpha - self.k * dT_dx[i]
    #             ) 

    #         dq_dx_ss[i] = (
    #             (self.rho * J ** 2. * (1. + self.ZT)) - J *
    #             self.alpha * q0[i] / self.k
    #             )

    #     # hot side BC, q_h
    #     q0[0] = self.U_hot * (self.T_h_conv - T[0]) 

    #     # cold side BC, q_c 
    #     q0[-1] = self.U_cold * (T[-1] - self.T_c_conv)

    #     dq_dx[1:-1] = (
    #         (q0[2:] - q0[:-2]) / (2. * self.delta_x)
    #         )
    #     dq_dx[0] = (
    #         (q0[1] - q0[0]) / self.delta_x
    #         )
    #     dq_dx[-1] = (
    #         (q0[-1] - q0[-2]) / self.delta_x
    #         )

    #     for i in range(self.nodes):

    #         T_props = T[i]  # i for central differencing
    #         self.set_TEproperties(T_props)
    #         self.set_ZT()

    #         dT_dt[i] = (
    #             1. / self.C * (-dq_dx[i] + dq_dx_ss[i])
    #             )

    #     # dVs_dx = self.alpha * dT_dx
    #     # dR_dx = self.rho * self.area
        
    #     # print "\ndT_dt is \n", dT_dt
    #     # print "\ndq_dx is \n", dq_dx
    #     # print "\ndVs_dx is \n", dVs_dx
    #     # print "\ndR_dx is \n", dR_dx

    #     # which dq_dx should I return??? both or just one?
    #     # need to return bunch of stuff from here
    #     # not only dT_dt

    #     # return self.dT_dt, self.dq_dx, self.dVs_dx, self.dR_dx
    #     # self.return_array = (
    #     #     np.array([dT_dt, dq_dx, dVs_dx, dR_dx]).flatten()
    #     #     )
    #     # print "return_array is ", self.return_array
    #     # return self.return_array
    #     # print "\ndT_dt is \n", dT_dt


    #     # self.return_array = (
    #     #     np.array([dT_dt, dq_dx, dVs_dx, dR_dx]).flatten()
    #     #     )
        
    #     # return self.return_array

    #     return dT_dt

    # def solve_leg_transient_once(self):

    #     """Solves leg based on array of transient BC's."""

    #     self.delta_x = self.x[1] - self.x[0]

    #     # self.y0 = np.array([self.T_x, self.q_x, self.Vs_x, self.R_x]).flatten()
    #     # print "\nAfter flatten, y0 reads like \n", self.y0
        
    #     self.y0 = self.T_x
    #     # need to use flatten here and provide more guess

    #     try: 
    #         self.T_xt

    #     except AttributeError:
    #         self.odeint_output = odeint(
    #             self.get_dTx_dt, y0=self.y0, t=self.t_array,
    #             full_output=1 
    #             )
    #         self.T_xt = self.odeint_output[0]

    #     else:
    #         self.y0 = self.T_xt[-1,:]
    #         self.odeint_output = odeint(
    #             self.get_dTx_dt, y0=self.y0, t=self.t_array,
    #             full_output=1 
    #             )
    #         self.T_xt = np.concatenate((self.T_xt, self.odeint_output[0]))


    # # def get_transient_error(self):
    # #     """ """
    # #     # find all the errors based on one run

    # #     # can also do something like following 
    # #     self.error = (
    # #         np.array([T_c_error, q_c_error, q_h_error, I_error]).flatten()
    # #         )
    # #     return self.error


    # #     return self.transient_error





#===========================================
#===========================================
# Upto here 
#===========================================
#===========================================













































# This is for steady state which I dont know works or nor, I do not
# really use it since I mostly run te_pair and not a single leg


    # def solve_leg_for_real(self):

    #     self.fsolve_output0 = fsolve(self.get_error_I, x0= self.I)
    #     print "Final I is ", self.I
        
    # def get_error_I(self, I):
        
    #     self.I = I
    #     self.solve_leg()
    #     print "I in this step is ", self.I

    #     # this is what I needs to be equal to 
    #     self.I_correct = self.Vs / (self.R_load + self.R_internal)

    #     # the difference between the correct J and calculated J
    #     self.I_error = self.I_correct - self.I

    #     print "Error in I is ", self.I_error
    #     return self.I_error

 
