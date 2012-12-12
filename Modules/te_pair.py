import numpy as np
from scipy.optimize import fsolve
import leg
reload(leg)

class TE_Pair(object):
    """ TE_pair
    """

    def __init__(self):

        """Sets attributes and instantiates classes.
        """

        self.R_load = 4.*1.0/256.0
        self.leg_area_ratio = 0.7
        self.fill_fraction = 0.03
        self.length = 1.e-3
        self.area = (1.5e-3) ** 2
        self.Vs = 1.64/256. # initial guess for Voc
        self.R_internal = 1./256. # initial guess for R_internal

        self.Ptype = leg.Leg()
        self.Ntype = leg.Leg()
        self.Ptype.material = 'HMS'
        self.Ntype.material = 'MgSi'
        self.nodes = 10
        self.set_J()
        self.set_constants()

    def set_J(self):

        """Sets a single J value for a TE pair
        """
        self.J = self.Vs / (self.R_load + self.R_internal)
        print "\nGuess for J is", self.J, "\n"

    def set_constants(self):

        """Sets a bunch of attributes that are usually held constant.
        """
        self.Ntype.length = self.length
        self.Ptype.length = self.length
        self.Ptype.nodes = self.nodes
        self.Ntype.nodes = self.nodes
        self.Ntype.area = self.area
        self.Ptype.area = self.area
        self.Ntype.set_constants()
        self.Ptype.set_constants()
        
        self.Ntype.Vs = self.Vs
        self.Ptype.Vs = self.Vs
        self.Ntype.R_internal = self.R_internal
        self.Ptype.R_internal = self.R_internal
        self.Ntype.J = - self.J
        self.Ptype.J = self.J

    def set_q_guess(self):

        """Sets cold side guess for both Ntype and Ptype legs.
        """
        self.Ntype.set_q_guess()
        self.Ptype.set_q_guess()

    def set_TEproperties(self, T_props):

        """Sets properties for both legs based on temperature.
        """
        self.Ntype.set_TEproperties(T_props)
        self.Ptype.set_TEproperties(T_props)

    def solve_te_pair_once(self):

        """Solves legs and combines results of leg pair.
        """
        self.Ntype.solve_leg_once(self.Ntype.q_h)
        self.Ptype.solve_leg_once(self.Ptype.q_h)
        
        self.T_c = self.Ntype.T_c

        # area averaged hot side heat flux (kW/m^2)
        self.q_h = (
            (self.Ptype.q_h * self.Ptype.area + self.Ntype.q_h *
             self.Ntype.area) / self.area * 0.001
            )
        print "self.q_h is ", self.q_h

        # area averaged cold side heat flux (kW/m^2)
        self.q_c = (
            (self.Ptype.q_c * self.Ptype.area + self.Ntype.q_c *
             self.Ntype.area) / self.area * 0.001
            )

        self.h_eff= self.q_h / (self.T_h - self.T_c)
        self.R_thermal = 1. / self.h_eff

    def get_error(self, knob_arr):

        """Returns BC error.
        """
        self.Ntype.q_h = knob_arr[0]
        self.Ptype.q_h = knob_arr[1]
        self.T_h = knob_arr[2]

        self.Ptype.T_h = self.T_h
        self.Ntype.T_h = self.T_h

        self.solve_te_pair_once()

        self.q_c_conv = self.U_cold * (self.T_c - self.T_c_conv)
        self.q_h_conv = self.U_hot * (self.T_h_conv - self.T_h)

        T_c_error = self.Ntype.T_c - self.Ptype.T_c
        q_c_error = self.q_c - self.q_c_conv
        q_h_error = self.q_h - self.q_h_conv

        self.error = (
            np.array([T_c_error, q_c_error, q_h_error]).flatten()
            )
        return self.error

    def solve_te_pair(self):

        """Solves legs and combines results of leg pair.
        """

        self.Ptype.T_h = self.T_h_conv 
        self.Ntype.T_h = self.T_h_conv
        self.Ptype.T_c = self.T_c_conv
        self.Ntype.T_c = self.T_c_conv

        self.set_q_guess()
        knob_arr0 = (
            np.array([self.Ntype.q_h_guess, self.Ptype.q_h_guess,
        self.T_h_conv])
            )

        self.Ptype.T_c_goal = None
        self.Ntype.T_c_goal = None

        self.fsolve_output = fsolve(self.get_error, x0=knob_arr0)

        self.P = (self.Ntype.P + self.Ptype.P) * 0.001
        # power for the entire leg pair(kW). Negative sign makes this
        # a +ve number. Heat flux is negative so efficiency needs
        # a negative sign also.
        self.P_flux = self.P / self.area
        self.Vs = -self.Ntype.Vs + self.Ptype.Vs
        self.V = self.J * self.R_load / self.area
        self.R_internal = ( 
            self.Ntype.R_internal + self.Ptype.R_internal
            )

    # def get_J_error(self, J):
    #     """Return the error in actual and guessed J value
    #     """
    #     self.solve_te_pair()
        
    #     self.J_correct = (
    #         self.Vs / (self.R_load + self.Ntype.R_internal +
    #         self.Ptype.R_internal)
    #         )
        
    #     self.J_error = self.J_correct - self.J
    #     return self.J_error

    # def solve_te_pair_for_real(self):
    #     """Solves TE pair iterating through different J values
    #     """
    #     self.fsolve_output0 = fsolve(self.get_J_error, x0= self.J)

    # def set_ZT(self):

    #     """Sets ZT based on whatever properties were used last."""

    #     self.ZT = ( ((self.Ptype.alpha - self.Ntype.alpha) /
    #     ((self.Ptype.rho * self.Ptype.k) ** 0.5 + (self.Ntype.rho *
    #     self.Ntype.k) ** 0.5)) ** 2. * self.T_props )

    # def set_eta_max(self):

    #     """Sets theoretical maximum efficiency.

    #     Methods:

    #     self.set_TEproperties(T_props)

    #     Uses material properties evaluated at the average temperature
    #     based on Sherman's analysis.

    #     """

    #     self.T_props = 0.5 * (self.T_h + self.T_c)
    #     self.set_TEproperties(T_props=self.T_props)
    #     self.set_ZT()
    #     delta_T = self.T_h - self.T_c
    #     self.eta_max = ( delta_T / self.T_h * ((1. + self.ZT) ** 0.5 -
    #     1.) / ((1. + self.ZT) ** 0.5 + self.T_c / self.T_h) )

    # def set_A_opt(self):

    #     """Sets Ntype / Ptype area that results in max efficiency.

    #     Methods:

    #     self.set_TEproperties(T_props)

    #     Based on material properties evaluated at the average
    #     temperature.

    #     """

    #     self.set_TEproperties(T_props=self.T_props)
    #     self.A_opt = np.sqrt(self.Ntype.rho * self.Ptype.k /
    #     (self.Ptype.rho * self.Ntype.k))

    # def set_power_max(self):

    #     """Sets power factor and maximum theoretical power.

    #     Methods:

    #     self.Ntype.set_power_factor
    #     self.Ptype.set_power_factor

    #     """

    #     self.Ntype.set_power_factor()
    #     self.Ptype.set_power_factor()
    #     self.power_max = self.Ntype.power_max + self.Ptype.power_max

    # def set_leg_areas(self):

    #     """Sets leg areas and void area.

    #     Based on leg area ratio and fill fraction.

    #     self.Ptype.area must be held constant.  self.Ntype.area and
    #     self.area_void are varied here.

    #     """

    #     leg_area_ratio = self.leg_area_ratio
    #     fill_fraction = self.fill_fraction

    #     self.Ntype.area = self.Ptype.area * leg_area_ratio
    #     self.area_void = (
    #         (1. - fill_fraction) / fill_fraction * (self.Ptype.area +
    #     self.Ntype.area)
    #         )
    #     self.area = self.Ntype.area + self.Ptype.area + self.area_void
        
    # def get_minpar(self, apar):

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
    #         print "\noptimizaton iteration", self.opt_iter
    #         print "leg length =", self.length, "m"
    #         print "fill fraction =", self.fill_fraction * 100., "%"
    #         print "current =", self.I, "A"
    #         print "area ratio =", self.leg_area_ratio
    #         print "power flux (kW/m^2)", self.P_flux
    #     apar = np.array(apar)

    #     self.length = apar[0]
    #     self.fill_fraction = apar[1]
    #     self.I = apar[2]
    #     self.leg_area_ratio = apar[3]

    #     # reset surrogate variables
    #     self.set_constants()

    #     self.solve_te_pair()

    #     if (apar <= 0.).any():
    #         minpar = np.abs(self.P_flux) ** 3. + 100
    #         print "Encountered impossible value."

    #     else:
    #         minpar = - self.P_flux

    #     return minpar

    # def optimize(self):

    #     """Minimizes self.get_minpar

    #     Methods:

    #     self.get_minpar

    #     self.x0 and self.xb must be defined elsewhere."""

    #     time.clock()

    #     # dummy function that might be used with minimization
    #     def fprime():
    #         return 1

    #     self.opt_iter = 0

    #     self.x0 = np.array([self.length, self.fill_fraction,
    #     self.I, self.leg_area_ratio])

    #     from scipy.optimize import fmin

    #     self.xmin = fmin(self.get_minpar, self.x0)

    #     t1 = time.clock()

    #     print '\n'

    #     print "Optimized parameters:"
    #     print "leg length =", self.length, "m"
    #     print "fill fraction =", self.fill_fraction * 100., "%"
    #     print "current =", self.I, "A"
    #     print "area ratio =", self.leg_area_ratio

    #     print "\npower:", self.P * 1000., 'W'
    #     print "power flux:", self.P_flux, "kW/m^2"

    #     print """Elapsed time solving xmin1 =""", t1
