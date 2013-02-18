# Haiyan Fateh
# Updated on 18 Feb, 2013
"""Module containing set properties function."""


import numpy as np


def import_raw_property_data(self):

    """Imports and sets values for material properties as a function
    of temperature.  These values come from literature, and they may
    come from experiments or curve fitting in the future.

    """

    print "running import_raw_property_data"

    if self.material == "HMS":
        # Raw data taken from Luo et al. HMS is p-type
        poly_deg = 3
        # print "Curve fitting for HMS"


        self.alpha_raw = np.array([[315.705,232.42625],
                            [330.76667, 231.47188],
                            [345.84167, 230.17526],
                            [360.81667, 235.14671],
                            [375.785, 233.48876],
                            [390.8333, 234.53112],
                            [405.84, 231.80171],
                            [420.8833, 231.16316],
                            [435.8567, 232.57991],
                            [450.85, 238.2938],
                            [465.8633, 238.88064],
                            [480.8667, 228.6744],
                            [495.8534, 227.07233],
                            [510.85, 228.60436],
                            [525.8617, 228.28546],
                            [540.85, 223.78528],
                            [555.885, 219.25532],
                            [570.895, 218.75806],
                            [585.8233, 211.26455],
                            [600.8817, 208.26649],
                            [615.875, 207.52015],
                            [630.85, 202.52498],
                            [645.8533, 194.35518],
                            [660.895, 184.67011],
                            [675.8417, 185.92441],
                            [690.8834, 186.66846],
                            [705.88, 175.34614],
                            [720.9167, 171.48397],
                            [735.9333, 162.04891],
                            [750.8651, 153.29339],
                            [765.8734, 152.25655],
                            [780.8667, 143.86649],
                            [795.8667, 139.59344],
                            [810.9, 131.00214],
                            [825.905, 128.76147],
                            [831.0483, 117.81014]])
        self.alpha_params = np.polyfit(
            self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
            )

        self.k_raw = np.array([[322.6, 2.5579899317],
                               [372.9, 2.51810604],
                               [423.3, 2.5371132561],
                               [473.3, 2.5417150908],
                               [523.3, 2.578562122],
                               [573.3, 2.5915271088],
                               [623.3, 2.5654607341],
                               [673.3, 2.5650773755],
                               [723.3, 2.5735915738],
                               [773.3, 2.6374224558],
                               [823.3, 2.7425225851],
                               [873.3, 2.9102448442]])
        self.k_params = np.polyfit(
            self.k_raw[:, 0], self.k_raw[:, 1], poly_deg
            )

        self.sigma_raw = np.array([[315.705, 447.48359],
                                   [330.76667, 427.31343],
                                   [345.84167, 409.70235],
                                   [360.81667, 392.78365],
                                   [375.785, 378.55348],
                                   [390.8333, 365.56717],
                                   [405.84, 354.14173],
                                   [420.8833, 343.39514],
                                   [435.8567, 335.56024],
                                   [450.85, 324.76268],
                                   [465.8633, 316.43513],
                                   [480.8667, 306.77187],
                                   [495.8534, 301.23106],
                                   [510.85, 294.17404],
                                   [525.8617, 290.47311],
                                   [540.85, 285.56701],
                                   [555.885, 280.54875],
                                   [570.895, 274.45275],
                                   [585.8233, 271.61472],
                                   [600.8817, 268.3056],
                                   [615.875, 271.1431],
                                   [630.85, 267.11241],
                                   [645.8533, 265.20722],
                                   [660.895, 264.45594],
                                   [675.8417, 264.3013],
                                   [690.8834, 270.92926],
                                   [705.88, 263.90337],
                                   [720.9167, 262.34766],
                                   [735.9333, 264.72492],
                                   [750.8651, 264.79832],
                                   [765.8734, 278.1682],
                                   [780.8667, 274.54929],
                                   [795.8667, 282.64445],
                                   [810.9, 284.56422],
                                   [825.905, 287.06549],
                                   [831.0483, 293.59086]])
        self.sigma_params = np.polyfit(
            self.sigma_raw[:, 0], self.sigma_raw[:, 1], poly_deg
            )

        self.density_raw = np.array([[322.6, 4786.],
                                     [372.9, 4786.],
                                     [423.3, 4786.],
                                     [473.3, 4786.],
                                     [523.3, 4786.],
                                     [573.3, 4786.],
                                     [623.3, 4786.],
                                     [673.3, 4786.],
                                     [723.3, 4786.],
                                     [773.3, 4786.],
                                     [823.3, 4786.],
                                     [873.3, 4786.]])
        self.density_params = np.polyfit(
            self.density_raw[:, 0], self.density_raw[:, 1], poly_deg
            )


        self.c_raw = np.array([[322.6, 582.85],
                               [372.9,	592.5],
                               [423.3,	604.46],
                               [473.3,	615.38],
                               [523.3,	631.62],
                               [573.3,	644.62],
                               [623.3,	652.11],
                               [673.3,	667.44],
                               [723.3,	678.1],
                               [773.3,	689.7],
                               [823.3,	696.27],
                               [873.3,	703.79]])
        self.c_params = np.polyfit(
            self.c_raw[:, 0], self.c_raw[:, 1], poly_deg
            )

    if self.material == "MgSi":
        # Raw data comes from Gao et al. MgSi is n-type
        poly_deg = 2

        self.alpha_raw = np.array([[311.289993567, -111.872146119],
                              [464.006967001, -141.552511416],
                              [644.121200709, -184.931506849],
                              [777.984904831, -207.762557078]])
        self.alpha_params = np.polyfit(
            self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
            )

        self.k_raw = np.array([[291.236965464, 2.80871520138],
                          [472.020791479, 2.62097005644],
                          [725.982971396, 2.38897924041],
                          [576.615963519, 2.50282215632]])
        self.k_params = np.polyfit(
            self.k_raw[:, 0], self.k_raw[:, 1], poly_deg
            )

        self.sigma_raw = np.array([[307.385007162, 13.156135604],
                              [456.638548464, 9.79627566449],
                              [574.442145472, 8.21502466974],
                              [722.524271845, 7.17849753303]])
        self.sigma_params = np.polyfit(self.sigma_raw[:, 0], self.sigma_raw[:, 1],
                              poly_deg)

        self.density_raw = np.array([[283.888641142, 4786.],
                              [396.056571319, 4786.],
                              [573.510861948, 4786.],
                              [786.035548194, 4786.],
                              [856.520354208, 4786.],
                              [901.20405173, 4786.]])
        self.density_params = np.polyfit(
            self.density_raw[:, 0], self.density_raw[:, 1], poly_deg
            )

        self.c_raw = np.array([[322.6, 582.85],
                               [372.9,	592.5],
                               [423.3,	604.46],
                               [473.3,	615.38],
                               [523.3,	631.62],
                               [573.3,	644.62],
                               [623.3,	652.11],
                               [673.3,	667.44],
                               [723.3,	678.1],
                               [773.3,	689.7],
                               [823.3,	696.27],
                               [873.3,	703.79]])
        self.c_params = np.polyfit(
            self.c_raw[:, 0], self.c_raw[:, 1], poly_deg
            )

def set_properties_v_temp(self, T_props):

    """ Sets properties based on polynomial fit values.
    """
    try:
        self.alpha_params
    except AttributeError:
        self.import_raw_property_data()

    # Seebeck coefficient (V/K)
    self.alpha = (np.polyval(self.alpha_params, T_props) * 1.e-6)

    # thermal conductivity (W/m-K)
    self.k = np.polyval(self.k_params, T_props)

    # electrical conductivity (S/m)
    self.sigma = np.polyval(self.sigma_params, T_props) * 1.e2

    # electrical resistivity (Ohm-m)
    self.rho = 1. / self.sigma

    # density
    self.density = np.polyval(self.density_params, T_props)

    # specific heat
    self.c = np.polyval(self.c_params, T_props)

    # effective thermal mass, density times c (J/m3-K)
    self.C = self.density * self.c

def set_TEproperties(self, T_props):

    """Sets TE properties

    Inputs:

    T_props : temperature (K) at which properties are to be evaluated

    This method exists to separater materials with constant properties
    from materials with temperature dependent properties.  It uses
    set_properties_v_temp for the latter type of materials.  

    """

    self.T_props = T_props

    # Materials with tabulated properties
    if self.material == 'HMS':
        self.set_properties_v_temp(T_props)

    elif self.material == 'MgSi':
        self.set_properties_v_temp(T_props)
