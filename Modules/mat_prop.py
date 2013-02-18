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

        self.alpha_raw = np.array([[304.97667, -76.77944],
                                   [320.13334, -78.54224],
                                   [335.07833, -90.72432],
                                   [350.11667, -90.27498],
                                   [365.15, -95.12417],
                                   [380.1783, -94.27527],
                                   [395.085, -102.2901],
                                   [410.1, -99.26342],
                                   [425.1283, -106.00969],
                                   [440.1167, -105.58557],
                                   [455.1167, -107.40208],
                                   [470.1333, -110.31651],
                                   [485.185, -114.40619],
                                   [500.1, -123.45917],
                                   [515.095, -121.51015],
                                   [530.0334, -126.49965],
                                   [545.0667, -132.28749],
                                   [560.1333, -128.94191],
                                   [575.1034, -133.0895],
                                   [590.1, -137.9517],
                                   [605.1933, -143.75023],
                                   [620.15, -147.52096],
                                   [635.1684, -157.04423],
                                   [650.0533, -149.56359],
                                   [665.0683, -150.01531],
                                   [680.1833, -158.23988],
                                   [695.1117, -174.23052],
                                   [710.0667, -175.58246],
                                   [725.0367, -170.67081],
                                   [740.1167, -167.63654],
                                   [755.12, -186.00031],
                                   [770.0717, -185.53065],
                                   [785.0617, -181.88676],
                                   [800.05, -197.40704],
                                   [815.1167, -200.92754]])
        self.alpha_params = np.polyfit(
            self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
            )



        self.k_raw = np.array([[323.15, 3.13731],
                               [373.15, 3.00759],
                               [423.15, 2.91603],
                               [473.15, 2.84598],
                               [523.15, 2.83228],
                               [573.15, 2.77579],
                               [623.15, 2.75821],
                               [673.15, 2.6791],
                               [723.15, 2.61775],
                               [773.15, 2.5552],
                               [803.15, 2.51657]])
        self.k_params = np.polyfit(
            self.k_raw[:, 0], self.k_raw[:, 1], poly_deg
            )

        self.sigma_raw = np.array([[304.97667,2025.4303064],
                                   [320.13334,1953.0404522],
                                   [335.07833,1929.119928],
                                   [350.11667,1893.2014767],
                                   [365.15,1813.2854751],
                                   [380.1783,1849.1840193],
                                   [395.085,1674.9071539],
                                   [410.1,1731.5919452],
                                   [425.1283,1657.625887],
                                   [440.1167,1596.9892047],
                                   [455.1167,1575.9967393],
                                   [470.1333,1476.6277729],
                                   [485.185,1417.9797169],
                                   [500.1,1398.2341828],
                                   [515.095,1374.0772588],
                                   [530.0334,1332.5900323],
                                   [545.0667,1319.6739576],
                                   [560.1333,1285.4625348],
                                   [575.1034,1288.3484849],
                                   [590.1,1207.2220289],
                                   [605.1933,1194.4769917],
                                   [620.15,1132.0255237],
                                   [635.1684,1135.2257844],
                                   [650.0533,1223.7744333],
                                   [665.0683,1109.2040017],
                                   [680.1833,1047.5815571],
                                   [695.1117,1056.8158787],
                                   [710.0667,1008.1394804],
                                   [725.0367,936.9948972],
                                   [740.1167,952.0983141],
                                   [755.12,856.7463382],
                                   [770.0717,904.9065428],
                                   [785.0617,847.1685289],
                                   [800.05,764.5969332],
                                   [815.1167,743.0947173]])
        self.sigma_params = np.polyfit(self.sigma_raw[:, 0], self.sigma_raw[:, 1],
                              poly_deg)

        self.density_raw = np.array([[283.888641142, 2900.],
                                     [396.056571319, 2900.],
                                     [573.510861948, 2900.],
                                     [786.035548194, 2900.],
                                     [856.520354208, 2900.],
                                     [901.20405173, 2900.]])
        self.density_params = np.polyfit(
            self.density_raw[:, 0], self.density_raw[:, 1], poly_deg
            )


        self.c_raw = np.array([[323.15, 1905.],
                               [373.15, 1823.],
                               [423.15, 1752.],
                               [473.15, 1698.],
                               [523.15, 1644.],
                               [573.15, 1587.],
                               [623.15, 1526.],
                               [673.15, 1468.],
                               [723.15, 1423.],
                               [773.15, 1389.],
                               [803.15, 1368.]])
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

    # density (kg/m3)
    self.density = np.polyval(self.density_params, T_props)

    # specific heat (J/kg-K)
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
