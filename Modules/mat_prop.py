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

        self.alpha_raw = np.array([[296.89119171, 138.265544041],
                              [380.829015544, 140.620466321],
                              [561.139896373, 176.845854922],
                              [701.03626943, 206.270725389],
                              [806.735751295, 217.652849741],
                              [900., 205.769430052]])
        self.alpha_params = np.polyfit(
            self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
            )

        self.k_raw = np.array([[300, 2.40620446533],
                          [485.869565217, 2.20460634548],
                          [593.47826087, 2.1252173913],
                          [707.608695652, 2.07168037603],
                          [815.217391304, 2.09607520564],
                          [900.0, 2.12944770858]])
        self.k_params = np.polyfit(
            self.k_raw[:, 0], self.k_raw[:, 1], poly_deg
            )

        self.sigma_raw = np.array([[283.888641142, 6.55346563038],
                              [396.056571319, 6.22507485507],
                              [573.510861948, 4.86979996178],
                              [786.035548194, 3.5398961585],
                              [856.520354208, 3.34810791871],
                              [901.20405173, 3.34610116583]])
        self.sigma_params = np.polyfit(
            self.sigma_raw[:, 0], self.sigma_raw[:, 1], poly_deg
            )

        self.density_raw = np.array([[283.888641142, 4786.],
                              [396.056571319, 4786.],
                              [573.510861948, 4786.],
                              [786.035548194, 4786.],
                              [856.520354208, 4786.],
                              [901.20405173, 4786.]])
        self.density_params = np.polyfit(
            self.density_raw[:, 0], self.density_raw[:, 1], poly_deg
            )

        self.c_raw = np.array([[283.888641142,  582.85],
                              [396.056571319,  582.85],
                              [573.510861948,  582.85],
                              [786.035548194,  582.85],
                              [856.520354208,  582.85],
                              [901.20405173,  582.85]])
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

        self.c_raw = np.array([[283.888641142, 582.85],
                              [396.056571319,  582.85],
                              [573.510861948,  582.85],
                              [786.035548194,  582.85],
                              [856.520354208,  582.85],
                              [901.20405173,  582.85]])
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
    self.sigma = np.polyval(self.sigma_params, T_props) * 1.e4

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
