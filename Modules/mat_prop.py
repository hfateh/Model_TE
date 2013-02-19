# =================================================
# THERE IS A MUCH CLEANER MAT_PROP.PY AT THE BOTTOM
# =================================================


"""Module containing set properties function."""


import numpy as np


def import_raw_property_data(self):

    """Imports and sets values for material properties as a function
    of temperature.  These values come from literature, and they may
    come from experiments or curve fitting in the future.

    """

    print "running import_raw_property_data"

    if self.material == "marlow p-type":
        # added on 10/03/2012
        # we measured Seebeck coefficient for n and p-type
        # sigma and k are from literature right now
        # sigma and k are better performance than actual Marlow
        poly_deg = 3
        
        self.alpha_raw = np.array([[305.88834, 184.4201],
                                   [315.955, 180.9814],
                                   [325.96334, 183.9347],
                                   [335.91166, 187.8893],
                                   [346.045, 184.0316],
                                   [356.065, 185.297],
                                   [365.95001, 188.8359],
                                   [375.9667, 182.5712],
                                   [385.95, 187.3913],
                                   [396.0383, 185.5793],
                                   [405.9867, 179.1448],
                                   [415.9667, 178.0592],
                                   [426.0117, 180.6971],
                                   [435.995, 170.9204],
                                   [446.0167, 175.7425],
                                   [455.95, 165.7763],
                                   [465.9333, 165.8615],
                                   [476.0367, 158.0667],
                                   [485.955, 145.0191],
                                   [495.885, 144.997],
                                   [505.9884, 133.7025],
                                   [515.995, 125.7248],
                                   [525.9117, 122.7499],
                                   [535.98, 115.2979],
                                   [545.935, 105.3918],
                                   [555.855, 104.7148],
                                   [565.9833, 101.7378],
                                   [575.9683, 95.72422]])
        self.alpha_params = np.polyfit(
            self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
            )

        self.k_raw = np.array([[299.5228426396, 1.3873417722],
                          [321.8578680203, 1.3265822785],
                          [345.5888324873, 1.3164556962],
                          [369.3197969543, 1.3569620253],
                          [391.654822335, 1.3873417722],
                          [416.781725888, 1.4683544304],
                          [440.512690355, 1.6],
                          [462.847715736, 1.7620253165],
                          [474.015228426, 1.8329113924],
                          [497.746192893, 2.035443038],
                          [522.873096447, 2.2075949367]])
        self.k_params = np.polyfit(
            self.k_raw[:, 0], self.k_raw[:, 1], poly_deg
            )

        self.sigma_raw = np.array([[299.4305754926, 9.746835443],
                              [323.4029100874, 8.6835443038],
                              [344.5312214537, 7.5443037975],
                              [369.9479968681, 6.7088607595],
                              [392.5476641, 6.0253164557],
                              [416.58280047, 5.4936708861],
                              [440.626908521, 5.0379746835],
                              [463.271434164, 4.7341772152],
                              [473.158227848, 4.4303797468],
                              [497.229250946, 4.2025316456],
                              [522.744714864, 4.2025316456]])
        self.sigma_params = np.polyfit(
            self.sigma_raw[:, 0], self.sigma_raw[:, 1], poly_deg
            )

    if self.material == "marlow n-type":
        # added on 10/03/2012
        # we measured Seebeck coefficient for n and p-type
        # sigma and k are from literature right now
        # sigma and k are better performance than actual Marlow
        poly_deg = 3

        self.alpha_raw = np.array([[305.88834, -184.4201],
                                   [316.02167, -174.2899],
                                   [325.995, -169.4341],
                                   [335.97834, -168.759],
                                   [346.02834, -170.3785],
                                   [355.98333, -172.1414],
                                   [365.98334, -174.0112],
                                   [375.9, -169.3056],
                                   [385.9517, -171.4511],
                                   [395.9383, -165.3796],
                                   [405.905, -163.3501],
                                   [415.95, -159.9725],
                                   [425.9117, -161.7775],
                                   [435.9783, -151.9636],
                                   [445.9483, -156.6737],
                                   [455.9667, -146.8893],
                                   [465.905, -142.9715],
                                   [475.9833, -135.9831],
                                   [485.8583, -133.8981],
                                   [496.005, -130.982],
                                   [505.8667, -127.4579],
                                   [515.9284, -121.2369],
                                   [525.9783, -118.1372],
                                   [535.915, -109.9285],
                                   [545.955, -105.78],
                                   [555.9167, -105.4163],
                                   [565.905, -100.953],
                                   [575.9317, -96.50673],
                                   [575.9384, -97.34524]])
        self.alpha_params = np.polyfit(
            self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
            )

        self.k_raw = np.array([[299.5228426396, 1.3873417722],
                               [321.8578680203, 1.3265822785],
                               [345.5888324873, 1.3164556962],
                               [369.3197969543, 1.3569620253],
                               [391.654822335, 1.3873417722],
                               [416.781725888, 1.4683544304],
                               [440.512690355, 1.6],
                               [462.847715736, 1.7620253165],
                               [474.015228426, 1.8329113924],
                               [497.746192893, 2.035443038],
                               [522.873096447, 2.2075949367]])
        self.k_params = np.polyfit(
            self.k_raw[:, 0], self.k_raw[:, 1], poly_deg
            )

        self.sigma_raw = np.array([[299.4305754926, 9.746835443],
                                   [323.4029100874, 8.6835443038],
                                   [344.5312214537, 7.5443037975],
                                   [369.9479968681, 6.7088607595],
                                   [392.5476641, 6.0253164557],
                                   [416.58280047, 5.4936708861],
                                   [440.626908521, 5.0379746835],
                                   [463.271434164, 4.7341772152],
                                   [473.158227848, 4.4303797468],
                                   [497.229250946, 4.2025316456],
                                   [522.744714864, 4.2025316456]])
        self.sigma_params = np.polyfit(
            self.sigma_raw[:, 0], self.sigma_raw[:, 1], poly_deg
            )

    if self.material == "BiTe variable n-type":
        poly_deg = 3

        self.alpha_raw = np.array([[112.301006188, -43.0457272427],
                              [137.129765223, -53.5256743351],
                              [156.431574472, -62.2549482734],
                              [179.857338561, -73.5981860267],
                              [210.184704048, -87.5634182059],
                              [239.165761255, -98.9227171808],
                              [272.270773301, -112.895979971],
                              [308.195946903, -124.275355473],
                              [339.954650668, -135.642685058],
                              [373.158864377, -143.544806084]])
        self.alpha_params = np.polyfit(
            self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
            )

        self.k_raw = np.array([[100, 29.467689848],
                          [130, 27.8005198844],
                          [160, 26.4696917156],
                          [190, 24.9615826351],
                          [220, 23.1842086837],
                          [250, 21.0620306825],
                          [280, 18.4187680451],
                          [310, 14.8860234335],
                          [340, 9.7526188408],
                          [370, 1.6249327924]])
        self.k_params = np.polyfit(self.k_raw[:, 0], self.k_raw[:, 1], poly_deg)

        self.sigma_raw = np.array([[109.722222222, 64.802436126],
                              [125, 59.3315508021],
                              [138.888888889, 52.8995840761],
                              [158.333333333, 46.7825311943],
                              [179.166666667, 41.3057040998],
                              [201.388888889, 36.1482471777],
                              [229.166666667, 30.3431372549],
                              [259.722222222, 25.4976232917],
                              [284.722222222, 22.2623291741],
                              [309.722222222, 19.3478906714],
                              [336.111111111, 16.7528223411],
                              [359.722222222, 15.1232917409],
                              [376.388888889, 13.8220439691]])
        self.sigma_params = np.polyfit(
            self.sigma_raw[:, 0], self.sigma_raw[:, 1], poly_deg
            )

    if self.material == "BiTe variable p-type":
        # properties trial 2 - need to add a comment
        poly_deg = 3

        self.alpha_raw = np.array([[112.301006188, 43.0457272427],
                              [137.129765223, 53.5256743351],
                              [156.431574472, 62.2549482734],
                              [179.857338561, 73.5981860267],
                              [210.184704048, 87.5634182059],
                              [239.165761255, 98.9227171808],
                              [272.270773301, 112.895979971],
                              [308.195946903, 124.275355473],
                              [339.954650668, 135.642685058],
                              [373.158864377, 143.544806084]])
        self.alpha_params = np.polyfit(
            self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
            )

        self.k_raw = np.array([[100, 29.467689848],
                          [130, 27.8005198844],
                          [160, 26.4696917156],
                          [190, 24.9615826351],
                          [220, 23.1842086837],
                          [250, 21.0620306825],
                          [280, 18.4187680451],
                          [310, 14.8860234335],
                          [340, 9.7526188408],
                          [370, 1.6249327924]])
        self.k_params = np.polyfit(self.k_raw[:, 0], self.k_raw[:, 1], poly_deg)

        self.sigma_raw = np.array([[109.722222222, 64.802436126],
                              [125, 59.3315508021],
                              [138.888888889, 52.8995840761],
                              [158.333333333, 46.7825311943],
                              [179.166666667, 41.3057040998],
                              [201.388888889, 36.1482471777],
                              [229.166666667, 30.3431372549],
                              [259.722222222, 25.4976232917],
                              [284.722222222, 22.2623291741],
                              [309.722222222, 19.3478906714],
                              [336.111111111, 16.7528223411],
                              [359.722222222, 15.1232917409],
                              [376.388888889, 13.8220439691]])
        self.sigma_params = np.polyfit(
            self.sigma_raw[:, 0], self.sigma_raw[:, 1], poly_deg
            )

    if self.material == "typical BiTe n-type":
        # Extracted from Bed Poudel et al, Science 320, 634 (2008)
        # This was the properties used for first trial of validation
        # process.
        poly_deg = 3

        self.alpha_raw = np.array([[297.3450032873, -213.717948718],
                              [321.1747205786, -223.974358974],
                              [343.6845825115, -227.820512821],
                              [367.6137409599, -231.025641026],
                              [391.61522025, -229.102564103],
                              [414.22452334, -225.897435897],
                              [439.726660092, -217.564102564],
                              [462.462524655, -205.384615385],
                              [472.47896121, -195.128205128],
                              [496.72452334, -175.897435897],
                              [522.425542406, -153.46153846]])
        self.alpha_params = np.polyfit(
            self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
            )

        self.k_raw = np.array([[299.5228426396, 1.3873417722],
                          [321.8578680203, 1.3265822785],
                          [345.5888324873, 1.3164556962],
                          [369.3197969543, 1.3569620253],
                          [391.654822335, 1.3873417722],
                          [416.781725888, 1.4683544304],
                          [440.512690355, 1.6],
                          [462.847715736, 1.7620253165],
                          [474.015228426, 1.8329113924],
                          [497.746192893, 2.035443038],
                          [522.873096447, 2.2075949367]])
        self.k_params = np.polyfit(
            self.k_raw[:, 0], self.k_raw[:, 1], poly_deg
            )

        self.sigma_raw = np.array([[299.4305754926, 9.746835443],
                              [323.4029100874, 8.6835443038],
                              [344.5312214537, 7.5443037975],
                              [369.9479968681, 6.7088607595],
                              [392.5476641, 6.0253164557],
                              [416.58280047, 5.4936708861],
                              [440.626908521, 5.0379746835],
                              [463.271434164, 4.7341772152],
                              [473.158227848, 4.4303797468],
                              [497.229250946, 4.2025316456],
                              [522.744714864, 4.2025316456]])
        self.sigma_params = np.polyfit(
            self.sigma_raw[:, 0], self.sigma_raw[:, 1], poly_deg
            )

    if self.material == "typical BiTe p-type":

        poly_deg = 3
        # Extracted from Bed Poudel et al, Science 320, 634 (2008)
        # This was the properties used for first trial of validation
        # process.
        self.alpha_raw = np.array([[297.3450032873, 213.717948718],
                              [321.1747205786, 223.974358974],
                              [343.6845825115, 227.820512821],
                              [367.6137409599, 231.025641026],
                              [391.61522025, 229.102564103],
                              [414.22452334, 225.897435897],
                              [439.726660092, 217.564102564],
                              [462.462524655, 205.384615385],
                              [472.47896121, 195.128205128],
                              [496.72452334, 175.897435897],
                              [522.425542406, 153.46153846]])
        self.alpha_params = np.polyfit(
            self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
            )

        self.k_raw = np.array([[299.5228426396, 1.3873417722],
                          [321.8578680203, 1.3265822785],
                          [345.5888324873, 1.3164556962],
                          [369.3197969543, 1.3569620253],
                          [391.654822335, 1.3873417722],
                          [416.781725888, 1.4683544304],
                          [440.512690355, 1.6],
                          [462.847715736, 1.7620253165],
                          [474.015228426, 1.8329113924],
                          [497.746192893, 2.035443038],
                          [522.873096447, 2.2075949367]])
        self.k_params = np.polyfit(
            self.k_raw[:, 0], self.k_raw[:, 1], poly_deg
            )

        self.sigma_raw = np.array([[299.4305754926, 9.746835443],
                              [323.4029100874, 8.6835443038],
                              [344.5312214537, 7.5443037975],
                              [369.9479968681, 6.7088607595],
                              [392.5476641, 6.0253164557],
                              [416.58280047, 5.4936708861],
                              [440.626908521, 5.0379746835],
                              [463.271434164, 4.7341772152],
                              [473.158227848, 4.4303797468],
                              [497.229250946, 4.2025316456],
                              [522.744714864, 4.2025316456]])
        self.sigma_params = np.polyfit(
            self.sigma_raw[:, 0], self.sigma_raw[:, 1], poly_deg
            )

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

    if self.material == "MgSi":
        # Raw data comes from Gao et al. MgSi is n-type
        poly_deg = 2
        # print "Curve fitting for MgSi"

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


def set_properties_v_temp(self, T_props):

    """ Sets properties based on polynomial fit values.

    Used by set_TEproperties to set the temperature-dependent
    properties of materials for which polynomial curve fits have been
    done.

    This may need to changed to spline fitting at some point for a
    more accurate fit.

    Inputs:

    T_props : temperature (K) at which properties are to be evaluated

    """
    try:
        self.alpha_params
    except AttributeError:
        self.import_raw_property_data()
    self.alpha = (np.polyval(self.alpha_params, T_props) * 1.e-6)
    # Seebeck coefficient (V/K)
    self.k = np.polyval(self.k_params, T_props)
    # thermal conductivity (W/m-K)
    self.sigma = np.polyval(self.sigma_params, T_props) * 1.e4
    # electrical conductivity (S/m)
    self.rho = 1. / self.sigma
    # electrical resistivity (Ohm-m)


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

    elif self.material == 'BiTe variable n-type':
        self.set_properties_v_temp(T_props)

    elif self.material == 'BiTe variable p-type':
        self.set_properties_v_temp(T_props)

    elif self.material == 'typical BiTe n-type':
        self.set_properties_v_temp(T_props)

    elif self.material == 'typical BiTe p-type':
        self.set_properties_v_temp(T_props)

    elif self.material == 'marlow p-type':
        self.set_properties_v_temp(T_props)

    elif self.material == 'marlow n-type':
        self.set_properties_v_temp(T_props)

    # Material properties for validation trial. These properties are
    # for typical BiTe materials at 423 K. The properties for a range
    # of temperature are given above as a poly curve.
    elif self.material == 'BiTe variable n-type':
        self.k = 1.54  # Thermal conductivity (W/m-K)
        self.alpha = -150.e-6  # Seebeck coefficient (V/K)
        # I made this negative even though it's for a p-type
        # material.  This is just for a hypothetical model.
        self.rho = 9.0909 * 1.e-6
        # electrical resistivity (Ohm-m)

    elif self.material == 'BiTe variable p-type':
        self.k = 1.54  # Thermal conductivity (W/m-K)
        self.alpha = 150.e-6  # Seebeck coefficient (V/K)
        self.rho = 9.0909 * 1.e-6
        # electrical resistivity (Ohm-m)

    # Materials with properties that are either constant or dependent
    # on temperature by some convenient function.
    # from CRC TE Handbook Table 12.1
    elif self.material == 'ex1 n-type':
        self.k = 54. / T_props * 100.
        # thermal conductivity (W/m-K)
        self.alpha = (0.268 * T_props - 329.) * 1.e-6
        # Seebeck coefficient (V/K)
        self.sigma = (T_props - 310.) / 0.1746
        # electrical conductivity (S/cm)
        self.rho = 1. / self.sigma / 100.
        # electrical resistivity (Ohm-m)

    # from CRC TE Handbook Table 12.1
    elif self.material == 'ex1 p-type':
        self.k = 3.194 / T_props * 100.
        # thermal conductivity (W/m-K)
        self.alpha = (0.150 * T_props + 211.) * 1.e-6
        # Seebeck coefficient (V/K)
        self.sigma = 25.
        # electrical conductivity (S/cm)
        self.rho = 1. / self.sigma / 100.
        # electrical resistivity (Ohm-m)

    # from CRC TE Handbook Table 12.1
    elif self.material == 'ex2 n-type':
        self.k = 3. / T_props * 100.
        # thermal conductivity (W/m-K)
        self.alpha = (0.20 * T_props - 400.) * 1.e-6
        # Seebeck coefficient (V/K)
        self.sigma = 1.e5 / T_props
        # electrical conductivity (S/cm)
        self.rho = 1. / self.sigma / 100.
        # electrical resistivity (Ohm-m)

    # from CRC TE Handbook Table 12.1
    elif self.material == 'ex2 p-type':
        self.k = 10. / T_props * 100.
        # thermal conductivity (W/m-K)
        self.alpha = (200.) * 1.e-6
        # Seebeck coefficient (V/K)
        self.sigma = T_props
        # electrical conductivity (S/cm)
        self.rho = 1. / self.sigma / 100.
        # electrical resistivity (Ohm-m)

    # from CRC TE Handbook Table 12.1
    elif self.material == 'ex3 n-type':
        self.k = 3. / T_props * 100.
        # thermal conductivity (W/m-K)
        self.alpha = 0.20 * T_props * 1.e-6
        # Seebeck coefficient (V/K)
        self.sigma = 1000.
        # electrical conductivity (S/cm)
        self.rho = 1. / self.sigma / 100.
        # electrical resistivity (Ohm-m)

    # from CRC TE Handbook Table 12.1
    elif self.material == 'ex3 p-type':
        self.k = 10. / T_props * 100.
        # thermal conductivity (W/m-K)
        self.alpha = 200. * 1.e-6
        # Seebeck coefficient (V/K)
        self.sigma = T_props
        # electrical conductivity (S/cm)
        self.rho = 1. / self.sigma / 100.
        # electrical resistivity (Ohm-m)

    # From CRC TE Handbook Table 27.7
    elif self.material == 'constant BiTe n-type':
        self.k = 1.5  # Thermal conductivity (W/m-K)
        self.alpha = -206.e-6  # Seebeck coefficient (V/K)
        # I made this negative even though it's for a p-type
        # material.  This is just for a hypothetical model.
        self.rho = 8.89 * 1.e-6
        # electrical resistivity (Ohm-m)

    # From CRC TE Handbook Table 27.7
    elif self.material == 'constant BiTe p-type':
        self.k = 1.5  # Thermal conductivity (W/m-K)
        self.alpha = 206.e-6  # Seebeck coefficient (V/K)
        self.rho = 8.89 * 1.e-6
        # electrical resistivity (Ohm-m)

    # Alumina with pure conduction
    elif self.material == 'alumina':
        self.k = 1.5  # Thermal conductivity (W/m-K) this needs to be
                     # updated to what it actually is
        self.alpha = 1.e-9  # Seebeck coefficient (V/K)
        # I made this negative even though it's for a p-type
        # material.  This is just for a hypothetical model.
        self.rho = 1.
        # dummy value since I don't know what it
        # actually is and it doesn't matter.
        # electrical resistivity (Ohm-m)

    # Direct contact between hot and cold side
    elif self.material == 'none':
        self.k = 1.e9  # really high thermal condutivity (W/m-K) so
                      # that resistance is zero
        self.alpha = 1.e-9  # dummy value
        self.rho = 1.
        # dummy value
























































# # Haiyan Fateh
# # Updated on 18 Feb, 2013
# """Module containing set properties function."""


# import numpy as np


# def import_raw_property_data(self):

#     """Imports and sets values for material properties as a function
#     of temperature.  These values come from literature, and they may
#     come from experiments or curve fitting in the future.

#     """

#     print "running import_raw_property_data"

#     if self.material == "HMS":
#         # Raw data taken from Luo et al. HMS is p-type
#         poly_deg = 3
#         # print "Curve fitting for HMS"


#         self.alpha_raw = np.array([[315.705,232.42625],
#                             [330.76667, 231.47188],
#                             [345.84167, 230.17526],
#                             [360.81667, 235.14671],
#                             [375.785, 233.48876],
#                             [390.8333, 234.53112],
#                             [405.84, 231.80171],
#                             [420.8833, 231.16316],
#                             [435.8567, 232.57991],
#                             [450.85, 238.2938],
#                             [465.8633, 238.88064],
#                             [480.8667, 228.6744],
#                             [495.8534, 227.07233],
#                             [510.85, 228.60436],
#                             [525.8617, 228.28546],
#                             [540.85, 223.78528],
#                             [555.885, 219.25532],
#                             [570.895, 218.75806],
#                             [585.8233, 211.26455],
#                             [600.8817, 208.26649],
#                             [615.875, 207.52015],
#                             [630.85, 202.52498],
#                             [645.8533, 194.35518],
#                             [660.895, 184.67011],
#                             [675.8417, 185.92441],
#                             [690.8834, 186.66846],
#                             [705.88, 175.34614],
#                             [720.9167, 171.48397],
#                             [735.9333, 162.04891],
#                             [750.8651, 153.29339],
#                             [765.8734, 152.25655],
#                             [780.8667, 143.86649],
#                             [795.8667, 139.59344],
#                             [810.9, 131.00214],
#                             [825.905, 128.76147],
#                             [831.0483, 117.81014]])
#         self.alpha_params = np.polyfit(
#             self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
#             )

#         self.k_raw = np.array([[322.6, 2.5579899317],
#                                [372.9, 2.51810604],
#                                [423.3, 2.5371132561],
#                                [473.3, 2.5417150908],
#                                [523.3, 2.578562122],
#                                [573.3, 2.5915271088],
#                                [623.3, 2.5654607341],
#                                [673.3, 2.5650773755],
#                                [723.3, 2.5735915738],
#                                [773.3, 2.6374224558],
#                                [823.3, 2.7425225851],
#                                [873.3, 2.9102448442]])
#         self.k_params = np.polyfit(
#             self.k_raw[:, 0], self.k_raw[:, 1], poly_deg
#             )

#         self.sigma_raw = np.array([[315.705, 447.48359],
#                                    [330.76667, 427.31343],
#                                    [345.84167, 409.70235],
#                                    [360.81667, 392.78365],
#                                    [375.785, 378.55348],
#                                    [390.8333, 365.56717],
#                                    [405.84, 354.14173],
#                                    [420.8833, 343.39514],
#                                    [435.8567, 335.56024],
#                                    [450.85, 324.76268],
#                                    [465.8633, 316.43513],
#                                    [480.8667, 306.77187],
#                                    [495.8534, 301.23106],
#                                    [510.85, 294.17404],
#                                    [525.8617, 290.47311],
#                                    [540.85, 285.56701],
#                                    [555.885, 280.54875],
#                                    [570.895, 274.45275],
#                                    [585.8233, 271.61472],
#                                    [600.8817, 268.3056],
#                                    [615.875, 271.1431],
#                                    [630.85, 267.11241],
#                                    [645.8533, 265.20722],
#                                    [660.895, 264.45594],
#                                    [675.8417, 264.3013],
#                                    [690.8834, 270.92926],
#                                    [705.88, 263.90337],
#                                    [720.9167, 262.34766],
#                                    [735.9333, 264.72492],
#                                    [750.8651, 264.79832],
#                                    [765.8734, 278.1682],
#                                    [780.8667, 274.54929],
#                                    [795.8667, 282.64445],
#                                    [810.9, 284.56422],
#                                    [825.905, 287.06549],
#                                    [831.0483, 293.59086]])
#         self.sigma_params = np.polyfit(
#             self.sigma_raw[:, 0], self.sigma_raw[:, 1], poly_deg
#             )

#         self.density_raw = np.array([[322.6, 4786.],
#                                      [372.9, 4786.],
#                                      [423.3, 4786.],
#                                      [473.3, 4786.],
#                                      [523.3, 4786.],
#                                      [573.3, 4786.],
#                                      [623.3, 4786.],
#                                      [673.3, 4786.],
#                                      [723.3, 4786.],
#                                      [773.3, 4786.],
#                                      [823.3, 4786.],
#                                      [873.3, 4786.]])
#         self.density_params = np.polyfit(
#             self.density_raw[:, 0], self.density_raw[:, 1], poly_deg
#             )


#         self.c_raw = np.array([[322.6, 582.85],
#                                [372.9,	592.5],
#                                [423.3,	604.46],
#                                [473.3,	615.38],
#                                [523.3,	631.62],
#                                [573.3,	644.62],
#                                [623.3,	652.11],
#                                [673.3,	667.44],
#                                [723.3,	678.1],
#                                [773.3,	689.7],
#                                [823.3,	696.27],
#                                [873.3,	703.79]])
#         self.c_params = np.polyfit(
#             self.c_raw[:, 0], self.c_raw[:, 1], poly_deg
#             )

#     if self.material == "MgSi":
#         # Raw data comes from Gao et al. MgSi is n-type
#         poly_deg = 2

#         self.alpha_raw = np.array([[304.97667, -76.77944],
#                                    [320.13334, -78.54224],
#                                    [335.07833, -90.72432],
#                                    [350.11667, -90.27498],
#                                    [365.15, -95.12417],
#                                    [380.1783, -94.27527],
#                                    [395.085, -102.2901],
#                                    [410.1, -99.26342],
#                                    [425.1283, -106.00969],
#                                    [440.1167, -105.58557],
#                                    [455.1167, -107.40208],
#                                    [470.1333, -110.31651],
#                                    [485.185, -114.40619],
#                                    [500.1, -123.45917],
#                                    [515.095, -121.51015],
#                                    [530.0334, -126.49965],
#                                    [545.0667, -132.28749],
#                                    [560.1333, -128.94191],
#                                    [575.1034, -133.0895],
#                                    [590.1, -137.9517],
#                                    [605.1933, -143.75023],
#                                    [620.15, -147.52096],
#                                    [635.1684, -157.04423],
#                                    [650.0533, -149.56359],
#                                    [665.0683, -150.01531],
#                                    [680.1833, -158.23988],
#                                    [695.1117, -174.23052],
#                                    [710.0667, -175.58246],
#                                    [725.0367, -170.67081],
#                                    [740.1167, -167.63654],
#                                    [755.12, -186.00031],
#                                    [770.0717, -185.53065],
#                                    [785.0617, -181.88676],
#                                    [800.05, -197.40704],
#                                    [815.1167, -200.92754]])
#         self.alpha_params = np.polyfit(
#             self.alpha_raw[:, 0], self.alpha_raw[:, 1], poly_deg
#             )



#         self.k_raw = np.array([[323.15, 3.13731],
#                                [373.15, 3.00759],
#                                [423.15, 2.91603],
#                                [473.15, 2.84598],
#                                [523.15, 2.83228],
#                                [573.15, 2.77579],
#                                [623.15, 2.75821],
#                                [673.15, 2.6791],
#                                [723.15, 2.61775],
#                                [773.15, 2.5552],
#                                [803.15, 2.51657]])
#         self.k_params = np.polyfit(
#             self.k_raw[:, 0], self.k_raw[:, 1], poly_deg
#             )

#         self.sigma_raw = np.array([[304.97667,2025.4303064],
#                                    [320.13334,1953.0404522],
#                                    [335.07833,1929.119928],
#                                    [350.11667,1893.2014767],
#                                    [365.15,1813.2854751],
#                                    [380.1783,1849.1840193],
#                                    [395.085,1674.9071539],
#                                    [410.1,1731.5919452],
#                                    [425.1283,1657.625887],
#                                    [440.1167,1596.9892047],
#                                    [455.1167,1575.9967393],
#                                    [470.1333,1476.6277729],
#                                    [485.185,1417.9797169],
#                                    [500.1,1398.2341828],
#                                    [515.095,1374.0772588],
#                                    [530.0334,1332.5900323],
#                                    [545.0667,1319.6739576],
#                                    [560.1333,1285.4625348],
#                                    [575.1034,1288.3484849],
#                                    [590.1,1207.2220289],
#                                    [605.1933,1194.4769917],
#                                    [620.15,1132.0255237],
#                                    [635.1684,1135.2257844],
#                                    [650.0533,1223.7744333],
#                                    [665.0683,1109.2040017],
#                                    [680.1833,1047.5815571],
#                                    [695.1117,1056.8158787],
#                                    [710.0667,1008.1394804],
#                                    [725.0367,936.9948972],
#                                    [740.1167,952.0983141],
#                                    [755.12,856.7463382],
#                                    [770.0717,904.9065428],
#                                    [785.0617,847.1685289],
#                                    [800.05,764.5969332],
#                                    [815.1167,743.0947173]])
#         self.sigma_params = np.polyfit(self.sigma_raw[:, 0], self.sigma_raw[:, 1],
#                               poly_deg)

#         self.density_raw = np.array([[283.888641142, 2900.],
#                                      [396.056571319, 2900.],
#                                      [573.510861948, 2900.],
#                                      [786.035548194, 2900.],
#                                      [856.520354208, 2900.],
#                                      [901.20405173, 2900.]])
#         self.density_params = np.polyfit(
#             self.density_raw[:, 0], self.density_raw[:, 1], poly_deg
#             )


#         self.c_raw = np.array([[323.15, 1905.],
#                                [373.15, 1823.],
#                                [423.15, 1752.],
#                                [473.15, 1698.],
#                                [523.15, 1644.],
#                                [573.15, 1587.],
#                                [623.15, 1526.],
#                                [673.15, 1468.],
#                                [723.15, 1423.],
#                                [773.15, 1389.],
#                                [803.15, 1368.]])
#         self.c_params = np.polyfit(
#             self.c_raw[:, 0], self.c_raw[:, 1], poly_deg
#             )

# def set_properties_v_temp(self, T_props):

#     """ Sets properties based on polynomial fit values.
#     """
#     try:
#         self.alpha_params
#     except AttributeError:
#         self.import_raw_property_data()

#     # Seebeck coefficient (V/K)
#     self.alpha = (np.polyval(self.alpha_params, T_props) * 1.e-6)

#     # thermal conductivity (W/m-K)
#     self.k = np.polyval(self.k_params, T_props)

#     # electrical conductivity (S/m)
#     self.sigma = np.polyval(self.sigma_params, T_props) * 1.e2

#     # electrical resistivity (Ohm-m)
#     self.rho = 1. / self.sigma

#     # density (kg/m3)
#     self.density = np.polyval(self.density_params, T_props)

#     # specific heat (J/kg-K)
#     self.c = np.polyval(self.c_params, T_props)

#     # effective thermal mass, density times c (J/m3-K)
#     self.C = self.density * self.c

# def set_TEproperties(self, T_props):

#     """Sets TE properties

#     Inputs:

#     T_props : temperature (K) at which properties are to be evaluated

#     This method exists to separater materials with constant properties
#     from materials with temperature dependent properties.  It uses
#     set_properties_v_temp for the latter type of materials.  

#     """

#     self.T_props = T_props

#     # Materials with tabulated properties
#     if self.material == 'HMS':
#         self.set_properties_v_temp(T_props)

#     elif self.material == 'MgSi':
#         self.set_properties_v_temp(T_props)
