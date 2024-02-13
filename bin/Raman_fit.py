# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 12:30:38 2023

Library with all functions to fit Raman data once baseline is removed

@author:  J.Anaya


# Copyright (C) 2023 Julian anaya Calvo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

# Raman_fit
import numpy as np
from scipy.special import wofz
import lmfit
from lmfit import Parameters, Minimizer
from bin import Raman_plot



def f_0(x):
    return gauss_f(x, a_1=0, a_2=1, a_3=0)

def gauss_f(x, a_1=0, a_2=1, a_3=1):
    """
    Calculate the Gaussian shape.

    Parameters:
        x (array-like): The input variable.
        a_1 (float): The center of the Gaussian.
        a_2: The Full Width at Half Maximum (FWHM) of the Gaussian.
        a_3(float): The amplitude of the Gaussian.

    Returns:
        array-like: The Gaussian shape evaluated at x.
    """
    sigma = a_2 / (2 * np.sqrt(2 * np.log(2)))
    return a_3 * np.exp(-((x - a_1)**2) / (2 * sigma**2))


def lorentz_f(x, a_1=0, a_2=1, a_3=1):
    """
   Calculate the Lorentzian shape.

   Parameters:
       x (array-like): The input variable.
      a_1 (float): The center of the Lorentzian.
      a_2: The Full Width at Half Maximum (FWHM) of the Lorentzian.
      a_3(float): The amplitude of the Gaussian.

   Returns:
       array-like: The Lorentzian shape evaluated at x.
   """
    gamma = a_2/ 2
    return (a_3* gamma**2) / ((x - a_1)**2 + gamma**2)


def gauss_lorentz_f(x, a_0=1, a_1=0, a_2=1, a_3=0):
    """
   Calculate the pseudo Voigt as a Gauss-Lorentzian shape.

   Parameters:
       x (array-like): The input variable.
       a_0 (float): The componenet of Gaussian in the mix.
       a_1 (float): The center of the Gaussian.
       a_2: The Full Width at Half Maximum (FWHM) of the Gaussian.
       a_3(float): The amplitude of the Gaussian.
       From a_0 is possible to approximate  the FWHM of
       the lorentzian (fₗ) in the Voigh to ~1% from:
           a_0= 1.36603(fₗ/a_2) - 0.47719(fₗ/a_2)² + 0.11116(fₗ/a_2)³
       andthe Gaussian part (f₉) from:
          a_0 = [f₉⁵ + 2.69269f₉⁴fₗ + 2.42843f₉³fₗ² + 4.47163f₉²fₗ³ + 0.07842f₉fₗ⁴ + fₗ⁵]^(1/5)
   Returns:
       array-like: The Lorentzian shape evaluated at x.
   """

    f1 = a_0 *gauss_f(x, a_1, a_2, a_3)
    f2 = (1 - a_0) * lorentz_f(x,a_1, a_2, a_3)
    return  f1 + f2

def voigt_f_not_normalised(x, x0, g_FWHM, l_FWHM):
    """
    Calculate the Voigt function.

    The Voigt function represents the convolution of a Gaussian and a Lorentzian function.

    Parameters:
        x (array-like): The input variable.
        x0 (float): The center of the Voigt function.
        g_FWHM (float): The full width at half maximum (FWHM) of the Gaussian component.
        l_FWHM (float): The full width at half maximum (FWHM) of the Lorentzian component.

    Returns:
        array-like: The Voigt function evaluated at x.

    Notes:
        The Voigt function is computed using the Faddeeva function (wofz) from the scipy library.

    References:
        - Faddeeva, W., & Faddeeva, S. (2010). Algorithm 916: Evaluation of the Complex Error Function. ACM Transactions on Mathematical Software (TOMS), 40(3), 15:1-15:16. doi:10.1145/1824820.1824821

    """

    x_i = x - x0
    alpha = np.maximum(g_FWHM / 2,1e-4)
    gamma =np.maximum(l_FWHM / 2,1e-4)
    sigma = alpha / np.sqrt(2 * np.log(2))

    # Calculate the complex argument for the Faddeeva function
    z = (x_i + 1j * gamma) / sigma / np.sqrt(2)

    # Evaluate the Faddeeva function to compute the Voigt function
    faddeeva = wofz(z)

    # Compute the real part of the Faddeeva function and normalize by the standard deviation
    voigt = np.real(faddeeva) / sigma / np.sqrt(2 * np.pi)

    return voigt
def voigt_f(x, x0, g_FWHM, l_FWHM,amplitude):
    """
    Calculate the normalized Voigt function by normalizing the peak.

    Parameters:
        x (array-like): The input variable.
        x0 (float): The center of the Voigt function.
        g_FWHM (float): The full width at half maximum (FWHM) of the Gaussian component.
        l_FWHM (float): The full width at half maximum (FWHM) of the Lorentzian component.
        amplitude (float): The amplitude of the lineshape.
    Returns:
        array-like: The normalized Voigt function evaluated at x.

    """

    # Calculate the unnormalized Voigt function
    voigt = voigt_f_not_normalised(x, x0, g_FWHM, l_FWHM)

    # Find the maximum value of the Voigt function at x0
    max_value = np.maximum(voigt_f_not_normalised(x0, x0, g_FWHM, l_FWHM),1e-20)

    # Normalize the Voigt function by dividing by the maximum value
    voigt_normalized = voigt / max_value

    return voigt_normalized*amplitude


def fano_f_not_normalised(x, a_1=0, a_2=1, q=1e10):
    """
    Calculate the Fano lineshape without Voigt

    Parameters:
        x (array-like): The input variable.
        a_1 (float): The center of the Fano lineshape.
        a_2 (float): The Full Width at Half Maximum (FWHM) of the Fano lineshape.
        q (float): The asymmetry parameter (Fano factor) determining the shape of the line.

    Returns:
        array-like: The Fano lineshape evaluated at x.
    """
    gamma = a_2 / 2
    epsilon = (x - a_1) / gamma

    fano_par=-1*q


    unscaled=((1-fano_par**2)-2*fano_par*epsilon)/(1+epsilon**2)


    return unscaled

def fano_f(x, a_1=0, a_2=1, a_3=1, q=1e10):
    """
    Calculate the Fano lineshape without Voigt

    Parameters:
        x (array-like): The input variable.
        a_1 (float): The center of the Fano lineshape.
        a_2 (float): The Full Width at Half Maximum (FWHM) of the Fano lineshape.
        a_3 (float): The amplitude of the lineshape.
        q (float): The asymmetry parameter (Fano factor) determining the shape of the line.

    Returns:
        array-like: The Fano lineshape evaluated at x.
    """
    # Calculate the unnormalized Voigt function
    fano = fano_f_not_normalised(x, a_1, a_2,q)

    # Find the maximum value of the Voigt function at x0
    max_value = np.maximum(fano_f_not_normalised(a_1, a_1, a_2,q),1e-20)

    # Normalize the Voigt function by dividing by the maximum value
    fano_normalized = fano / max_value

    return fano_normalized*a_3


def voigt_fano_not_normalised(x, x0, g_FWHM, l_FWHM, q=0):
    """
    Calculate the Voigt profile using the convolution of the Fano lineshape and a Gaussian.

    Parameters:
        x (array-like): The input variable.
        peak_pos (float): The center of the profile.
        g_FWHM (float): The Gaussian component's full width at half maximum (FWHM).
        l_FWHM (float): The Fano lineshape component's full width at half maximum (FWHM).
        q (float): The asymmetry parameter (Fano factor) determining the shape of the line.

    Returns:
        array-like: The Voigt profile evaluated at x.
    """


    x_i = x - x0
    alpha = np.maximum(g_FWHM / 2,1e-4)
    gamma =np.maximum(l_FWHM / 2,1e-4)
    sigma = alpha / np.sqrt(2 * np.log(2))

    # Calculate the complex argument for the Faddeeva function
    z = (x_i + 1j * gamma) / sigma / np.sqrt(2)

    # Evaluate the Faddeeva function to compute the Voigt function
    faddeeva = wofz(z)

    # Compute the real part of the Faddeeva function and normalize by the standard deviation
    voigt_r = np.real(faddeeva)/ sigma / np.sqrt(2 * np.pi)
    voigt_imag = np.imag(faddeeva) / sigma / np.sqrt(2 * np.pi)

    fano_par=q

    unscaled=((1-fano_par**2)*voigt_r+2*fano_par*voigt_imag)


    return unscaled

def voigt_fano_f(x, x0, g_FWHM, l_FWHM, amplitude, q=0):
     """
     Calculate the normalised Voigt profile using the convolution of the Fano lineshape and a Gaussian.

     Parameters:
         x (array-like): The input variable.
         peak_pos (float): The center of the profile.
         g_FWHM (float): The Gaussian component's full width at half maximum (FWHM).
         l_FWHM (float): The Fano lineshape component's full width at half maximum (FWHM).
         amplitude (float): The amplitude of the Voigt profile.
         q (float): The asymmetry parameter (Fano factor) determining the shape of the line.

     Returns:
         array-like: The Voigt profile evaluated at x.
     """
     # Calculate the unnormalized Voigt function
     voight_fano = voigt_fano_not_normalised(x, x0, g_FWHM, l_FWHM, q)

     # Find the maximum value of the Voigt function at x0
     max_value = np.maximum(voigt_fano_not_normalised(x0, x0, g_FWHM, l_FWHM, q),1e-20)

     # Normalize the Voigt function by dividing by the maximum value
     voigt_fano_normalized = voight_fano / max_value

     return voigt_fano_normalized*amplitude



def model_f(params, x, peaks,model_type=None):
    """
    Calculate the simple model function consisting of set of peak functions without baseline.

    Parameters:
        params (lmfit.Parameters): The parameters object containing the fitting parameters.
        x (array-like): The input variable.
        peaks (list): The list of peak positions, and intensities.
        model_type (list, optional): The type of peak model to use for each peak, the array must have the same size as peaks.
                                   Options are 'Gaussian', 'Lorentz', 'Gauss-Lorentz', 'Voigt',
                                   'Fano-Simply', and 'Fano-Full'. Default is ' all Gaussian'.

    Returns:
        array-like: The calculated model function evaluated at x.
    """


    function_composed=[]
    if model_type==None:
        model_type=["Gaussian" for item in peaks]

    for  item in range(len(peaks)):
        # Load the parameters for the peaks
        if model_type[item]=="Gaussian":
            function_composed.append(gauss_f(x,
                                      params['Peak_'+str(item+1)+'_Center'],
                                      params['Peak_'+str(item+1)+'_FWHM'],
                                      params['Peak_'+str(item+1)+'_Intensity'])
                                     )
        elif model_type[item]=="Lorentz" :
            function_composed.append(lorentz_f(x,
                                      params['Peak_'+str(item+1)+'_Center'],
                                      params['Peak_'+str(item+1)+'_FWHM'],
                                      params['Peak_'+str(item+1)+'_Intensity'])
                                     )
        elif model_type[item]=="Gauss-Lorentz":
            function_composed.append(gauss_lorentz_f(x,
                                      params['Peak_'+str(item+1)+'_Gaussian_component'],
                                      params['Peak_'+str(item+1)+'_Center'],
                                      params['Peak_'+str(item+1)+'_FWHM'],
                                      params['Peak_'+str(item+1)+'_Intensity'])
                                     )
        elif model_type[item]=="Voigt":
            function_composed.append(voigt_f(x,
                                      params['Peak_'+str(item+1)+'_Center'],
                                      params['Peak_'+str(item+1)+'_Gauss_FWHM'],
                                      params['Peak_'+str(item+1)+'_Lorentz_FWHM'],
                                      params['Peak_'+str(item+1)+'_Intensity'])
                                     )
        elif model_type[item]=="Fano-Simply":
            function_composed.append(fano_f(x,
                                      params['Peak_'+str(item+1)+'_Center'],
                                      params['Peak_'+str(item+1)+'_Fano_FWHM'],
                                      params['Peak_'+str(item+1)+'_Intensity'],
                                      params['Peak_'+str(item+1)+'_Fano_Asymmetry'])
                                     )
        elif model_type[item]=="Fano-Full":
            function_composed.append(voigt_fano_f(x,
                                      params['Peak_'+str(item+1)+'_Center'],
                                      params['Peak_'+str(item+1)+'_Gauss_FWHM'],
                                      params['Peak_'+str(item+1)+'_Fano_FWHM'],
                                      params['Peak_'+str(item+1)+'_Intensity'],
                                      params['Peak_'+str(item+1)+'_Fano_Asymmetry'])
                                     )
        elif model_type[item]=='Not used':

            f_0(x)

        else:
            raise ValueError("Wrong model.")

    peak_term=np.sum(function_composed,
                     axis=0)

    return peak_term


def params_f(peaks, model_type=None):
    """
   Set up the parameters for a fit model.

   Args:
       peaks (list): List of peak positions and intensities.
       model_type (list, optional): The type of peak model to use for each peak, the array must have the same size as peaks.
                                  Options are 'Gaussian', 'Lorentz', 'Gauss-Lorentz', 'Voigt',
                                  'Fano-Simply', and 'Fano-Full'. Default is ' all Gaussian'.

   Returns:
       params (Parameters): Parameters object for the fit.
   """
    # Set up the parameters for the fit
    params = Parameters()



    if model_type==None:
        model_type=["Gaussian" for item in peaks]


    for item in range(len(peaks)):
        # Load the parameters for the peaks

        if model_type[item]=="Gaussian":
            params.add('Peak_'+str(item+1)+'_Center', value=peaks[item][0],min=0)
            params.add('Peak_'+str(item+1)+'_FWHM', value=5,min=1)
            params.add('Peak_'+str(item+1)+'_Intensity', value=peaks[item][1],min=0)
        elif model_type[item]=="Lorentz" :
            params.add('Peak_'+str(item+1)+'_Center', value=peaks[item][0],min=0)
            params.add('Peak_'+str(item+1)+'_FWHM', value=5,min=1)
            params.add('Peak_'+str(item+1)+'_Intensity', value=peaks[item][1],min=0)
        elif model_type[item]=="Gauss-Lorentz":
            params.add('Peak_'+str(item+1)+'_Gaussian_component', value=0.5,min=0,max=1)
            params.add('Peak_'+str(item+1)+'_Center', value=peaks[item][0],min=0)
            params.add('Peak_'+str(item+1)+'_FWHM', value=5,min=1)
            params.add('Peak_'+str(item+1)+'_Intensity', value=peaks[item][1],min=0)
        elif model_type[item]=="Voigt":
            params.add('Peak_'+str(item+1)+'_Center', value=peaks[item][0],min=0)
            params.add('Peak_'+str(item+1)+'_Gauss_FWHM', value=5,min=0.1)
            params.add('Peak_'+str(item+1)+'_Lorentz_FWHM', value=5,min=0.1)
            params.add('Peak_'+str(item+1)+'_Intensity', value=peaks[item][1],min=0)
        elif model_type[item]=="Fano-Simply":
            params.add('Peak_'+str(item+1)+'_Center', value=peaks[item][0],min=0)
            params.add('Peak_'+str(item+1)+'_Fano_FWHM', value=5,min=0.1)
            params.add('Peak_'+str(item+1)+'_Intensity', value=peaks[item][1],min=0)
            params.add('Peak_'+str(item+1)+'_Fano_Asymmetry', value=0,min=-1e4,max=1e4)
        elif model_type[item]=="Fano-Full":
            params.add('Peak_'+str(item+1)+'_Center', value=peaks[item][0],min=0)
            params.add('Peak_'+str(item+1)+'_Gauss_FWHM', value=5,min=0.1)
            params.add('Peak_'+str(item+1)+'_Fano_FWHM', value=5,min=0.1)
            params.add('Peak_'+str(item+1)+'_Intensity', value=peaks[item][1],min=0)
            params.add('Peak_'+str(item+1)+'_Fano_Asymmetry',  value=0,min=-1e4,max=1e4)
        elif model_type[item]=='Not used':
            print("not used")
        else:
            raise ValueError("Wrong model.")

    return params


def objective(params, x, data=None,peaks=None,model_type=None):
    """
   Objective function for fitting a model to data.

   Parameters:
       params (lmfit.Parameters): Model parameters to be optimized.
       x (array-like): Independent variable data.
       data (array-like): Dependent variable data to fit the model to.
       peaks (list): List of peak positions and intensities.
       model_type (list, optional): The type of peak model to use for each peak, the array must have the same size as peaks.
                                  Options are 'Gaussian', 'Lorentz', 'Gauss-Lorentz', 'Voigt',
                                  'Fano-Simply', and 'Fano-Full'. Default is ' all Gaussian'.

   Returns:
       array-like: Difference between the model values and the data.
   """

    model_values = model_f(params, x, peaks,model_type)


    return model_values - data


def fit(x,y,peaks,model_type=None):
    """
   Fits the given data with the specified model.

   Parameters:
       x (array-like): The x-values of the data.
       y (array-like): The y-values of the data.
       order (int): The order of the polynomial baseline.
       peaks (list): List of peak positions and intensities.
       model_type (list, optional): The type of peak model to use for each peak, the array must have the same size as peaks.
                                  Options are 'Gaussian', 'Lorentz', 'Gauss-Lorentz', 'Voigt',
                                  'Fano-Simply', and 'Fano-Full'. Default is ' all Gaussian'.

   Returns:
       object: The result of the fitting process.

   """
    pars=params_f(peaks,model_type)


    minimizer = Minimizer(objective, pars, fcn_args=(x, y,peaks,model_type))
    result = minimizer.least_squares(**{'xtol': 1e-6,
                                   'gtol': 1e-6,
                                   'ftol':1e-6,
                                   'max_nfev':5000})
    return result

def fit_info(fit):
    return(lmfit.fit_report(fit))



def test():


    # Generate example data
    x = np.linspace(0, 2000, 1000)
    y =voigt_f(x, 500, 50,50,1)#-voigt_fano_f(x,500, 50,50,1, -1/5) #fano_f(x,40, 4, 1,1/8)+voigt_fano_f(x, 80, 4,4, 1,1/8) # True underlying function

    peaks=[[500,10]]#,[60.5,40000],[80,40000]]

    models=["Lorentz","Gaussian","Fano-Simply"]
    fitting=fit(x,y,
                peaks,
                models)
    print(lmfit.fit_report(fitting))

    new_y=model_f(fitting.params, x,peaks,models)

    Raman_plot.plotter([[x,y],[x,new_y]],
                   ['x','y'],
                   'none',
                   leyends=['raw','fit'],
                    lines=True)



#test()

#test("Gaussian")
#test("Lorentz")
#test("Gauss-Lorentz")
#test("Voigt")
#test2()
