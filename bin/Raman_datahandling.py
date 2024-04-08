# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 16:13:16 2023
 Library to remove baselines and process the data
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

# Raman_datahandling

from numpy import abs, mean, sort, searchsorted, minimum, max, column_stack, asarray, concatenate, ones_like
from numpy import argsort, ones, std, clip, exp, diff, eye, allclose, sqrt, array, empty_like, transpose, expand_dims, where, sign, append
from scipy.signal import savgol_filter
from scipy.signal import find_peaks as sp_find_peaks
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve
from scipy.interpolate import PchipInterpolator
from scipy import sparse
from scipy.sparse import linalg
from numpy.linalg import norm
from scipy.interpolate import interp1d


from bin import Raman_plot
from bin import Raman_dataloader
from bin import Raman_fit

def smooth_spectra(x, y, window_size, poly_order,deriv_order=0):
    """
   Apply the Savitzky-Golay filter to smooth the spectra.

   Parameters:
   - x: The x-coordinates of the spectra.
   - y: The y-coordinates of the spectra.
   - window_size: The size of the smoothing window.
   - poly_order: The order of the polynomial to fit the data in the smoothing window.

   Returns:
   A list containing the x-coordinates and the smoothed y-coordinates.
   """
    # Apply the Savitzky-Golay filter to smooth the spectra
    smoothed_y = savgol_filter(y, window_size, poly_order,deriv=deriv_order, delta=y[1] - y[0])
    return [x, smoothed_y]
def set_below_threshold(arr, threshold):
    arr[abs(arr) < threshold] = 0
    return arr
def local_avg(y,peak,window):
    if peak-window>=0 and peak+window<=len(y)-1 :
        avg_l=mean(y[peak-window:peak])
        avg_r=mean(y[peak:peak+window])
    elif  peak-window<0 and peak+window<=len(y)-1:
        avg_l=mean(y[0-window:peak])
        avg_r=mean(y[peak:peak+window])
    elif peak-window>=0 and peak+window>len(y)-1 :
        avg_l=mean(y[peak-window:peak])
        avg_r=mean(y[peak:len(y)-1])
    return avg_l,avg_r

def neighbours(y, array, n):
    """
    Find the closest smaller and closest bigger elements in a sorted array compared to a given value.

    Parameters:
        y (ndarray): The original array from which `array` is derived.
        array (ndarray): A sorted array.
        n (float): The value to compare with.

    Returns:
        list: A list containing the closest smaller and closest bigger elements.

    """
    sorted_arr = sort(array)

    # Find index for insertion
    insert_index = searchsorted(sorted_arr, n)

    # Find closest smaller and closest bigger elements
    if insert_index - 1 > 0 and insert_index + 1 < len(array) - 1:
        closest_smaller = sorted_arr[insert_index - 1]
        closest_bigger = sorted_arr[insert_index + 1]
    elif insert_index - 1 <= 0:
        closest_smaller = sorted_arr[0]
        closest_bigger = sorted_arr[insert_index + 1]
    elif insert_index + 1 >= len(array) - 1:
        closest_smaller = sorted_arr[insert_index - 1]
        closest_bigger = sorted_arr[len(array) - 1]

    min_extremes = 3
    if closest_smaller < min_extremes:
        closest_smaller = 3
    if closest_bigger > len(y) - 1 - min_extremes:
        closest_bigger = len(y) - 1 - min_extremes

    return [closest_smaller, closest_bigger]


def peak_classification(y, peak, window):
    """
    Classify a peak as a local maximum or minimum based on its comparison with local average values.

    Parameters:
        y (ndarray): The input array.
        peak (int): The index of the peak to classify.
        window (int): The size of the window for calculating local averages.

    Returns:
        int: 1 if the peak is a local maximum, -1 if it is a local minimum.

    """
    avg_l, avg_r = local_avg(y, peak, window)
    min_avd = minimum(avg_l, avg_r)

    if y[peak] > min_avd:
        return 1  # Return 1 if it is a local maximum
    else:
        return -1  # Return -1 if it is a local minimum


def find_maxmin(x, y, window):
    """
    Find the indices of local maxima and minima in a given spectrum.

    Parameters:
        x (ndarray): The x-values of the spectrum.
        y (ndarray): The y-values of the spectrum.
        window (int): The window size for smoothing.

    Returns:
        list: A list of indices representing the local maxima and minima.

    """
    x_d, y_d = smooth_spectra(x, y, window, 2, 1)  # Derivative
    #tol = max(y_d) / 10000  # All with a very small derivative is filtered
    #y_filter = set_below_threshold(y_d, tol)

    peaks = where(diff(sign(y_d)) != 0)[0] + 1

    return [peak for peak in peaks]

def find_peaks_alt(x, y, window,threshold):
    """
    Find the indices of peaks and valleys in a given spectrum.

    Parameters:
        x (ndarray): The x-values of the spectrum.
        y (ndarray): The y-values of the spectrum.
        window (int): The window size for smoothing.

    Returns:
        tuple: A tuple containing two lists. The first list represents the indices of peaks, and the second list
               represents the indices of valleys.

    """
    max_min_index = find_maxmin(x, y, window)
    peak_class = [peak_classification(y, peak, 3) for peak in max_min_index]  # Min of close neighbourhood
    peaks = [max_min_index[item] for item in range(len(max_min_index)) if peak_class[item] > 0]
    valleys = [max_min_index[item] for item in range(len(max_min_index)) if peak_class[item] < 0]

    #Get maximum value for peak
    peak_values=[y[peak] for peak in peaks]
    max_peak=max(peak_values)
    filter_peaks=[peaks[item] for item in range(len(peak_values)) if peak_values[item]/(max_peak+1e-12)>threshold]
    return filter_peaks, valleys

def find_peaks(x, y, window, threshold, prominence=None,distance=None):
    """
    Find peaks in a spectrum based on smoothed data and optional threshold, prominence, and distance.

    Parameters:
        x (array-like): The x-coordinates of the spectrum.
        y (array-like): The y-coordinates of the spectrum.
        window (int): The size of the window used for smoothing.
        threshold (float): The threshold for peak detection.
        prominence (float): The minimum prominence required for a peak (optional).
        distance (float): The minimum distance between peaks (optional).

    Returns:
        tuple: A tuple containing two lists:
            - A list of indices representing the filtered peaks.
            - An empty list (dummy variable for compatibility with previous method).

    """

    # Smooth the spectrum
    x_d, smoothed = smooth_spectra(x, y, window, 2, deriv_order=0)

    smoothed=smoothed/max(smoothed)

    # Find peaks in the smoothed data
    peaks, _ = sp_find_peaks(smoothed, prominence=prominence, distance=distance)

    # Get maximum value for peak
    peak_values = [y[peak] for peak in peaks]

    max_peak = max(peak_values)

    # Filter peaks based on the threshold
    filter_peaks = [peaks[item] for item in range(len(peak_values)) if (peak_values[item] / (max_peak + 1e-12) > threshold) ]

    # Order the peaks:

    valleys = []  # Dummy variable to maintain compatibility with the previous method

    return filter_peaks, valleys




def peak_asymmetry_index(x, y_wobaseline, peak, height_to_check):
    """
    Calculate the asymmetry of a peak.

    Note: This function assumes that the baseline has already been removed from the data.

    Parameters:
        x (array-like): The x-values of the data.
        y_wobaseline (array-like): The y-values of the data without the baseline.
        peak (int): The index of the peak in the data.
        height_to_check (float): The fraction of the peak's height to descend and check for symmetry.

    Returns:
        list: A list containing the indices of the left and right boundaries of the peak's Full Width
        at the especified heigth.
    """

    # Left side, descend to half the value and check symmetry
    control_l = True
    iterator = 0
    while control_l:
        if y_wobaseline[peak - iterator] > y_wobaseline[peak] * height_to_check:
            if peak - iterator > 0:
                iterator += 1
            else:
                control_l = False
                left_FWHM = peak - iterator
        else:
            left_FWHM = peak - iterator
            control_l = False

    # Right side, descend to half the value and check symmetry
    control_r = True
    iterator = 0
    while control_r:
        if y_wobaseline[peak + iterator] > y_wobaseline[peak] * height_to_check:
            if peak + iterator < len(y_wobaseline) - 1:
                iterator += 1
                right_FWHM = peak + iterator
            else:
                control_r = False
        else:
            right_FWHM = peak + iterator
            control_r = False

    return [left_FWHM, right_FWHM]
def peak_side_ratio(x,peak,boundaries):
    left_side=x[peak]-x[boundaries[0]]
    right_side=x[boundaries[1]]-x[peak]
    ratio=right_side/(left_side+1e-20)
    return ratio

def peak_type_clasification(x, y_wobaseline,peaks):
    def inner_class(val):
        val_co=abs(1-val)
        if val_co<0.05:
            return False
        else:
            return True

    peak_boun_a=[inner_class(peak_side_ratio(x,peak,peak_asymmetry_index(x, y_wobaseline,peak,0.25)))
                 for peak in peaks]
    peak_boun_b=[inner_class(peak_side_ratio(x,peak,peak_asymmetry_index(x, y_wobaseline,peak,0.5)))
                 for peak in peaks]
    peak_boun_c=[inner_class(peak_side_ratio(x,peak,peak_asymmetry_index(x, y_wobaseline,peak,0.75)))
                 for peak in peaks]
    peak_boun_d=[inner_class(peak_side_ratio(x,peak,peak_asymmetry_index(x, y_wobaseline,peak,0.9)))
                 for peak in peaks]
    peak_boun_e=[inner_class(peak_side_ratio(x,peak,peak_asymmetry_index(x, y_wobaseline,peak,0.95)))
                 for peak in peaks]
    output=column_stack((asarray(peak_boun_a),
                             asarray(peak_boun_b),
                             asarray(peak_boun_c),
                             asarray(peak_boun_d),
                             asarray(peak_boun_e)
                             ))
    print(output)








def find_window(x, y, peak, window):
    """
    Find the indices of the closest smaller and closest bigger elements to a given peak index.

    Parameters:
        x (ndarray): The x-values of the spectrum.
        y (ndarray): The y-values of the spectrum.
        peak (int): The index of the peak.
        window (int): The window size for smoothing.

    Returns:
        list: A list containing the indices of the closest smaller and closest bigger elements.

    """
    x_d, y_d = smooth_spectra(x, y, window, 2, 1)  # Derivative
    peaks, valleys = find_peaks(x_d, y_d, window)
    combined = concatenate((peaks, valleys))
    out = neighbours(y, combined, peak)

    return out


def exclude_indexes(x, y, exclude_indexes):
    """
    Exclude specific index ranges from the given x and y arrays.

    Parameters:
        x (ndarray): The x-values of the spectrum.
        y (ndarray): The y-values of the spectrum.
        exclude_indexes (list): A list of index ranges to exclude. Each range is defined as a tuple of start and end indices.

    Returns:
        tuple: A tuple containing the filtered x and y arrays.

    """
    mask = ones_like(x, dtype=bool)

    for exclusion_range in exclude_indexes:
        start_index, end_index = exclusion_range
        mask[start_index:end_index + 1] = False

    return x[mask], y[mask]


def fit_piecewise_spline_with_interpolation(x, y, ignore_regions):
    """
    Perform piecewise-cubic interpolation with continuous first derivatives.


    Parameters:
        x (array-like): The x-values of the data.
        y (array-like): The y-values of the data.
        ignore_regions (list of tuples): List of regions [x_min, x_max] to be ignored during the fitting process.

    Returns:
        PchipInterpolator: The piecewise-cubic interpolation object.
    """

    filter_x,filter_y=exclude_indexes(x, y, ignore_regions)
    # Sort the data points by x-values (if not already sorted)
    sorted_indices = argsort(filter_x)
    x_sorted = filter_x[sorted_indices]
    y_sorted = filter_y[sorted_indices]

    # Perform piecewise-cubic interpolation
    interpolator = PchipInterpolator(x_sorted, y_sorted)

    return interpolator



#DCO implementation of als
def baseline_arPLS(y, ratio=1e-6, lam=1e5, niter=1000, full_output=False):
    """
    Perform asymmetrically reweighted penalized least squares (arPLS) baseline correction on a given spectrum.
    Sung-June Baek, Aaron Park, Young-Jin Ahna and Jaebum Choo: "Baseline correction using asymmetrically reweighted penalized least squares smoothing", Analyst, 2015,140, 250-257
    Here's how it works:

    Initially, the algorithm assumes that the input spectrum y consists of both the baseline and the peaks.
    
    In each iteration, the algorithm estimates the baseline z by solving a linear system (W + H)z = Wy, where W is the weight matrix and H is the smoothing matrix. This estimation of the baseline is based on the current weights assigned to the data points.
    
    Once the baseline z is estimated, the algorithm calculates the residuals d by subtracting the estimated baseline from the input spectrum: d = y - z.
    
    The residuals represent the deviations from the estimated baseline. If a point in the spectrum corresponds to a peak or feature, its value in the residuals would be significantly different from zero, as it deviates from the estimated baseline.
    
    The algorithm focuses on the negative residuals dn, which are the points in the residuals that fall below the estimated baseline. These negative residuals are likely to correspond to the peaks in the original spectrum.
    
    The algorithm calculates the mean m and standard deviation s of the negative residuals dn. These statistical measures provide information about the magnitude and distribution of the peaks in relation to the baseline.
    
    Based on the exponential decay weighting scheme, the algorithm assigns new weights w_new to the data points. The calculation of the weights depends on the values of the residuals d, the mean m, and the standard deviation s. The idea is that points with large negative residuals (corresponding to peaks) receive smaller weights, while points closer to the baseline receive higher weights.
    
    The convergence criterion is computed based on the difference between the new weights w_new and the previous weights w. If the weights change negligibly, it indicates that the algorithm has successfully captured the baseline and is effectively ignoring the peaks.
    
    By updating the weights in each iteration, the algorithm gradually adapts to the presence of the peaks and assigns lower weights to those data points, effectively ignoring their influence on the estimated baseline.
    
    The iterative process continues until the convergence criterion falls below a specified threshold or the maximum number of iterations is reached. At this point, the estimated baseline represents the smooth trend of the spectrum, excluding the peaks or features.

    Parameters:
        y (ndarray): The input spectrum.
        ratio (float, optional): Convergence threshold for stopping iterations. Defaults to 1e-6.
        lam (float, optional): Smoothing parameter. Higher values result in smoother baselines. Defaults to 100.
        niter (int, optional): Maximum number of iterations. Defaults to 10.
        full_output (bool, optional): Flag to indicate whether to return additional information or just the baseline.
                                      Defaults to False.

    Returns:
        ndarray or tuple: If `full_output` is False, returns the estimated baseline as a 1D array.
                         If `full_output` is True, returns a tuple containing the estimated baseline, the residuals,
                         and additional information as a dictionary.

    Raises:
        ValueError: If the input spectrum `y` is empty or has a length less than 3.

    """

    L = len(y)

    if L < 3:
        raise ValueError("Input spectrum `y` must have a length greater than or equal to 3.")

    # Create the difference matrix D for calculating the second derivative
    diag = ones(L - 2)
    D = sparse.spdiags([diag, -2 * diag, diag], [0, -1, -2], L, L - 2)

    # Calculate the smoothing matrix H using the difference matrix D
    H = lam * D.dot(D.T)

    # Initialize the weights and the weight matrix W
    w = ones(L)
    W = sparse.spdiags(w, 0, L, L)

    crit = 1
    count = 0

    while crit > ratio:
        # Solve the linear system (W + H)z = Wy to estimate the baseline z
        z = linalg.spsolve(W + H, W * y)

        # Calculate the residuals d = y - z
        d = y - z

        # Calculate the negative residuals dn
        dn = d[d < 0]

        # Calculate the mean and standard deviation of the negative residuals
        m = mean(dn)
        s = std(dn)

       # Calculate the updated weights based on the exponential decay weighting
        exp_input = 2 * (d - (2 * s - m)) / s
        exp_input = clip(exp_input, -500, 500)  # Clip the values to a manageable range
        w_new = 1 / (1 + exp(exp_input))

        # Calculate the convergence criterion
        crit = norm(w_new - w) / norm(w)

        # Update the weights and the weight matrix
        w = w_new
        W.setdiag(w)

        count += 1

        if count > niter:
            print('Maximum number of iterations exceeded')
            break

    if full_output:
        # Return the baseline, residuals, and additional information
        info = {'num_iter': count, 'stop_criterion': crit}
        return z, d, info
    else:
        # Return only the baseline
        return z


def whittaker_smoother(y, lambda_, differences=2):
    """

    Parameters
    ----------
    y : double/float
       y data of spectra
    lambda_ : TYPE
       smoothing parameter
    differences : TYPE, optional
        DESCRIPTION. The default is 2.

    Returns
    -------
    z : double/float
        baseline
        The algorithm is based on a penalized least squares approach, where a balance is struck between fitting the data and achieving smoothness in the resulting baseline estimate.

Here's a more detailed explanation of the Whittaker Smoother algorithm:

Input Data: The algorithm takes as input a one-dimensional array representing the spectrum to be processed. The spectrum may contain both peaks and a baseline component.

Initial Estimate: The algorithm starts with an initial estimate of the baseline, which can be a simple linear fit or any other reasonable initial guess.

Iteration: The Whittaker Smoother algorithm performs an iterative process to refine the baseline estimate. In each iteration:

a. Smoothing Step: The current estimate of the baseline is smoothed using a smoothing parameter (λ) and a second-order difference penalty term. This step helps in achieving a smooth baseline while preserving the underlying features.

b. Weighted Least Squares Step: The smoothed baseline estimate is combined with the original spectrum to create a weighted least squares problem. The spectrum values are given higher weights near the peaks to ensure that the peaks are not significantly affected during the baseline estimation process.

c. Solving the Least Squares Problem: The weighted least squares problem is solved to obtain an updated estimate of the baseline.

Convergence: The iteration process continues until the difference between consecutive baseline estimates falls below a predefined threshold or until a maximum number of iterations is reached.

Output: The final result of the Whittaker Smoother algorithm is an estimate of the baseline that has been smoothed while preserving the underlying peaks.

The Whittaker Smoother algorithm is effective in removing baseline variations in Raman spectra, especially when the peaks have irregular shapes and the baseline exhibits complex variations. The smoothing parameter (λ) controls the trade-off between smoothness and fitting to the data, and its value needs to be adjusted based on the characteristics of the spectrum.

    """
    n = len(y)
    d = diff(eye(n), differences)
    w = ones(n)

    for i in range(10):
        W = spdiags(w, 0, n, n)
        Z = W + lambda_ * d.dot(d.transpose())
        z = spsolve(Z, w * y)
        w_new = 1.0 / (abs(y - z) + 1e-6)

        if allclose(w, w_new):
            break
        w = w_new

    return z


def snr(y):
    """
   Calculate the Signal-to-Noise Ratio (SNR) of a signal.

   Parameters:
   - y: The input signal.

   Returns:
   The calculated SNR value.
   """
    avg = mean(y)
    dif = sqrt((y-avg)**2)
    snr = mean(dif)

    return snr

def compute_snr(data, window_size):
    """
    Compute the Signal-to-Noise Ratio (SNR) for windows of a given size and find the window with the smallest SNR.

    Parameters:
        data (array-like): The input data.
        window_size (int): The size of the sliding window.

    Returns:
        min_snr_window (array): The window with the smallest SNR.
        min_snr (float): The smallest SNR value.
    """
    # Calculate the number of windows
    if window_size >= len(data):
        num_windows=1
        window_size=len(data)-1
    else:
        num_windows = len(data) - window_size + 1

    # Calculate the SNR for each window
    window_snr = array([snr(data[i:i+window_size]) for i in range(num_windows)], dtype=object)
    window_snr = append(window_snr, snr(data[num_windows:]))

    # Get the smallest SNR value
    min_snr = min(window_snr)
    max_snr = max(window_snr)

    return min_snr, max_snr

def select_data_region(x, y, window=[0, -1]):
    """
    Select a specific region of data from the given x and y arrays.

    Parameters:
        x (ndarray): The x-values of the spectrum.
        y (ndarray): The y-values of the spectrum.
        window (list): A list defining the start and end indices of the desired region. The default value is [0, -1],
                       which includes the entire range.

    Returns:
        tuple: A tuple containing the selected x and y arrays.

    """
    try:
        if window[0] >= 0 and window[1] <= len(x) - 1 and window[1] > window[0]:
            return x[window[0]:window[1]], y[window[0]:window[1]]
        elif window[0] < 0 and window[1] <= len(x) - 1 and window[1] > window[0]:
            return x[0:window[1]], y[window[0]:window[1]]
        elif window[0] >= 0 and window[1] > len(x) - 1 and window[1] > window[0]:
            return x[0:window[1]], y[window[0] > len(x) - 1]
        else:
            raise ValueError("Invalid window range.")
    except ValueError:
        print("Invalid window range. Adjusting the window in steps of 1...")
        start = max(0, min(window[0], len(x) - 1))
        end = min(len(x), max(window[1], 0))

        for i in range(start, len(x)):
            for j in range(end, -1, -1):
                try:
                    if i < j:
                        return x[i:j], y[i:j]
                except ValueError:
                    pass

        print("Unable to find a valid window.")
        return x, y

def lin_baseline(points, x_range, model='linear'):

    """
    Perform piecewise linear interpolation of Y values within the specified X range using SciPy.

    Parameters:
    points (numpy.ndarray): An array of X, Y data points.
    x_range (tuple): A tuple containing the start and end X values for interpolation.

    Returns:
    numpy.ndarray: An array of interpolated Y values within the specified X range.
    """
    x_values = points[:, 0]
    y_values = points[:, 1]
    # Create a linear interpolation function for the data within the specified X range
    interp_func = interp1d(x_values,y_values, kind=model, fill_value='extrapolate')

    # Interpolate values within the specified X range

    interpolated_y_values = interp_func(x_range)

    return interpolated_y_values




def wavenumber_to_index(y, wn, tol):
    """
    Find the index of the closest wavenumber to a given value.

    Parameters:
        y (ndarray): The array of wavenumbers.
        wn (float): The wavenumber value to find.
        tol (float): The tolerance for the closeness of the wavenumber.

    Returns:
        int: The index of the closest wavenumber.

    Raises:
        ValueError: If no matching wavenumber is found within the specified tolerance.

    """
    max_attempts = 20
    index = where(abs(y - wn) < tol)[0]

    attempts = 0
    while len(index) == 0 and attempts < max_attempts:
        wn += 1  # Increase the wavenumber value
        index = where(abs(y - wn) < tol)[0]
        attempts += 1

    if len(index) == 0:
        raise ValueError("No matching wavenumber found within the specified tolerance.")

    return index[0]

path2=r"H:\OneDrive - UVa\Equipamiento- Raman Soleil\Datos Raw\Julian Anaya\2023\May\MA2TEC\Cellulose\MA2TEC\MCC Type I Powder\OptronLab\MCC Type I Powder _1_20 s_785nm_Edge_600 (500nm)_100x_500 µm_100% (39mW)_c_04.txt"
path1=r"H:\OneDrive - UVa\Equipamiento- Raman Soleil\Datos Raw\Julian Anaya\2023\May\MA2TEC\Cellulose\MA2TEC\MCC3pcNaOH Type II Powder\OptronLab\MCC3pcNaOH Type II I Powder _1_20 s_785nm_Edge_600 (500nm)_100x_500 µm_100% (39mW)_02.txt"


def test_fit(path,window):
    dat=Raman_dataloader.load_spectra_data(path)
    info=Raman_dataloader.load_spectra_info(path)
    x,y=select_data_region(dat[:,0],dat[:,1],[wavenumber_to_index(dat[:,0],window[0],1),wavenumber_to_index(dat[:,0],window[1],1)])
    y_baseline=baseline_arPLS(y, lam=1E5, niter=200,full_output=False)

    y_fix=y-y_baseline

    if norm:
        y_fix=y_fix/max(y_fix)+0

    peaks,valleys=find_peaks(x,y_fix,10,0.05)
    peak_val=[[x[peak],y_fix[peak]] for peak in peaks]
    print(peak_val)
    #Fit
    fit=Raman_fit.fit(x, y_fix, peak_val)

    y_fit=Raman_fit.model_f(fit.params, x, peak_val)
    print(Raman_fit.fit_info(fit))
    Raman_plot.plotter([[x,y_fix],[x,y_fit]
                        ],
                   [
                     "Raman shift "+"(1/cm)",
                     "Normalised Intens (A.U.)"],
                   "Luminescence in samples-20 s_785nm_Edge_600 (500nm)_100x_500 µm",
                   leyends= ['raw_norm','fit'],
                   lines=True,
                   res=1200,
                   size="double_size_double_heigh",
                   leyend_frame=True
                   )


#test_fit(path1,[400,1300])


def test_prep(iterator,norm=False):
    path1=r"H:\OneDrive - UVa\Equipamiento- Raman Soleil\Datos Raw\Julian Anaya\2023\May\MA2TEC\Cellulose\MA2TEC\MCC Type I Powder\OptronLab\MCC Type I Powder _1_20 s_785nm_Edge_600 (500nm)_100x_500 µm_100% (39mW)_c_03.txt"
    path2=r"H:\OneDrive - UVa\Equipamiento- Raman Soleil\Datos Raw\Julian Anaya\2023\May\MA2TEC\Cellulose\MA2TEC\MCC3pcNaOH Type II Powder\OptronLab\MCC3pcNaOH Type II I Powder _1_20 s_785nm_Edge_600 (500nm)_100x_500 µm_100% (39mW)_01.txt"
    path3=r"H:\OneDrive - UVa\Equipamiento- Raman Soleil\Datos Raw\Julian Anaya\2023\May\MA2TEC\Cellulose\MA2TEC\MCC20pcNaOH Type II Powder\OptronLab\MCC20pcNaOH Type II Powder _1_20 s_785nm_Edge_600 (500nm)_100x_500 µm_100% (39mW)_02.txt"
    path4=r"H:\OneDrive - UVa\Equipamiento- Raman Soleil\Datos Raw\Julian Anaya\2023\May\MA2TEC\Cellulose\MA2TEC\MW24 Type II Powder\OptronLab\MW24 Type II Powder _1_20 s_785nm_Edge_600 (500nm)_100x_500 µm_100% (39mW)_04.txt"
    path5=r"H:\OneDrive - UVa\Equipamiento- Raman Soleil\Datos Raw\Julian Anaya\2023\May\MA2TEC\Cellulose\MA2TEC\MW25 Type II Powder\OptronLab\MW25 Type II Powder _1_20 s_785nm_Edge_600 (500nm)_100x_500 µm_100% (39mW)_04.txt"
    path6=r"H:\OneDrive - UVa\Equipamiento- Raman Soleil\Datos Raw\Julian Anaya\2023\May\MA2TEC\Cellulose\MA2TEC\MW26 Type II Powder\OptronLab\MW26 Type II Powder _1_20 s_785nm_Edge_600 (500nm)_100x_500 µm_100% (39mW)_03.txt"
    path7=r"H:\OneDrive - UVa\Equipamiento- Raman Soleil\Datos Raw\Julian Anaya\2023\May\MA2TEC\Cellulose\MA2TEC\MW27 Type II Powder\OptronLab\MW27 Type II Powder _1_20 s_785nm_Edge_600 (500nm)_100x_500 µm_100% (39mW)_02.txt"
    path8=r"H:\OneDrive - UVa\Equipamiento- Raman Soleil\Datos Raw\Julian Anaya\2023\May\MA2TEC\Cellulose\MA2TEC\MW28 Type II Powder\OptronLab\MW28 Type II Powder _1_20 s_785nm_Edge_600 (500nm)_100x_500 µm_100% (39mW)_05.txt"
    path9=r"H:\OneDrive - UVa\Equipamiento- Raman Soleil\Datos Raw\Julian Anaya\2023\May\MA2TEC\Cellulose\MA2TEC\SHC Type II PowderY\OptronLab\SHC Type II Powder _1_20 s_785nm_Edge_600 (500nm)_100x_500 µm_10% (3_9mW)_01.txt"

    if iterator==0:

        path=path1
    elif iterator==1:

        path=path2
    elif iterator==2:

        path=path3
    elif iterator==3:

        path=path4
    elif iterator==4:

        path=path5
    elif iterator==5:

        path=path6
    elif iterator==6:

        path=path7
    elif iterator==7:

        path=path8
    elif iterator>7:

        path=path9

    dat=Raman_dataloader.load_spectra_data(path)
    info=Raman_dataloader.load_spectra_info(path)

    x,y=select_data_region(dat[:,0],dat[:,1],[wavenumber_to_index(dat[:,0],400,1),wavenumber_to_index(dat[:,0],1200,1)])

    y_baseline=baseline_arPLS(y, lam=1E5, niter=200,full_output=False)

    y_fix=y-y_baseline

    if norm:
        y_fix=y_fix/max(y_fix)+iterator

    peaks,valleys=find_peaks(x,y_fix,5,0.05)


    print(x[peaks])


    return [x,y_fix],[x,y_baseline]

def test1():
    all_dat=[test_prep(item)[1] for item in range(9)]
    legends=["MCC type I-100% power","MCC 3% NaOH-100% power","MCC 20% NaOH-100% power","MW24-100% power","MW25-100% power","MW26-100% power","MW27-100% power","MW28-100% power","SHC-10% power"]

    Raman_plot.plotter(all_dat,
                   [
                     "Raman shift "+"(1/cm)",
                     "Normalised Intens (A.U.)"],
                   "Luminescence in samples-20 s_785nm_Edge_600 (500nm)_100x_500 µm",
                   leyends= legends,
                   lines=True,
                   res=1200,
                   size="double_size_double_heigh",
                   leyend_frame=True
                   )
def test2():
    all_dat=[test_prep(item,True)[0] for item in range(9)]
    legends=["MCC type I-100% power","MCC 3% NaOH-100% power","MCC 20% NaOH-100% power","MW24-100% power","MW25-100% power","MW26-100% power","MW27-100% power","MW28-100% power","SHC-10% power"]

    Raman_plot.plotter(all_dat,
                   [
                     "Raman shift "+"(1/cm)",
                     "Normalised Intens (A.U.)"],
                   "Luminescence in samples-20 s_785nm_Edge_600 (500nm)_100x_500 µm",
                   leyends= legends,
                   lines=True,
                   res=1200,
                   size="double_size_double_heigh",
                   leyend_frame=True
                   )

def test3():
    path=r"E:\OneDrive - UVa\Univerdad Valladolid\Proyectos\Miscelanea\2023\Medidas Carina\Espectros\1st layer_785nm_Edge_600 (500nm)_500 µm_100x_10% (3_9mW)_30 s_1 a_01.txt"

    dat=Raman_dataloader.load_spectra_data(path)
    info=Raman_dataloader.load_spectra_info(path)

    x,y=select_data_region(dat[:,0],dat[:,1],[wavenumber_to_index(dat[:,0],500,1),wavenumber_to_index(dat[:,0],3500,1)])

    y_baseline=baseline_arPLS(y, lam=1E5, niter=200,full_output=False)

    y_fix=y-y_baseline
    peaks,valleys=find_peaks(x,y_fix,10,0.1)
    peak_val=[[x[peak],y_fix[peak]] for peak in peaks]
    print(peak_val)

    peak_asy=[peak_asymmetry_index(x,y_fix,peak,0.5) for peak in peaks]
    print(peak_asy)
    pos=[[x[item[0]],x[item[1]]] for item in peak_asy]
    print(pos)
    peak_type_clasification(x, y_fix,peaks)
    model_list=['Voigt' for item in peaks]
    fitting=Raman_fit.fit(x,y_fix,
                peak_val,model_list)
    print(Raman_fit.fit_info(fitting))

    new_y=Raman_fit.model_f(fitting.params, x,peak_val,model_list)

    #labelling of peaks
    arrow_data=[(
        peaks[item],
        peak_val[item][1]
        ) for item in range(len(peaks))
        ]
    Raman_plot.plotter([[x,y_fix],[x,new_y]],
                    [
                     "Raman shift"+"("+info['AxisUnit[0]']+")",
                     info['AxisType[1]']+" ("+info['AxisUnit[1]']+")"
                     ],
                     info['Title'],
                     leyends= ['Raw data'],
                   lines=True,
                   res=150,
                   size="double_size_double_heigh",
                   leyend_frame=[True,'b'],
                   arrow=arrow_data
                   )

def test4():
    path=r"E:\OneDrive - UVa\Univerdad Valladolid\Proyectos\Miscelanea\2023\Medidas Carina\Espectros\Resin high_785nm_Edge_600 (500nm)_500 µm_100x_10% (3_9mW)_30 s_1 a_01.txt"
    path=r"E:\OneDrive - UVa\Univerdad Valladolid\Proyectos\Miscelanea\2023\Medidas Carina\Espectros\1st layer_785nm_Edge_600 (500nm)_500 µm_100x_10% (3_9mW)_30 s_1 a_01.txt"
    path=r"E:\OneDrive - UVa\Univerdad Valladolid\Proyectos\Miscelanea\2023\Medidas Carina\Espectros\2nd layer_785nm_Edge_600 (500nm)_500 µm_100x_10% (3_9mW)_30 s_1 a_01.txt"
    path=r"E:\OneDrive - UVa\Univerdad Valladolid\Proyectos\Miscelanea\2023\Medidas Carina\Espectros\3nd layer_785nm_Edge_600 (500nm)_500 µm_100x_10% (3_9mW)_30 s_1 a_01.txt"
    path=r"E:\OneDrive - UVa\Univerdad Valladolid\Proyectos\Miscelanea\2023\Medidas Carina\Espectros\4nd layer_785nm_Edge_600 (500nm)_500 µm_100x_10% (3_9mW)_30 s_1 a_01.txt"
    path=r"E:\OneDrive - UVa\Univerdad Valladolid\Proyectos\Miscelanea\2023\Medidas Carina\Espectros\5nd layer_785nm_Edge_600 (500nm)_500 µm_100x_10% (3_9mW)_30 s_1 a_01.txt"
    path=r"E:\OneDrive - UVa\Univerdad Valladolid\Proyectos\Miscelanea\2023\Medidas Carina\Espectros\6nd layer_785nm_Edge_600 (500nm)_500 µm_100x_10% (3_9mW)_30 s_1 a_01.txt"
    path=r"E:\OneDrive - UVa\Univerdad Valladolid\Proyectos\Miscelanea\2023\Medidas Carina\Espectros\7nd layer_785nm_Edge_600 (500nm)_500 µm_100x_10% (3_9mW)_30 s_1 a_01.txt"
    path=r"E:\OneDrive - UVa\Univerdad Valladolid\Proyectos\Miscelanea\2023\Medidas Carina\Espectros\8nd layer_785nm_Edge_600 (500nm)_500 µm_100x_10% (3_9mW)_30 s_1 a_01.txt"
    path=r"E:\OneDrive - UVa\Univerdad Valladolid\Proyectos\Miscelanea\2023\Medidas Carina\Espectros\Resin low_785nm_Edge_600 (500nm)_500 µm_100x_10% (3_9mW)_30 s_1 a_01.txt"
    dat=Raman_dataloader.load_spectra_data(path)
    info=Raman_dataloader.load_spectra_info(path)

    x,y=select_data_region(dat[:,0],dat[:,1],[wavenumber_to_index(dat[:,0],500,1),wavenumber_to_index(dat[:,0],3500,1)])

    y_baseline=baseline_arPLS(y, lam=1E5, niter=200,full_output=False)

    y_fix=y-y_baseline
    peaks,valleys=find_peaks(x,y_fix,10,0.05)
    peak_val=[[x[peak],y_fix[peak]] for peak in peaks]
    print(peak_val)

    peak_asy=[peak_asymmetry_index(x,y_fix,peak,0.5) for peak in peaks]
    print(peak_asy)
    pos=[[x[item[0]],x[item[1]]] for item in peak_asy]
    print(pos)
    peak_type_clasification(x, y_fix,peaks)
    #model_list=['Voigt' for item in peaks]
    # fitting=Raman_fit.fit(x,y_fix,
    #             peak_val,model_list)
    # print(Raman_fit.fit_info(fitting))

    # new_y=Raman_fit.model_f(fitting.params, x,peak_val,model_list)

    #labelling of peaks
    norm_y=y_fix/max(y_fix)
    arrow_data=[(
        peaks[item],
        peak_val[item][1]
        ) for item in range(len(peaks))
        ]
    Raman_plot.plotter([[x,norm_y]],#,[x,new_y]],
                 [
                  "Raman shift"+"("+info['AxisUnit[1]']+")",
                  info['AxisType[0]']+" ("+info['AxisUnit[0]']+")"
                  ],
                  info['Title'],
                  leyends= ['Raw data'],
                   lines=True,
                   res=300,
                   size="double_size",
                   leyend_frame=[True,'b'],
                   #arrow=arrow_data
                   )
#test4()

