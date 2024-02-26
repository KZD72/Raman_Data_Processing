# -*- coding: utf-8 -*-
"""
Created on Fri May 19 13:33:04 2023

Module to create a lateral panel in the single peak GUI to analyse the data

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

# Raman_lateralpanel
import numpy as np
import tkinter as tk
from tkinter import ttk, filedialog
from tkinter import messagebox
import os
from tqdm import tqdm
import time

from bin import Raman_dataloader
from bin import Raman_datahandling
from bin import Raman_plot
from bin import Raman_single_peak_fit_GUI  # Import the FDTR_plot module for data plotting
from bin import Raman_single_peak_fit_dashboard


def error(type):
    """
    Display an error message based on the given type.

    Args:
        type (int): The error type.

    Raises:
        tkinter.TclError: If the type is not recognized.

    """
    if type == 0:
        messagebox.showerror(
            "Error", "The field must be within the spectral range")
    elif type == 1:
        messagebox.showerror(
            "Error", "The upper limit must be above the lower limit")
    elif type == 2:
        messagebox.showerror("Error", "The spectral window is too small")
    elif type == 3:
        messagebox.showerror("Error", "Range from 10e-9 to 10e9")
    elif type == 4:
        messagebox.showerror(
            "Error", "Must be an integer bigger than 3 and smaller than the length of the data")
    elif type == 5:
        messagebox.showerror("Error", "Not a number detected")
    elif type == 6:
        messagebox.showerror("Error", "Range from 10e-9 to 10e9")
    elif type == 7:
        messagebox.showerror("Error", "Range from 1 to 25")
    else:
        messagebox.showerror("Error", "Invalid data in field")
       #raise ttk.TclError("Unknown error type")


def gen_error(text):
    """
    Display an error message based on the given text.

    Args:
        text (String): The error message.

    Raises:
        tkinter.TclError: If the type is not recognized.

    """

    messagebox.showerror("Error", text)
    #raise ttk.TclError("Unknown error type")


def check_spectral_window(initial_window, new_window):
    """
    Check if the new spectral window is valid based on the initial window.

    Args:
        initial_window (list): The initial spectral window [lower_limit, upper_limit].
        new_window (list): The new spectral window [lower_limit, upper_limit].

    Returns:
        bool: True if the new window is valid, False otherwise.

    """
    check_0 = check_number_f(new_window[0]) and check_number_f(new_window[1])

    if check_0:
        new_window = [float(new_window[0]), float(new_window[1])]
        check_1 = new_window[0] >= initial_window[0] and new_window[1] <= initial_window[1]
        check_2 = new_window[0] < new_window[1]
        check_3 = new_window[0] < new_window[1] + 10
        if check_1 and check_2 and check_3:
            return True
        elif not check_1:
            messagebox.showerror(
                "Error", f"The field must be within the spectral range [{initial_window[0]}-{initial_window[1]}]")
            return False
        elif check_1 and not check_2:
            error(1)
            return False
        elif check_1 and check_2 and not check_3:
            error(2)
            return False
    else:
        error(5)
        return False


def check_baseline(n):
    """
    Check if the baseline number is valid.

    Args:
        n (float). The bvalue passed to the fucntion

    Returns:
        bool: True if the n is valid, False otherwise.

    """
    check_0 = check_number_f(n)
    if check_0:
        n = float(n)
        check_1 = n >= 1e-9 and n <= 1e9
        if check_1:
            return True
        else:
            error(3)
            return False
    else:
        error(5)
        return False


def check_man_baseline(n, window):

    check_0 = check_number_f(n)
    if check_0:
        n = float(n)
        check_1 = n >= window[0] and n <= window[1]
        if check_1:
            return True
        else:
            error(3)
            return False
    else:
        error(5)
        return False


def check_smoothing(n, l):
    """
    Check if the smoothing number is valid.

    Args:
        n (int). The bvalue passed to the function
        l (int). The legth of the data

    Returns:
        bool: True if the n is valid, False otherwise.

    """
    check_0 = check_number_i(n)
    if check_0:
        n = int(n)
        check_1 = n >= 3 and n <= l
        if check_1:
            return True
        else:
            error(4)
            return False
    else:
        error(4)
        return False


def check_peak(n):
    """
    Check if the peak number is valid.

    Args:
        n (float). The bvalue passed to the function

    Returns:
        bool: True if the n is valid, False otherwise.

    """
    check_0 = check_number_f(n)
    if check_0:
        n = float(n)
        check_1 = n > 0 and n < 1
        if check_1:
            return True
        else:
            error(6)
            return False
    else:
        error(5)
        return False


def check_n_peak(n):
    """
    Check if the peak number is valid.

    Args:
        n (float). The bvalue passed to the function

    Returns:
        bool: True if the n is valid, False otherwise.

    """
    check_0 = check_number_f(n)
    if check_0:
        n = float(n)
        check_1 = n > 0 and n < 25
        if check_1:
            return True
        else:
            error(7)
            return False
    else:
        error(5)
        return False


def check_number_f(n):
    """
    Checks if a given value can be converted to a float.

    Args:
        n (str or numeric): The value to be checked.

    Returns:
        bool: True if the value can be converted to a float, False otherwise.
    """
    try:
        float_value = float(n)  # Attempt to convert the value to a float
        return True
    except ValueError:
        return False


def check_number_i(n):
    """
    Checks if a given value can be converted to an integer.

    Args:
        n (str or numeric): The value to be checked.

    Returns:
        bool: True if the value can be converted to an integer, False otherwise.
    """
    try:
        float_value = int(n)  # Attempt to convert the value to an integer
        return True
    except ValueError:
        return False


def data_clipper(dat, spectral_window):
    """
    Clip the data based on the given spectral window.

    Args:
        dat (ndarray): The data array.
        spectral_window (list): The spectral window [lower_limit, upper_limit].

    Returns:
        tuple: A tuple containing the clipped x and y arrays.

    """
    x, y = Raman_datahandling.select_data_region(dat[:, 0], dat[:, 1], [
        Raman_datahandling.wavenumber_to_index(
            dat[:, 0], spectral_window[0], 1),
        Raman_datahandling.wavenumber_to_index(
            dat[:, 0], spectral_window[1], 1)
    ])
    return x, y


def calc_baseline(dat, spectral_window, lam=1e5):
    """
    Calculate the baseline of the data within the spectral window.

    Args:
        dat (ndarray): The data array.
        spectral_window (list): The spectral window [lower_limit, upper_limit].
        lam (float): The lambda value for baseline calculation. Default is 1e5.

    Returns:
        ndarray: The baseline-corrected y array.

    """
    x, y = data_clipper(dat, spectral_window)
    y_baseline = Raman_datahandling.baseline_arPLS(
        y, lam=lam, niter=1000, full_output=False)
    return x, y_baseline


def on_popup_close(popup):
    popup.grab_release()  # Release the grab
    popup.destroy()


def create_lateral_panel(canvas, canvas_panel, main_window, path, figure, raw_dat_a, info_a, data_type):

    # Global variables to the panel
    global raw_dat, x_raw, y_raw, x_baseline, y_baseline, spectral_window, info
    global field_low_k, field_upp_k, initial_spectral_window, lamb_value
    global peak_value, peak_val, smoothing_window, peak_sep, normalize_var, peak_number_control
    global bas_man_points, point_baseline
    global valid_files, normalise_check
    global baseline_type
    global dict_container
    global new_peak_entry
    global peak_finding_tab

    raw_dat = raw_dat_a
    x_raw = raw_dat[:, 0]
    y_raw = raw_dat[:, 1]
    x_baseline = x_raw
    y_baseline = y_raw
    initial_spectral_window = [x_raw[1], x_raw[-2]]  # intial spectral window
    spectral_window = [x_raw[1], x_raw[-2]]  # intial spectral window
    info = info_a
    lamb_value = 1  # initial value for baseline algorithm
    peak_value = 0.01
    peak_val = []
    peaks = []
    bas_man_points = np.zeros((2, 2))
    smoothing_window = 10
    peak_sep = 4
    peak_promi = 0.05
    peak_number_control = 8
    valid_files = []
    normalise_check = False
    baseline_type = "Auto"
    dict_container={}

    # Internal functions to clickbutton

    def button_clipper_clicked(silent=True):
        global info, field_low_k, field_upp_k, raw_dat, initial_spectral_window
        c1 = check_spectral_window(initial_spectral_window, [
                                   field_low_k.get(), field_upp_k.get()])
        if c1:
            spectral_window[0] = float(field_low_k.get())
            spectral_window[1] = float(field_upp_k.get())

            raw_x, raw_y = data_clipper(raw_dat, spectral_window)
            dat = [[raw_x, raw_y]]
            if silent:
                fig, ax = Raman_plot.plotter(dat,
                                             ["Wavenumber (1/cm)",
                                              "Intensity (A.U.)"],
                                             info['Title'],
                                             leyends=['Raw data'],
                                             lines=True,
                                             res=150,
                                             # size="double_size_double_heigh",
                                             leyend_frame=[True, 'b'],
                                             )
                Raman_plot.update_plot(canvas, canvas_panel, fig, ax, dat)

                button_substract_baseline.config(state="disabled")
                button_peak_detection.config(state="disabled")
                button_manual_peak_adding.config(state="disabled")
                button_peak_processing.config(state="disabled")
                button_load_batch_folder.config(state="disabled")
                button_batch.config(state="disabled")

    def button_baseline_clicked(silent=True):
        global x_baseline, y_baseline, info, field_low_k, field_upp_k, raw_dat, initial_spectral_window, lamb_value, baseline_type
    # Call the metadata function when the button is clicked
        c1 = check_spectral_window(initial_spectral_window, [
                                   field_low_k.get(), field_upp_k.get()])
        c2 = check_baseline(field_lambda.get())

        if c1 and c2:
            spectral_window[0] = float(field_low_k.get())
            spectral_window[1] = float(field_upp_k.get())
            lamb_value = float(field_lambda.get())

            raw_x, raw_y = data_clipper(raw_dat, spectral_window)
            x, y = calc_baseline(raw_dat, spectral_window, lam=1e4*lamb_value)
            dat = [[raw_x, raw_y], [x, y]]
            if silent:
                fig, ax = Raman_plot.plotter(dat,
                                             ["Wavenumber (1/cm)",
                                              "Intensity (A.U.)"],
                                             info['Title'],
                                             leyends=['Raw data', 'Baseline'],
                                             lines=True,
                                             res=150,
                                             # size="double_size_double_heigh",
                                             leyend_frame=[True, 'b'],
                                             )
                Raman_plot.update_plot(canvas, canvas_panel, fig, ax, dat)
            x_baseline = x
            y_baseline = raw_y-y
            baseline_type = "Auto"
            button_substract_baseline.config(state="normal")
            button_peak_detection.config(state="disabled")
            button_manual_peak_adding.config(state="disabled")
            button_peak_processing.config(state="disabled")
            button_load_batch_folder.config(state="disabled")
            button_batch.config(state="disabled")
            button_dashboard.config(state="disabled")

    def update_baseline(raw_x, raw_y, bas_man_points, model_type, silent=True):
        global x_baseline, y_baseline
        # guardar en algun lado
        if len(bas_man_points) > 5:
            y = Raman_datahandling.lin_baseline(
                bas_man_points, raw_x, model=model_type)
        elif len(bas_man_points) > 4:
            if model_type != "cubic":
                y = Raman_datahandling.lin_baseline(
                    bas_man_points, raw_x, model=model_type)
            else:
                y = Raman_datahandling.lin_baseline(
                    bas_man_points, raw_x, model='linear')
        elif len(bas_man_points) > 3:
            if model_type != "cubic" or model_type != "quadratic":
                y = Raman_datahandling.lin_baseline(
                    bas_man_points, raw_x, model=model_type)
            else:
                y = Raman_datahandling.lin_baseline(
                    bas_man_points, raw_x, model='linear')
        else:
            if model_type == "linear" or model_type != "slinear":
                y = Raman_datahandling.lin_baseline(
                    bas_man_points, raw_x, model=model_type)
            else:
                y = Raman_datahandling.lin_baseline(
                    bas_man_points, raw_x, model='linear')
        x = raw_x

        # get arrows:

        point_index_list = [np.where(raw_x == bas_man_points[item, 0])[0][0] for item in range(
            len(bas_man_points[:, 0])) if len(np.where(raw_x == bas_man_points[item, 0])[0]) > 0]

        arrow_data = [(
            point_index_list[item],
            y[point_index_list[item]]
        ) for item in range(len(point_index_list))
        ]

        dat = [[raw_x, raw_y], [x, y]]
        if silent:
            fig, ax = Raman_plot.plotter(dat,
                                         ["Wavenumber (1/cm)",
                                          "Intensity (A.U.)"],
                                         info['Title'],
                                         leyends=['Raw data', 'Baseline'],
                                         lines=True,
                                         res=150,
                                         arrow=arrow_data,
                                         # size="double_size_double_heigh",
                                         leyend_frame=[True, 'b'],
                                         )
            Raman_plot.update_plot(canvas, canvas_panel, fig, ax, dat)
        x_baseline = x
        y_baseline = raw_y-y
        point_baseline.delete(0, tk.END)  # Clear the current content
        point_baseline.insert(0, 0.0)  # Clear the current content

    def button_man_baseline_clicked(silent=True):
        global x_baseline, y_baseline, bas_man_points, baseline_type
        global info, field_low_k, field_upp_k, raw_dat, initial_spectral_window, point_baseline
    # Call the metadata function when the button is clicked
        c1 = check_spectral_window(initial_spectral_window, [
                                   field_low_k.get(), field_upp_k.get()])
        c2 = check_man_baseline(point_baseline.get(), [x_raw[0], x_raw[-1]])
        model_type = str(baseline_model.get())

        if c1 and c2:
            spectral_window[0] = float(field_low_k.get())
            spectral_window[1] = float(field_upp_k.get())
            point = float(point_baseline.get())
            if point >= spectral_window[0] and point <= spectral_window[1]:

                raw_x, raw_y = data_clipper(raw_dat, spectral_window)
                bas_man_points[0] = [raw_x[0], raw_y[0]]
                bas_man_points[-1] = [raw_x[-1], raw_y[-1]]

                # Find the index where the point should be inserted based on x

                res = np.abs((raw_x[0]-raw_x[1])/2)
                index = Raman_datahandling.wavenumber_to_index(
                    raw_x, point, res)

                insert_point = [raw_x[index], raw_y[index]]

                insert_index = np.searchsorted(
                    bas_man_points[:, 0], insert_point[0])
                bas_man_points = np.insert(
                    bas_man_points, insert_index, insert_point, axis=0)

                baseline_type = "Manual"

                update_baseline(raw_x, raw_y, bas_man_points,
                                model_type, silent)
                button_substract_baseline.config(state="normal")
                button_peak_detection.config(state="disabled")
                button_manual_peak_adding.config(state="disabled")
                button_peak_processing.config(state="disabled")
                button_load_batch_folder.config(state="disabled")
                button_batch.config(state="disabled")
                button_dashboard.config(state="disabled")
            else:
                gen_error(
                    f"Point out of data range [{spectral_window[0]}-{spectral_window[1]}]")

    def button_man_baseline_reset():
        global x_baseline, y_baseline, bas_man_points
        global info, field_low_k, field_upp_k, raw_dat, initial_spectral_window, point_baseline
        c1 = check_spectral_window(initial_spectral_window, [
                                   field_low_k.get(), field_upp_k.get()])
        point_baseline.delete(0, tk.END)  # Clear the current content
        point_baseline.insert(0, 0.0)  # Clear the current content

        if c1:
            spectral_window[0] = float(field_low_k.get())
            spectral_window[1] = float(field_upp_k.get())

            raw_x, raw_y = data_clipper(raw_dat, spectral_window)
            bas_man_points[0] = [raw_dat[0, 0], raw_dat[0, 1]]
            bas_man_points[-1] = [raw_dat[-1, 0], raw_dat[-1, 1]]

            raw_x, raw_y = data_clipper(raw_dat, spectral_window)

            bas_man_points = np.zeros((2, 2))
            bas_man_points[0] = [raw_dat[0, 0], raw_dat[0, 1]]
            bas_man_points[-1] = [raw_dat[-1, 0], raw_dat[-1, 1]]

            update_baseline(raw_x, raw_y, bas_man_points, 'linear')
            button_substract_baseline.config(state="disabled")
            button_peak_detection.config(state="disabled")
            button_manual_peak_adding.config(state="disabled")
            button_peak_processing.config(state="disabled")
            button_load_batch_folder.config(state="disabled")
            button_batch.config(state="disabled")
            button_dashboard.config(state="disabled")
            

    def button_add_baseline():
        # Prompt user to select a file from disk
        global x_baseline, y_baseline, bas_man_points
        global info, field_low_k, field_upp_k, raw_dat, initial_spectral_window, point_baseline

        model_type = str(baseline_model.get())

        filepath = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv"), (
            "TXT Files", "*.txt")])
        if filepath != '':
            data = np.genfromtxt(filepath, dtype=float, delimiter=' ')

        c1 = check_spectral_window(initial_spectral_window, [
                                   field_low_k.get(), field_upp_k.get()])
        if isinstance(data, np.ndarray):
            c2 = np.array([not (str(val).isalpha() or np.isnan(val)
                                or np.isinf(val)) for val in data])
        else:
            data = np.array([data])
            c2 = np.array([not (str(val).isalpha() or np.isnan(val)
                                or np.isinf(val)) for val in data])

        # Missing all error control
        if c1 and c2.all():
            c3 = np.min(data) >= float(field_low_k.get()) and np.max(
                data) <= float(field_upp_k.get())

            if c3:
                raw_x, raw_y = data_clipper(raw_dat, spectral_window)
                bas_man_points = np.zeros((2, 2))
                bas_man_points[0] = [raw_x[0], raw_y[0]]
                bas_man_points[-1] = [raw_x[-1], raw_y[-1]]

                res = np.abs((raw_x[0]-raw_x[1])/2)
                index = [Raman_datahandling.wavenumber_to_index(
                    raw_x, point, res) for point in data]

                for item in index:
                    insert_point = [raw_x[item], raw_y[item]]
                    insert_index = np.searchsorted(
                        bas_man_points[:, 0], insert_point[0])
                    bas_man_points = np.insert(
                        bas_man_points, insert_index, insert_point, axis=0)

                update_baseline(raw_x, raw_y, bas_man_points, model_type)
            else:
                gen_error(
                    "Check the range of the provided data points, it is out of range")
        else:
            gen_error(
                "There are non numerical values, or wrong format in the provided file")

        button_substract_baseline.config(state="normal")
        button_peak_detection.config(state="disabled")
        button_manual_peak_adding.config(state="disabled")
        button_peak_processing.config(state="disabled")
        button_load_batch_folder.config(state="disabled")
        button_batch.config(state="disabled")
        button_dashboard.config(state="disabled")

    def button_save_baseline():
        global x_baseline, y_baseline, bas_man_points
        global info, field_low_k, field_upp_k, raw_dat, initial_spectral_window, point_baseline
        filepath = filedialog.asksaveasfilename(filetypes=[("CSV Files", "*.csv"), (
            "TXT Files", "*.txt")])
        if filepath != '':
            if len(bas_man_points[:, 0]) > 2:
                # Save the vector to a text file
                np.savetxt(
                    filepath, bas_man_points[1:-1, 0], fmt="%d", delimiter=" ")
                print(bas_man_points[1:-1, 0])
            else:
                np.savetxt(filepath, [], fmt="%d", delimiter=" ")

    def baseline_model_selection(event=None):
        global x_baseline, y_baseline, bas_man_points
        global info, field_low_k, field_upp_k, raw_dat, initial_spectral_window, point_baseline
        model_type = str(baseline_model.get())
        c1 = check_spectral_window(initial_spectral_window, [
                                   field_low_k.get(), field_upp_k.get()])
        c2 = len(bas_man_points) > 2
        if c1 and c2:
            raw_x, raw_y = data_clipper(raw_dat, spectral_window)
            update_baseline(raw_x, raw_y, bas_man_points, model_type)

    def button_baseline_removed_clicked(silent=True):
        global x_baseline, y_baseline, info

        dat = [[x_baseline, y_baseline]]
        if silent:
            fig, ax = Raman_plot.plotter(dat,
                                         ["Wavenumber (1/cm)",
                                          "Intensity (A.U.)"],
                                         info['Title'],
                                         leyends=['Baseline_corrected'],
                                         lines=True,
                                         res=150,
                                         # size="double_size_double_heigh",
                                         leyend_frame=[True, 'b'],
                                         )
            Raman_plot.update_plot(canvas, canvas_panel, fig, ax, dat)
        button_peak_detection.config(state="normal")
        button_manual_peak_adding.config(state="normal")
        button_peak_processing.config(state="disabled")
        button_normalise.config(state="normal")
        button_load_batch_folder.config(state="disabled")
        button_batch.config(state="disabled")
        button_dashboard.config(state="disabled")

    def button_normalise_clicked(silent=True):
        global x_baseline, y_baseline, info, normalise_check
        normalise_check = True

        y_baseline = y_baseline/np.max(y_baseline)
        dat = [[x_baseline, y_baseline]]
        if silent:
            fig, ax = Raman_plot.plotter(dat,
                                         ["Wavenumber (1/cm)",
                                          "Intensity (A.U.)"],
                                         info['Title'],
                                         leyends=[
                                             'Normalised baseline corrected'],
                                         lines=True,
                                         res=150,
                                         # size="double_size_double_heigh",
                                         leyend_frame=[True, 'b'],
                                         )
            Raman_plot.update_plot(canvas, canvas_panel, fig, ax, dat)
        button_peak_detection.config(state="normal")
        button_peak_processing.config(state="disabled")
        button_load_batch_folder.config(state="disabled")
        button_batch.config(state="disabled")

    def button_peak_detection_clicked(silent=True):
        global x_baseline, y_baseline, info, peaks, peak_value, peak_val, smoothing_window, peak_sep, peak_promi, peak_number_control

       
        c1 = check_peak(field_peaks.get())
        c2 = check_smoothing(field_smooth_window.get(), len(x_baseline))
        c3 = check_smoothing(field_sep.get(), len(x_baseline))
        c4 = check_peak(field_promi.get())
        c5 = check_n_peak(field_peak_number_control.get())
        if c1 and c2 and c3 and c4 and c5:
            peak_value = float(field_peaks.get())
            smoothing_window = int(field_smooth_window.get())
            peak_sep = float(field_sep.get())
            peak_promi = float(field_promi.get())
            peak_number_control = float(field_peak_number_control.get())

            peaks, valleys = Raman_datahandling.find_peaks(x_baseline, y_baseline,
                                                           smoothing_window,
                                                           peak_value,
                                                           prominence=peak_promi,
                                                           distance=peak_sep)

            peak_mag = [y_baseline[peak] for peak in peaks]
            sorted_indices = np.argsort(peak_mag)[::-1]
            peaks = [peaks[item] for item in sorted_indices]
            peak_val = [[x_baseline[peaks[item]], y_baseline[peaks[item]]]
                        for item in range(len(peaks)) if item < peak_number_control]
            peaks = [peaks[item]
                     for item in range(len(peaks)) if item < peak_number_control]
            arrow_data = [(
                peaks[item],
                peak_val[item][1]
            ) for item in range(len(peaks))
            ]
            dat = [[x_baseline, y_baseline]
                   ]
            if silent:
                fig, ax = Raman_plot.plotter(dat,
                                             ["Wavenumber (1/cm)",
                                              "Intensity (A.U.)"],
                                             info['Title'],
                                             leyends=['Baseline_corrected'],
                                             lines=True,
                                             res=150,
                                             # size="double_size_double_heigh",
                                             leyend_frame=[True, 'b'],
                                             arrow=arrow_data
                                             )
                Raman_plot.update_plot(canvas, canvas_panel, fig, ax, dat)
            button_peak_processing.config(state="normal")
            button_peak_adding.config(state="normal")
            button_load_batch_folder.config(state="normal")
            button_batch.config(state="disabled")
    
    
    def add_peak(new_peak):
        """
        This function adds a new peak to the global peak_val list and updates the plot.

        Parameters:
        new_peak: The new peak to be added. It should be a Tkinter StringVar object that contains a string representation of a float.

        Returns:
        None
        """
        # Access the global variables
        global x_baseline, y_baseline, info, peak_val, peaks

        # Check if the new peak is a number
        c1 = check_number_f(new_peak.get())

        if c1:
            # Convert the new peak to a float
            p = float(new_peak.get())

            # Check if the new peak is within the range of x_baseline
            c2 = p > x_baseline[0] and p < x_baseline[-1]

            if c2:
                # Convert peak_val to a NumPy array
                np_peaks = np.asarray(peak_val)

                # Check if np_peaks is not empty
                if np_peaks.size:
                    # Get the x values of the peaks
                    x_peak = np_peaks[:, 0]

                    # Check if the new peak is already in x_peak
                    c3 = np.isin(p, x_peak)
                else:
                    c3 = False

                # If the new peak is not already in x_peak
                if not c3:
                    # Find the index of the new peak in x_baseline
                    index = np.searchsorted(np.asarray(x_baseline), p)

                    # Add the new peak to peak_val and peaks
                    peak_val.append([p, y_baseline[index]])
                    peaks.append(index)

                    # Reorder the peaks by y value in descending order
                    sorted_indexes = sorted(range(len(peak_val)),
                                            key=lambda i: peak_val[i][1],
                                            reverse=True)
                    peak_val = [peak_val[item] for item in sorted_indexes]
                    peaks = [peaks[item] for item in sorted_indexes]

                    # Update the plot
                    arrow_data = [(peaks[item], peak_val[item][1]) for item in range(len(peak_val))]
                    dat = [[x_baseline, y_baseline]]
                    fig, ax = Raman_plot.plotter(dat,
                                                    ["Wavenumber (1/cm)", "Intensity (A.U.)"],
                                                    info['Title'],
                                                    leyends=['Baseline_corrected'],
                                                    lines=True,
                                                    res=150,
                                                    leyend_frame=[True, 'b'],
                                                    arrow=arrow_data)
                    Raman_plot.update_plot(canvas, canvas_panel, fig, ax, dat)

                else:
                    error(7)  # The new peak is already in peak_val
            else:
                error(7)  # The new peak is not within the range of x_baseline
        else:
            error(7)  # The new peak is not a number

    def button_adding_clicked():
        # Access the global variables
        global x_baseline, y_baseline, info, peak_val, peaks

        # deactivate buttons
        button_clipper.config(state="disabled")
        button_baseline.config(state="disabled")
        button_substract_baseline.config(state="disabled")
        button_peak_detection.config(state="disabled")
        button_peak_processing.config(state="disabled")
        button_peak_adding.config(state="disabled")
        button_normalise.config(state="disabled")
        button_batch.config(state="disabled")
        button_dashboard.config(state="disabled")

        def on_closing():
            button_clipper.config(state="normal")
            button_substract_baseline.config(state="normal")
            button_baseline.config(state="normal")
            button_peak_detection.config(state="normal")
            button_peak_processing.config(state="normal")
            button_peak_adding.config(state="normal")
            button_normalise.config(state="normal")
            secondary_window.quit()
            secondary_window.destroy()

        def add_peak_local():
            add_peak(new_peak)

        # Area to create the peak fittingwindow
        secondary_window = tk.Toplevel(main_window)
        secondary_window.title('Add non-detected peak')
        secondary_window.geometry("400x150")
        secondary_window.resizable(False, False)  # Disable resizing
        secondary_window.grab_set()  # Grab the focus for the popup window
        secondary_window.protocol(
            "WM_DELETE_WINDOW", lambda: on_popup_close(secondary_window))

        secondary_panel = tk.Frame(secondary_window, bg='white')
        # Adjust the row and column values as needed
        secondary_panel.grid(row=0, column=0, padx=10, pady=2)

        # Create a style for the labels
        style = ttk.Style()
        style.configure(
            "Box.TLabel", background=secondary_window["background"])

        box_frame = ttk.Frame(secondary_panel, borderwidth=1, relief="groove")
        box_frame.grid(row=1, column=0, padx=10, pady=2, sticky='nsew')

        field_label_add_peak = ttk.Label(
            box_frame,
            text="Peak position in 1/cm:",
            anchor="center",
            justify="center",
            style="Box.TLabel"
        )

        field_label_add_peak.grid(
            row=0, column=0, padx=5, pady=2, sticky='nsew')
        field_new_peak = ttk.Entry(box_frame, justify='center')
        new_peak = field_new_peak
        field_new_peak.grid(row=0, column=1, padx=5, pady=2, sticky='nsew')

        # add button to add peak
        # Create button
        button_peak_add = tk.Button(box_frame, text='Add peak',
                                    command=add_peak_local
                                    )
        button_peak_add.grid(row=0, column=2, padx=10, pady=2, sticky='nsew')
        # Create button
        button_peak_close = tk.Button(box_frame, text='Close window',
                                      command=on_closing)
        button_peak_close.grid(row=1, column=0, padx=10, pady=2, sticky='nsew')
        secondary_window.protocol("WM_DELETE_WINDOW", on_closing)
        secondary_window.mainloop()

    def button_peak_analizer_clicked():
        global x_baseline, y_baseline, info, peak_value, peak_val

        def on_closing():
            button_clipper.config(state="normal")
            button_substract_baseline.config(state="normal")
            button_baseline.config(state="normal")
            button_peak_detection.config(state="normal")
            button_peak_processing.config(state="normal")
            button_peak_adding.config(state="normal")
            button_normalise.config(state="normal")
            button_load_batch_folder.config(state="normal")
            button_batch.config(state="disabled")
            fit_window.quit()
            fit_window.destroy()
        # deactivate buttons
        button_clipper.config(state="disabled")
        button_baseline.config(state="disabled")
        button_substract_baseline.config(state="disabled")
        button_peak_detection.config(state="disabled")
        button_manual_peak_adding.config(state="disabled")
        button_peak_processing.config(state="disabled")
        button_peak_adding.config(state="disabled")
        button_normalise.config(state="disabled")
        button_load_batch_folder.config(state="disabled")
        button_batch.config(state="disabled")
        button_dashboard.config(state="disabled")


        # Create window for fitting:
        # Area to create the peak fittingwindow
        fit_window = tk.Toplevel(main_window)
        fit_window.title('Raman peak analizer')
        fit_window.geometry("755x900")
        fit_window.resizable(False, False)  # Disable resizing
        fit_window.attributes("-topmost", True)
        # Grid layout configuration
        fit_window.grid_columnconfigure(0, weight=1)
        fit_window.grid_rowconfigure(0, weight=1)
        Raman_single_peak_fit_GUI.create_fit_panel(
            fit_window, canvas, canvas_panel, info, x_baseline, y_baseline, peak_val)

        fit_window.protocol("WM_DELETE_WINDOW", on_closing)
        fit_window.mainloop()

    def button_manual_adding_clicked():
        global peak_val, peaks, new_peak_entry  
        button_load_batch_folder.config(state="normal")
        button_peak_processing.config(state="normal")
        add_peak(new_peak_entry)

    def clean_peaks(event):
        global peak_val, peaks, x_baseline, y_baseline, info, selected_tab  

        if len(peak_val)>0:
            dat = [[x_baseline, y_baseline]]
            fig, ax = Raman_plot.plotter(dat,
                                            ["Wavenumber (1/cm)",
                                            "Intensity (A.U.)"],
                                            info['Title'],
                                            leyends=[
                                                'Baseline_corrected'],
                                            lines=True,
                                            res=150,
                                            # size="double_size_double_heigh",
                                            leyend_frame=[True, 'b'],                                        
                                            )
            Raman_plot.update_plot(
                canvas, canvas_panel, fig, ax, dat)
        peak_val=[]
        peaks=[]
        

    def load_batch_folder():
        """
        Opens a file dialog to select and load the data file.
        """
        global valid_files

        # Open the dialog to select a directory
        directory = filedialog.askdirectory()

        # Get a list of all files in the directory
        files = os.listdir(directory)

        # Filter the list to only include .txt or datfiles
        txt_files = [f for f in files if f.endswith('.txt') or f.endswith('.dat')]

        # Create a new Tkinter window
        root = tk.Tk()
        root.title("Checking files Progress")
           # Create a progress bar
          # Create a label for the estimated time remaining
        time_label = tk.Label(root, text="")
        time_label.pack()
        # Create a style for the progress bar
        style = ttk.Style()
        style.configure("TProgressbar", thickness=50)  # Adjust the thickness as needed
        progress = ttk.Progressbar(root, length=250, mode='determinate', style="TProgressbar") 
        progress.pack()

        # Valid data files
        valid_files = []
        numer_non_valid_files = 0
        for item, file in enumerate(txt_files):
            start_time = time.time()
            try:
                check, key = Raman_dataloader.load_spectra_data(
                    os.path.join(directory, file), data_type, silent=False)
                if key:
                    valid_files.append(os.path.join(directory, file))
                else:
                    numer_non_valid_files = numer_non_valid_files+1
            except:
                numer_non_valid_files = numer_non_valid_files+1
            
            progress['value'] = (item+1) / len(txt_files) * 100
            elapsed_time = time.time() - start_time
            estimated_time = elapsed_time * (len(txt_files) - item)
            
            time_label['text'] = "{:.2f} % completed".format((item+1) / len(txt_files) * 100)+"\nEstimated time remaining:\n{:.2f} seconds".format(estimated_time)
            root.update()
            root.update_idletasks()

        # Close the Tkinter window once the loop is finished
        root.destroy()
        messagebox.showinfo("Information", f"{len(valid_files)} files valid to be processed,{numer_non_valid_files} files to be excluded", parent=main_window)

        if len(valid_files) > 0:
            button_batch.config(state="normal")

    def button_batch_proccesing():
        global valid_files, normalise_check, raw_dat, info 
        global bas_man_points, x_baseline, y_baseline,dict_container
        global peak_finding_tab, peak_val, peaks

        new_x_baseline = [bas_man_points[i, 0]
                          for i in range(1, len(bas_man_points) - 1)]
        
        
       
        # Create a new Tkinter window
        root2 = tk.Tk()
        root2.title("Processing Progress")
           # Create a progress bar
          # Create a label for the estimated time remaining
        time_label = tk.Label(root2, text="")
        time_label.pack()
        # Create a style for the progress bar
        style = ttk.Style()
        style.configure("TProgressbar", thickness=50)  # Adjust the thickness as needed
        progress = ttk.Progressbar(root2, length=250, mode='determinate', style="TProgressbar") 
        progress.pack()

        dict_container={}

        for item,file in enumerate(valid_files):
            
            raw_dat, key = Raman_dataloader.load_spectra_data(file, data_type)
            info, key2 = Raman_dataloader.load_spectra_info(file, data_type)
            plotShow = False  # False To not show the images
            start_time = time.time()
            if file == valid_files[-1]:
                plotShow = False
            try:
                button_clipper_clicked(silent=False)
                try:
                    if baseline_type == "Auto":
                        button_baseline_clicked(silent=False)
                    else:
                        raw_x, raw_y = data_clipper(raw_dat, spectral_window)
                        # Look for the datapoints in ech set of data:

                        bas_man_points = np.zeros((2, 2))
                        bas_man_points[0] = [raw_x[0], raw_y[0]]
                        bas_man_points[-1] = [raw_x[-1], raw_y[-1]]
                        res = np.abs((raw_x[0]-raw_x[1])/4)
                        index = [Raman_datahandling.wavenumber_to_index(
                            raw_x, point, res) for point in new_x_baseline]

                        for element in index:
                            insert_point = [raw_x[element], raw_y[element]]
                            insert_index = np.searchsorted(
                                bas_man_points[:, 0], insert_point[0])
                            bas_man_points = np.insert(
                                bas_man_points, insert_index, insert_point, axis=0)

                        model_type = str(baseline_model.get())

                        update_baseline(
                            raw_x, raw_y, bas_man_points, model_type, silent=False)
                    try:
                        button_baseline_removed_clicked(silent=False)
                        if normalise_check:
                            button_normalise_clicked(silent=False)
                        try:
                            selected_tab = peak_finding_tab.tab(peak_finding_tab.select(), "text")
                            if selected_tab=="Auto":
                                button_peak_detection_clicked(silent=False)
                            else:
                                 peak_val= peak_val
                                 peaks=peaks
                            try:
                                dict_container.update(Raman_single_peak_fit_GUI.batch_fit(
                                    canvas, canvas_panel, info, x_baseline, y_baseline, peak_val, file, silent=plotShow)
                                )
                                progress['value'] = (item+1) / len(valid_files) * 100                               
                                elapsed_time = time.time() - start_time                                
                                estimated_time = elapsed_time * (len(valid_files) - item)                                
                                time_label['text'] = "{:.2f} % completed".format((item+1) / len(valid_files) * 100)+"\nEstimated time remaining:\n{:.2f} seconds".format(estimated_time)
                                
                            
                                root2.update()
                                root2.update_idletasks()
                            except:
                                messagebox.showinfo(
                                    "Error", f"{file} \nFit routine failed")
                        except:
                            messagebox.showinfo(
                                "Error", f"{file} \nPeak finding failed")

                    except:
                        messagebox.showinfo(
                            "Error", f"{file} \nBaseline substraction failed")

                except:
                    messagebox.showinfo(
                        "Error", f"{file} \nManual baseline failed")

            except:
                messagebox.showinfo("Error", f"{file} \nout of clipping range")
            
            

        # Close the Tkinter window once the loop is finished
        root2.destroy()
        button_dashboard.config(state="normal")
        

    def dashboard():
        global dict_container
        dashboard_window = tk.Toplevel(main_window)
        dashboard_window.title('Batch peak analizer dahsboard')
        dashboard_window.geometry("900x600")
        dashboard_window.resizable(False, False)  # Disable resizing
        Raman_single_peak_fit_dashboard.create_dashboard(dashboard_window, canvas, canvas_panel, dict_container)
     ###########################################################################
     ##### Main lateral panel definition                                    ####
     ###########################################################################

    main_panel = tk.Frame(main_window, bg='white')
    # Adjust the row and column values as needed
    main_panel.grid(row=0, column=0, padx=10, pady=2)

    # Create subpanel frames within the main panel
    subpanel1 = tk.Frame(main_panel, bg='lightblue')
    subpanel2 = tk.Frame(main_panel, bg='lightgrey')
    subpanel3 = tk.Frame(main_panel, bg='lightpink')
    subpanel4 = tk.Frame(main_panel, bg='light salmon')

    # Grid layout configuration
    main_panel.grid_columnconfigure(0, weight=1)
    main_panel.grid_rowconfigure(0, weight=1)
    main_panel.grid_rowconfigure(1, weight=1)
    main_panel.grid_rowconfigure(1, weight=1)
    main_panel.grid_rowconfigure(2, weight=1)
    main_panel.grid_rowconfigure(2, weight=1)
    main_panel.grid_rowconfigure(3, weight=1)
    main_panel.grid_rowconfigure(3, weight=1)
    # Grid placement of subpanel frames
    subpanel1.grid(row=0, column=0, padx=10, pady=2, sticky='nsew')
    subpanel2.grid(row=1, column=0, padx=10, pady=2, sticky='nsew')
    subpanel3.grid(row=2, column=0, padx=10, pady=2, sticky='nsew')
    subpanel4.grid(row=3, column=0, padx=10, pady=2, sticky='nsew')

    # Create the decoration of the subpanels
    # Create a style for the labels
    style = ttk.Style()
    style.configure("Box.TLabel", background=main_window["background"])

    ###########################################################################
    ##### Sub panel 1                                                      ####
    ###########################################################################

    # Subpanel 1

    subpanel1.grid_rowconfigure(1, weight=1)  # Configure row weight
    subpanel1.grid_columnconfigure(0, weight=1)  # Configure column weight

    box_frame = ttk.Frame(subpanel1, borderwidth=1, relief="groove")
    box_frame.grid(row=1, column=0, padx=10, pady=2, sticky='nsew')
    
      # Configure column weights of the box_frame
    box_frame.grid_columnconfigure(0, weight=1)
    box_frame.grid_columnconfigure(1, weight=1)
    # Add elements to each subpanel

   
    # Create the tabs 
    range_tab = ttk.Notebook(box_frame)
    range_tab.grid(row=0, column=0, padx=10, pady=2, sticky='nsew')
    tab0 = ttk.Frame(range_tab)
    tab01 = ttk.Frame(range_tab)
    range_tab.add(tab0, text="Simple ROI")
    range_tab.add(tab01, text="Segmented ROI")
    
    ### TAB0
    
     # First child frame
    frame1 = ttk.Frame(tab0, borderwidth=2, relief="groove")
    frame1.grid(row=1, column=0, padx=10, pady=2)

    # Second child frame
    frame2 = ttk.Frame(tab0, borderwidth=2, relief="groove")
    frame2.grid(row=1, column=1, padx=10, pady=2)
  

    # area to create the fields
    label_fields =  tk.Label(subpanel1, text="Select the spectral range:")
    label_fields.config(font=("bold", 14))
    label_fields.config(font=("underline"))
    # Create a bold and underlined font

    label_fields.grid(row=0, column=0, padx=5, pady=2)

    field_label_low_k = ttk.Label(
        frame1,
        text="Min. wavenumber\n(1/cm)",
        anchor="center",
        justify="center",
        style="Box.TLabel"
    )

    field_label_low_k.grid(row=0, column=0, padx=5, pady=2)
    field_low_k = ttk.Entry(frame1, justify='center')
    field_low_k.insert(0, spectral_window[0])
    field_low_k.grid(row=1, column=0, padx=5, pady=2)

    field_label_upp_k = ttk.Label(
        frame2,
        text="Max. wavenumber\n(1/cm)",
        anchor="center",
        justify="center",
        style="Box.TLabel"
    )
    field_label_upp_k.grid(row=0, column=1, padx=5, pady=2)

    field_upp_k = ttk.Entry(frame2, justify='center')
    field_upp_k.insert(0, spectral_window[1])
    field_upp_k.grid(row=1, column=1, padx=5, pady=2)

    spectral_window[0] = float(field_low_k.get())
    spectral_window[1] = float(field_upp_k.get())
    check_spectral_window(initial_spectral_window, spectral_window)

    button_clipper = tk.Button(tab0, text='Clip data',
                               command=button_clipper_clicked)
    button_clipper.grid(row=1, column=2, padx=10, pady=2)
    
    
    ### TAB01
    

    ###########################################################################
    ##### Sub panel 2                                                      ####
    ###########################################################################
   # Subpanel 2
    subpanel2.grid_rowconfigure(1, weight=1)  # Configure row weight
    subpanel2.grid_columnconfigure(0, weight=1)  # Configure column weight

    box_frame_2 = ttk.Frame(subpanel2, borderwidth=1, relief="groove")
    box_frame_2.grid(row=1, column=0, padx=10, pady=2, sticky='nsew')
    # Configure column weights of the box_frame
    box_frame_2.grid_columnconfigure(0, weight=1)
    box_frame_2.grid_columnconfigure(1, weight=1)
    box_frame_2.grid_columnconfigure(3, weight=1)

    label_baseline = tk.Label(subpanel2, text='Baseline Correction:')
    label_baseline.config(font=("bold", 14))
    label_baseline.config(font=("underline"))
    label_baseline.grid(row=0, column=0, padx=10, pady=2, sticky='nsew')

    # Second child frame
    frame4 = ttk.Frame(box_frame_2, borderwidth=2, relief="groove")
    frame4.grid(row=1, column=1, padx=10, pady=2, sticky='nsew')

    # third child frame
    frame5 = ttk.Frame(box_frame_2, borderwidth=2, relief="groove")
    frame5.grid(row=1, column=2, padx=10, pady=2, sticky='nsew')

    # Buttons subpanel2
    baseline_tab = ttk.Notebook(frame4)
    baseline_tab.grid(row=0, column=0, padx=10, pady=2, sticky='nsew')
    tab1 = ttk.Frame(baseline_tab)
    tab2 = ttk.Frame(baseline_tab)
    baseline_tab.add(tab1, text="Auto")
    baseline_tab.add(tab2, text="Manual")

    # First child frame
    frame3 = ttk.Frame(tab1, borderwidth=2, relief="groove")
    frame3.grid(row=1, column=0, padx=10, pady=2, sticky='nsew')

    ###########################################################################
    ##### Auto                                                             ####
    ###########################################################################

    # Control of baseline
    field_label_bs_control = ttk.Label(
        tab1,
        text="Baseline lambda:",
        anchor="center",
        justify="center",
        style="Box.TLabel"
    )

    field_label_bs_control.grid(row=0, column=0, padx=5, pady=2, sticky='nsew')
    field_lambda = ttk.Entry(tab1, justify='center')
    field_lambda.insert(0, 1)
    lamb_value = float(field_lambda.get())
    field_lambda.grid(row=1, column=0, padx=5, pady=2, sticky='nsew')

    button_baseline = tk.Button(tab1, text='Baseline fit',
                                command=button_baseline_clicked)
    button_baseline.grid(row=1, column=1, padx=10, pady=2, sticky='nsew')

    ###########################################################################
    ##### Manual                                                           ####
    ###########################################################################

    # Load baseline
    load_baseline = tk.Button(tab2, text='Load file',
                              command=button_add_baseline)
    load_baseline.grid(row=0, column=1, padx=10, pady=2, sticky='nsew')

    save_baseline = tk.Button(tab2, text='Save points',
                              command=button_save_baseline)
    save_baseline.grid(row=0, column=2, padx=10, pady=2, sticky='nsew')

    options_bl = ['linear', 'slinear', 'quadratic', 'cubic']

    combobox_baseline = ttk.Combobox(tab2, values=options_bl, width=6)
    combobox_baseline.set(options_bl[0])
    baseline_model = combobox_baseline
    combobox_baseline.bind("<<ComboboxSelected>>", baseline_model_selection)
    combobox_baseline.grid(row=0, column=0, padx=1, pady=1)

    # Control of baseline
    point_baseline_add = ttk.Label(
        tab2,
        text="Add point (x coord):",
        anchor="center",
        justify="center",
        style="Box.TLabel"
    )

    point_baseline_add.grid(row=1, column=0, padx=5, pady=2, sticky='nsew')
    point_baseline = ttk.Entry(tab2, justify='center', width=3)
    point_baseline.insert(0, 0.0)
    point_baseline.grid(row=2, column=0, padx=5, pady=2, sticky='nsew')

    button_baseline = tk.Button(tab2, text='Add',
                                command=button_man_baseline_clicked)
    button_baseline.grid(row=2, column=1, padx=10, pady=2, sticky='nsew')

    button_baseline_res = tk.Button(tab2, text='Reset',
                                    command=button_man_baseline_reset)
    button_baseline_res.grid(row=2, column=2, padx=10, pady=2, sticky='nsew')

    ###########################################################################
    ##### Common                                                           ####
    ###########################################################################

    button_substract_baseline = tk.Button(frame5, text='Substract Baseline',
                                          command=button_baseline_removed_clicked,
                                          state="disabled")
    button_substract_baseline.grid(
        row=1, column=2, padx=10, pady=2, sticky='nsew')

    # check to normalize
    # Create a variable to hold the checkbox state

    # Create Normalize
    button_normalise = tk.Button(frame5, text='Normalise',
                                 command=button_normalise_clicked,
                                 state="disabled")
    button_normalise.grid(row=2, column=2, padx=10, pady=2, sticky='nsew')

    ###########################################################################
    ##### Sub panel 3                                                      ####
    ###########################################################################

    # Subpanel 3
    subpanel3.grid_rowconfigure(1, weight=1)  # Configure row weight
    subpanel3.grid_columnconfigure(0, weight=1)  # Configure column weight

    box_frame_3 = ttk.Frame(subpanel3, borderwidth=1, relief="groove")
    box_frame_3.grid(row=1, column=0, padx=10, pady=2, sticky='nsew')
    # Configure column weights of the box_frame
    box_frame_3.grid_columnconfigure(0, weight=4)
   

           # Buttons subpanel2
    peak_finding_tab = ttk.Notebook(box_frame_3)
    peak_finding_tab.grid(row=0, column=0, padx=10, pady=2, sticky='nsew')
    tab3 = ttk.Frame(peak_finding_tab)
    tab4 = ttk.Frame(peak_finding_tab)
    peak_finding_tab.add(tab3, text="Auto")
    peak_finding_tab.add(tab4, text="Manual")
    peak_finding_tab.bind("<<NotebookTabChanged>>", clean_peaks)
    tab4.grid_columnconfigure(0, weight=4)
    tab4.grid_columnconfigure(1, weight=2)
    tab4.grid_columnconfigure(3, weight=6)
    tab3.grid_columnconfigure(0, weight=4)

    ###########################################################################
    ##### Auto                                                             ####
    ###########################################################################
    # First child frame
    frame6 = ttk.Frame(tab3, borderwidth=2, relief="groove")
    frame6.grid(row=0, column=0, padx=10, pady=2, sticky='nsew')
    
    # Second child frame
    frame7 = ttk.Frame(tab3, borderwidth=2, relief="groove")
    frame7.grid(row=0, column=1, padx=10, pady=2, sticky='nsew')

    # third child frame
    frame8 = ttk.Frame(tab3, borderwidth=2, relief="groove")    
    frame8.grid(row=1, column=1, padx=10, pady=2, sticky='nsew')
    # Buttons subpanel2

    label_peaks = tk.Label(subpanel3, text='Peak detection:')
    label_peaks.config(font=("bold", 14))
    label_peaks.config(font=("underline"))
    label_peaks.grid(row=0, column=0, padx=10, pady=2, sticky='nsew')

    # Control of smoothing window
    field_label_smooth_control = ttk.Label(
        frame7,
        text="Smoothing window:",
        anchor="center",
        justify="center",
        style="Box.TLabel"
    )

    field_label_smooth_control.grid(
        row=0, column=0, padx=5, pady=2, sticky='nsew')
    field_smooth_window = ttk.Entry(frame7, justify='center')
    field_smooth_window .insert(0, smoothing_window)
    smoothing_window = int(field_smooth_window.get())
    field_smooth_window .grid(row=1, column=0, padx=5, pady=2, sticky='nsew')

    # Control of peak threshold
    field_label_peak_control = ttk.Label(
        frame7,
        text="Peak threshold:",
        anchor="center",
        justify="center",
        style="Box.TLabel"
    )

    field_label_peak_control.grid(
        row=2, column=0, padx=5, pady=2, sticky='nsew')
    field_peaks = ttk.Entry(frame7, justify='center')
    field_peaks.insert(0, peak_value)
    peak_value = float(field_peaks.get())
    field_peaks.grid(row=3, column=0, padx=5, pady=2, sticky='nsew')
    # Control of distance between peaks
    field_label_sep_control = ttk.Label(
        frame7,
        text="Peak separation:",
        anchor="center",
        justify="center",
        style="Box.TLabel"
    )

    field_label_sep_control.grid(
        row=0, column=1, padx=5, pady=2, sticky='nsew')
    field_sep = ttk.Entry(frame7, justify='center')
    field_sep.insert(0, peak_sep)
    peak_sep = float(field_sep.get())
    field_sep.grid(row=1, column=1, padx=5, pady=2, sticky='nsew')

    # Control of peak prominence

    field_label_prominence_control = ttk.Label(
        frame7,
        text="Peak prominence:",
        anchor="center",
        justify="center",
        style="Box.TLabel"
    )

    field_label_prominence_control.grid(
        row=2, column=1, padx=5, pady=2, sticky='nsew')
    field_promi = ttk.Entry(frame7, justify='center')
    field_promi.insert(0, peak_promi)
    peak_promi = float(field_sep.get())
    field_promi.grid(row=3, column=1, padx=5, pady=2, sticky='nsew')

    # max number of peaks:

    field_label_peak_number_control_label = ttk.Label(
        frame6,
        text="Max. number \nof peaks:",
        anchor="center",
        justify="center",
        style="Box.TLabel"
    )

    field_label_peak_number_control_label.grid(
        row=0, column=0, padx=5, pady=2, sticky='nsew')
    field_peak_number_control = ttk.Entry(frame6, justify='center')
    field_peak_number_control.insert(0, peak_number_control)
    peak_number_control = float(field_peak_number_control.get())
    field_peak_number_control.grid(
        row=1, column=0, padx=5, pady=2, sticky='nsew')

    # add button
    button_peak_detection = tk.Button(frame6, text='Peak detection',
                                      command=button_peak_detection_clicked,
                                      state="disabled")
    button_peak_detection.grid(row=2, column=0, padx=10, pady=2, sticky='nsew')

    # Create the separator

    field_label_button_sep = ttk.Label(
        frame6,
        text="==============",
        anchor="center",
        justify="center",
        style="Box.TLabel"
    )
    field_label_button_sep.grid(row=3, column=0, padx=5, pady=2, sticky='nsew')
    # Create button
    button_peak_adding = tk.Button(frame6, text='Add peaks',
                                   command=button_adding_clicked,
                                   state="disabled")
    button_peak_adding.grid(row=4, column=0, padx=10, pady=2, sticky='nsew')

    ###########################################################################
    ##### Manual                                                           ####
    ###########################################################################


    field_label_add_peak = ttk.Label(
        tab4,
        text="Peak position in 1/cm:",
        anchor="center",
        justify="center",
        style="Box.TLabel"
    )

    field_label_add_peak.grid(
        row=0, column=0, padx=5, pady=2, sticky='nsew')
    new_peak_entry = ttk.Entry(tab4, justify='center')    
    new_peak_entry.grid(row=1, column=0, padx=5, pady=2, sticky='nsew')

    # add button to add peak
    # Create button
    button_manual_peak_adding = tk.Button(tab4, text='Add peaks',
                                command=button_manual_adding_clicked,
                                state="disabled")
    button_manual_peak_adding.grid(row=1, column=1, padx=10, pady=2, sticky='nsew')


    ###########################################################################
    ##### Sub panel 4                                                   ####
    ###########################################################################

    # Subpanel 4
    subpanel4.grid_rowconfigure(1, weight=1)  # Configure row weight
    subpanel4.grid_columnconfigure(0, weight=1)  # Configure column weight

    box_frame_4 = ttk.Frame(subpanel4, borderwidth=1, relief="groove")
    box_frame_4.grid(row=1, column=0, padx=10, pady=2, sticky='nsew')
    # Configure column weights of the box_frame
    box_frame_4.grid_columnconfigure(0, weight=1)
    box_frame_4.grid_columnconfigure(1, weight=1)

    # First child frame
    frame9 = ttk.Frame(box_frame_4, borderwidth=2, relief="groove")
    frame9.grid(row=1, column=0, padx=10, pady=2, sticky='nsew')
    frame9.grid_columnconfigure(0, weight=1)
   
    # subpanel4

    field_processing = tk.Label(subpanel4, text='Peak processing:')
    field_processing.config(font=("bold", 14))
    field_processing.config(font=("underline"))
    field_processing.grid(row=0, column=0, padx=10, pady=2, sticky='nsew')

      
    # Buttons subpanel4
    # Create a style
   
    style.configure('TNotebook.Tab', 
                font=('Helvetica', 12, 'bold'))  # Bold font with size 18
    fit_tab = ttk.Notebook(frame9)
    fit_tab.grid(row=0, column=0, padx=10, pady=2, sticky='nsew')
    tab5 = ttk.Frame(fit_tab)
    tab6 = ttk.Frame(fit_tab)
    fit_tab.add(tab5, text="Single peak processing")
    fit_tab.add(tab6, text="Batch proceesing")
    tab5.grid_columnconfigure(0, weight=1)
    tab6.grid_columnconfigure(0, weight=1)
    # Create button in tab 3
 
    button_peak_processing = tk.Button(tab5, text='Peak fitting menu',
                                       command=button_peak_analizer_clicked,
                                       state="disabled")
    button_peak_processing.grid(
        row=0, column=0, padx=10, pady=2, sticky='e')
    # Create button in tab 4
    button_load_batch_folder = tk.Button(tab6, text='Load folder',
                                         command=load_batch_folder,
                                         state="disabled"
                                         )
    button_load_batch_folder.grid(
        row=0, column=0, padx=10, pady=2, sticky='nsew')

    button_batch = tk.Button(tab6, text='Batch processing',
                             command=button_batch_proccesing,
                             state="disabled")
    button_batch.grid(row=0, column=1, padx=10, pady=2, sticky='e')

    button_dashboard = tk.Button(tab6, text='Summary Dashboard',
                             command=dashboard,
                             state="disabled")
    button_dashboard.grid(row=1, column=1, padx=10, pady=2, sticky='e')

    return main_panel
