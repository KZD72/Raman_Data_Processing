# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 12:53:27 2023

This creates a smple GUI for the single peak fitting panel

@author: J.Anaya

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

# Raman_single_peak_fit_GUI

import tkinter as tk
import numpy as np
import time
import re
from tkinter import ttk
from tkinter import messagebox
import os
import scipy.integrate as spi
import csv

from bin import Raman_dataloader
from bin import Raman_fit
from bin import Raman_plot


def save_text(text_to_save):
    """
    Saves the provided text to a file using the Raman_dataloader module.

    Args:
        text_to_save (str or list): The text content to be saved. It can be a single string or a list of strings.

    Returns:
        bool: True if the text was successfully saved, False otherwise.
    """
    if isinstance(text_to_save, str):
        # If the input is a single string, save it directly
        return Raman_dataloader.save_text(text_to_save)
    elif isinstance(text_to_save, list):
        # If the input is a list of strings, concatenate them into a single string
        concatenated_text = "\n".join(text_to_save)
        return Raman_dataloader.save_text(concatenated_text)
    else:
        # Handle invalid input types
        raise TypeError(
            "Invalid input type. 'text_to_save' must be a string or a list.")


def error(text):
    """
    Display an error message based on the given type.

    Args:
        text (str): The message to show.


    """
    messagebox.showerror("Error", text)


def check_range_f(x, new_peak):
    c0 = check_number_f(new_peak)
    c1 = float(new_peak) > x[0] and float(new_peak) < x[-1]
    if c0 and c1:
        return True
    else:
        error("Keep the fine tunning\nwithin the data range")
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
        error("Introduce a valid number")
        return False


def dict_to_list(dictionary):
    """
    Converts the values of a dictionary to a NumPy array.

    Args:
        dictionary (dict): The input dictionary.

    Returns:
        numpy.ndarray: The NumPy array containing the dictionary values.
    """
    return np.array(list(dictionary.values()), dtype=object)


def pad_dictionary(dictionary):
    """
    Fill missing keys in each dictionary entry with 'float("nan")' while maintaining alphabetical order.

    Parameters:
        dictionary (dict): Dictionary to be filled.

    Returns:
        filled_dict (dict): Dictionary with filled entries.

    Example:
        your_dict = {
            'aa': {'a': 2},
            'bb': {'a': 5, 'b': 6, 'c': 7},
            'cc': {'a': 8, 'b': 6, 'd': 7},
            'dd': {'a': 11, 'b': 12, 'c': 13},
            'ee': {'a': 14, 'c': 15}
        }

        filled_dict = fill_missing_keys(your_dict)
    """
    keys = list(dictionary.keys())
    all_keys = sorted(set().union(*[dictionary[key] for key in keys]))

    filled_dict = {}
    for key in keys:
        entry = dictionary[key]
        filled_entry = {**{k: float("nan") for k in all_keys}, **entry}
        filled_dict[key] = filled_entry

    return filled_dict


def calculate_quotients(raw_dictionary, param):
    """
    Calculate the quotients of a specific parameter with each key as the denominator.

    Parameters:
        raw_dictionary (dict): Dictionary containing the data.
        param (str): Parameter to calculate the quotients.

    Returns:
        quotients (ndarray): Matrix of quotients.
        keys (list): List of keys from the dictionary.
        params (list): List of parameters from the dictionary.

    Example:
        your_dict = {
            'aa': {'a': 2, 'b': 3, 'c': 4},
            'bb': {'a': 5, 'b': 6, 'c': 7},
            'cc': {'a': 8, 'c': 10},
            'dd': {'a': 11, 'b': 12, 'c': 13},
            'ee': {'a': 14, 'b': 14, 'c': 15}
        }

        quotients, keys, params = calculate_quotients(your_dict, 'a')
    """
    dictionary = pad_dictionary(raw_dictionary)
    keys = list(dictionary.keys())  # Get the keys from the dictionary

    # Get the parameters from the first key
    params = list(dictionary[keys[0]].keys())
    # Create a 2D array of values from the dictionary
    values = np.array([list(dictionary[key].values())
                      for key in keys], dtype='float64')

    # Find the index of the specified parameter
    param_index = params.index(param)

    # Initialize an array to store the quotients
    quotients = np.zeros((len(keys), len(keys)), dtype=object)

    # Calculate the quotients
    for i in range(len(keys)):
        if param in dictionary[keys[i]]:
            if values[i, param_index] > 1e-12:
                quotients[:, i] = values[:, param_index] / \
                    values[i, param_index]
            else:
                quotients[:, i] = np.nan
        else:
            quotients[:, i] = np.nan

    return quotients, keys, params


def create_table(root, data, header):
    """
    Create a table in a Tkinter window.

    Parameters:
        root (tk.Tk): Root Tkinter window object.
        data (list): 2D list of data to populate the table.
        header (list): List of column headers.

    Returns:
        str: Formatted table as a string.

    Example:
        root = tk.Tk()
        data = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        header = ['A', 'B', 'C']
        table_text = create_table(root, data, header)
        print(table_text)
        root.mainloop()
    """
    # Get the number of rows and columns in the data
    num_rows = len(data)
    num_columns = len(data[0])

    table_text = ''

    # Create labels for column headers
    for j in range(num_columns):
        label = tk.Label(root, text=f"{header[j]}", relief=tk.RIDGE, width=12)
        label.grid(row=0, column=j + 1)

    # Create labels for row headers and data cells
    for i in range(num_rows):
        # Create label for row header
        label = tk.Label(root, text=f"{header[i]}", relief=tk.RIDGE, width=8)
        label.grid(row=i + 1, column=0)

        # Create labels for data cells
        for j in range(num_columns):
            formatted_value = f"{data[i][j]:.3f}"
            label = tk.Label(root, text=formatted_value,
                             relief=tk.SOLID, width=12)
            label.grid(row=i + 1, column=j + 1)

            # Add the formatted value to the table text
            table_text += formatted_value

            # Add a comma separator if not the last column
            if j < num_columns - 1:
                table_text += ', '

        # Add a new line after each row
        table_text += '\n'

    return table_text


def create_plain_table(data, header):
    """
    Create a string representation of a table.

    Parameters:
        data (list): 2D list of data to populate the table.
        header (list): List of column headers.

    Returns:
        str: Formatted table as a string.
    """
    # Get the number of rows and columns in the data
    num_rows = len(data)
    num_columns = len(data[0])

    table_text = ''

    # Add column headers to the table text
    table_text += ', '.join(header) + '\n'

    # Add data cells to the table text
    for i in range(num_rows):
        row_values = [f"{data[i][j]:.3f}" for j in range(num_columns)]
        table_text += ', '.join(row_values) + '\n'

    return table_text


def expanded_dict(dictionary, entries, label):
    """
    Add new subentry to each entry in the dictionary.

    Parameters:
        dictionary (dict): Original dictionary.
        entries (list): List of new entries to be added.
        label (str): Label for the new subentry.

    Returns:
        dict: Updated dictionary with the new subentry.

    Example:
        dictionary = {
            1: {'Gaussian_component': 0.65899564, 'Center': 1406.31829, 'FWHM': 69.6858762, 'Intensity': 7152.90713},
            2: {'Gaussian_component': 0.36174242, 'Center': 437.299737, 'FWHM': 10.925679, 'Intensity': 6093.78661},
            ...
        }
        entries = [100, 200, ...]
        label = 'Integrated Intensity'
        expanded_dict(dictionary, entries, label)
    """
    iterator = 0
    for key, subentry in dictionary.items():
        subentry[label] = entries[iterator]
        iterator = iterator+1
    return dictionary


def voigt_fix_dic(dictionary):
    """
    Calculates and adds the 'FWHM' entry to the dictionary using the Voigt function approximation with ~0.02% accuracy.

    The Voigt function approximation is based on John F. Kielkopf's approach (1973),
    "New approximation to the Voigt function with applications to spectral-line profile analysis",
    Journal of the Optical Society of America, 63 (8): 987, doi:10.1364/JOSA.63.000987

    Parameters:
        dictionary (dict): The input dictionary containing the subentries.

    Returns:
        dict: The dictionary with the 'FWHM' entry added to the subentries where 'Gauss_FWHM' and 'Lorentz_FWHM' exist.

    Example:
        input_dict = {
            1: {'Gauss_FWHM': 368.614301, 'Lorentz_FWHM': 116.375262, 'Center': 1406.16995},
            2: {'Gauss_FWHM': 10.9247146, 'Lorentz_FWHM': 20.346789, 'Center': 437.299772},
            3: {'Gauss_FWHM': 25.7036236, 'Lorentz_FWHM': 50.987631, 'Center': 1251.33855}
        }
        output_dict = voigt_fix_dic(input_dict)
        print(output_dict)
    """

    for subentry in dictionary.values():
        if 'Gauss_FWHM' in subentry and 'Lorentz_FWHM' in subentry:
           
            # Calculate the FWHM using the Voigt function approximation
            fwhm = 0.5346 * subentry['Lorentz_FWHM'] + np.sqrt(
                subentry['Gauss_FWHM'] * subentry['Gauss_FWHM'] + 0.2166 * subentry['Lorentz_FWHM'] * subentry['Lorentz_FWHM'])
            # Insert the 'FWHM' subentry after 'Center'
            subentry_keys = list(subentry.keys())
            center_index = subentry_keys.index('Center')
            subentry_keys.insert(center_index + 1, 'FWHM')
            subentry_values = list(subentry.values())
            subentry_values.insert(center_index + 1, fwhm)
            subentry.clear()
            subentry.update(zip(subentry_keys, subentry_values))
        if 'Asymmetry' in subentry:
            #Correction of the FWHM from both functions in bimodal
            subentry['FWHM']=subentry['FWHM']/2+subentry['Asymmetry']*subentry['FWHM']/2

    return dictionary


def on_popup_close(popup):
    popup.grab_release()  # Release the grab
    popup.destroy()


def extract_peak_info(text):
    """
    Extracts peak information from the given text and organizes it into a dictionary.

    Args:
        text (str): The text containing peak information.

    Returns:
        dict: A dictionary containing peak information, organized by peak number and parameter.
    """
    peak_info = {}
    pattern = r'Peak_(\d+)_([^\s:]+):\s+([\d.+-]+)'

    # Find all matches of the pattern in the text
    matches = re.findall(pattern, text)

    for match in matches:
        peak_number = int(match[0])
        parameter = match[1]
        value = float(match[2])

        # Check if the peak number is already present in the dictionary
        if peak_number not in peak_info:
            peak_info[peak_number] = {}

        # Add the parameter and its corresponding value to the peak_info dictionary
        peak_info[peak_number][parameter] = value

    return peak_info


def batch_fit(canvas, canvas_panel, info, x, y, peaks, file_path,models=[], silent=False):
    # Extract peak positions and models from UI elements
    np_peaks = np.asarray(peaks, dtype='float64')
    peak_positions = [np_peaks[i, 0] for i in range(len(np_peaks))]

    # Only Gauss-Lorentz for now
    if len(models)!=len(peak_positions):
        models = ['Gauss-Lorentz' for item in peak_positions]
  
    # Find the index corresponding to each peak position
    index = [np.searchsorted(np.asarray(x, dtype='float64'), peak)
             for peak in peak_positions]

    # Get the y values at the index positions
    y_pos = [y[index_val] for index_val in index]

    # Create peak_info list with peak positions and corresponding y values
    peak_info = [[peak_positions[item], y_pos[item]]
                 for item in range(len(peak_positions))]

    # Launch the fitting process
    # Clear previous fit information

    fitting = Raman_fit.fit(x, y, peak_info, models)

    # Get fit info
    info_fit = Raman_fit.fit_info(fitting)

    # Get fit data
    y_fit = Raman_fit.model_f(fitting.params, x, peak_info, models)

    # Plots
    plots_to_show = []
    plots_to_show.append([x, y])
    plots_to_show.append([x, y_fit])

    leyend = ['Baseline_corrected', 'Model fit']
    # individual
    # Format and display fit information
    formated_info = extract_peak_info(info_fit)

    peak_final_label = [
        f'{params["Center"]:.1f}'for params in formated_info.values()]

    int_val = []
    inner_iter = 0
    for item in range(len(models)):
        new_model = ['Not used' for item in models]
        new_model[item] = models[item]

        if new_model[item] != 'Not used':
            plots_to_show.append(np.array([x, Raman_fit.model_f(
                fitting.params, x, peak_info, new_model)], dtype='object'))
            leyend.append("P"+str(item+1)+"_"+peak_final_label[inner_iter])
            inner_iter = inner_iter+1
            # Get the integral:
            # Define the function to be integrated

            def f(x):
                return Raman_fit.model_f(fitting.params, x, peak_info, new_model)
            result, error = spi.quad(f, x[0], x[-1])
            int_val.append(result)
            # Extract the FWHM of the voight profile:

    # add the integrated intesity to dictionary:

    formated_info = expanded_dict(
        formated_info, int_val, 'Integrated Intensity')
    # Add the FWHM if I have a voigt profile:
    formated_info = voigt_fix_dic(formated_info)
    # print(formated_info)
    # Update plot

    result_dict = {info['Title']: {'title':info['Title'],
                                   'leyends': leyend,
                                   'plots_to_show': plots_to_show, 
                                   'fit_results':formated_info, 
                                   'post_process_results':formated_info,
                                   'full_fit_info':info_fit}}
    
   
    if silent:        
        print("Updated graph")
        fig, ax = Raman_plot.plotter(plots_to_show,
                                ["Wavenumber (1/cm)",
                                "Intensity (A.U.)"],
                                info['Title'],
                                leyends=leyend,
                                lines=True,
                                res=150,
                                leyend_frame=[True, 'b']
                                )
        Raman_plot.update_plot(canvas, canvas_panel, fig, ax, plots_to_show)
   
    return result_dict


###############################################################################
def create_fit_panel(main_window, canvas, canvas_panel, info, x, y, peaks):
    global peak_list_entries, peak_labels, model_list, info_fit, time_stamp, formated_info, dialog
    np_peaks = np.asarray(peaks, dtype='float64')
    x_peak = np_peaks[:, 0]
    y_peak = np_peaks[:, 1]
    peak_labels = ["{:.2f}".format(peak) for peak in x_peak]
    peak_list_entries = []
    model_list = []
    info_fit = "No info"
    time_stamp = "0.0s"
    formated_info = "None"
    dialog = None

    def on_close():
        """
        Handles the closing event of the main window.
        Performs necessary cleanup actions before closing the application.
        """
        global model_list
        fit_window.quit()  # Quit the main window event loop
        fit_window.destroy()  # Destroy the main window
        


        
    def field_creator(frame, position, peak, peak_label):
        """
       Creates and configures a set of UI elements for displaying and editing a peak.

       Args:
           frame (ttk.Frame): The parent frame where the UI elements will be placed.
           position (int): The position of the peak in the list.
           peak (str): The initial peak value to be displayed.
           peak_label (str): The label associated with the peak.

       Returns:
           None
       """
        global peak_list_entries, peak_labels, model_list

        label = ttk.Label(
            frame,
            text=f"Peak_{position} detected at:",
            anchor="center",
            justify="center",
            style="Box.TLabel"
        )

        label.grid(row=position+1, column=0, padx=5, pady=5, sticky='nsew')

        label_entry = ttk.Entry(frame, justify='center')
        label_entry.insert(0, peak)  # Initial peak entry
        label_entry.grid(row=position+1, column=1,
                         padx=5, pady=5, sticky='nsew')

        # Store the Entry widget in the peak_list_entries dictionary using the peak_label as the key
        peak_list_entries.append(label_entry)

        label_units = ttk.Label(
            frame,
            text="1/cm",
            anchor="center",
            justify="center",
            style="Box.TLabel"
        )
        label_units.grid(row=position+1, column=2,
                         padx=5, pady=5, sticky='nsew')

        # Add the combo box to select the fitting method
        label_model = ttk.Label(frame, text="Select the peak physics model:")
        label_model .grid(row=0, column=3, padx=5, pady=5)

        options = ['Not used', 'Gaussian', 'Lorentz',
                   'Gauss-Lorentz','Asy-Gauss-Lorentz', 'Voigt', 'Fano-Simply', 'Fano-Voigt']

        combobox = ttk.Combobox(frame, values=options)
        combobox.set(options[3])
        model_list.append(combobox)
        combobox.grid(row=position+1, column=3, padx=5, pady=5)

    def button_fit_clicked():
        """
         Handles the event when the fit button is clicked.
         Performs peak fitting on the data and updates the plot and fit information.

         Returns:
             None
             """
        global peak_list_entries, x_peak, y_peak, info_fit, time_stamp, formated_info

        # Extract peak positions and models from UI elements
        peak_positions = [float(peak_list_entrie.get(
        )) for peak_list_entrie in peak_list_entries if check_range_f(x, peak_list_entrie.get())]
        if len(peak_positions) == len(peak_list_entries):
            models = [model.get() for model in model_list]

            # Find the index corresponding to each peak position
            index = [np.searchsorted(np.asarray(x, dtype='float64'), peak)
                     for peak in peak_positions]

            # Get the y values at the index positions
            y_pos = [y[index_val] for index_val in index]

            # Create peak_info list with peak positions and corresponding y values
            peak_info = [[peak_positions[item], y_pos[item]]
                         for item in range(len(peak_positions))]

            # Launch the fitting process
            # Clear previous fit information
            text_widget_2.delete('1.0', 'end')
            text_widget_2.insert("end", "Computing...", "bold")
            main_window.update_idletasks()  # Update the display
            start = time.perf_counter()
            fitting = Raman_fit.fit(x, y, peak_info, models)
            end = time.perf_counter()
            time_stamp = "\nThe execution time of fit_method is {} seconds.\n".format(
                end - start)
            text_widget_2.insert("end", f"{time_stamp}\n", "bold")

            # Get fit info
            info_fit = Raman_fit.fit_info(fitting)

            # Format and display fit information
            formated_info = info_display()
            # print(formated_info)
            # Get fit data
            y_fit = Raman_fit.model_f(fitting.params, x, peak_info, models)
            # r2:
            r_2 = 1-(fitting.residual**2).sum() / \
                (sum(np.power(y-np.mean(y), 2)))
            text_widget_2.insert(
                "end", f"Fit r²={r_2}, check fit details for more info.\n\n", "bold")
            # Plots
            plots_to_show = []
            plots_to_show.append([x, y])
            plots_to_show.append([x, y_fit])
            leyend = ['Baseline_corrected', 'Model fit']
            # individual

            peak_final_label = [
                f'{params["Center"]:.1f}'for params in formated_info.values()]

            int_val = []
            inner_iter = 0
            for item in range(len(models)):
                new_model = ['Not used' for item in models]
                new_model[item] = models[item]

                if new_model[item] != 'Not used':
                    plots_to_show.append(np.array([x, Raman_fit.model_f(
                        fitting.params, x, peak_info, new_model)], dtype='object'))
                    leyend.append("P"+str(item+1)+"_" +
                                  peak_final_label[inner_iter])
                    inner_iter = inner_iter+1
                    # Get the integral:
                    # Define the function to be integrated

                    # def f(x):
                    #     return Raman_fit.model_f(fitting.params, x, peak_info, new_model)

                    # result, error = spi.quad(f, x[0], x[-1])
                    result=np.trapz(Raman_fit.model_f(fitting.params, x, peak_info, new_model))
                    int_val.append(result)
                    # Extract the FWHM of the voight profile:

            # add the integrated intesity to dictionary:

            formated_info = expanded_dict(
                formated_info, int_val, 'Integrated Intensity')
            # Add the FWHM if I have a voigt profile:
            formated_info = voigt_fix_dic(formated_info)
            # print(formated_info)
            # Update plot
            fig, ax = Raman_plot.plotter(plots_to_show,
                                         ["Wavenumber (1/cm)",
                                          "Intensity (A.U.)"],
                                         info['Title'],
                                         leyends=leyend,
                                         lines=True,
                                         res=150,
                                         leyend_frame=[True, 'b']
                                         )
            Raman_plot.update_plot(canvas, canvas_panel,
                                   fig, ax, plots_to_show)

            for peak_number, parameters in formated_info.items():
                text_widget_2.insert(
                    "end", f"Peak {peak_number}:" + "\n", "bold")
                for parameter, value in parameters.items():
                    text_widget_2.insert(
                        "end", f"{parameter}: {value}" + "\n", "bold")
                text_widget_2.insert("end", "\n", "bold")

            # Activate fit info button
            button_fit_info.config(state="normal")
            button_fit_info_save.config(state="normal")
            button_postpro.config(state="normal")
            return model_list
        else:
            error("Introduce a valid number")

    def fit_info_clicked():
        """
        Handles the event when the fit info button is clicked.
        Displays fit information in a popup window and provides an option to save the information.

        Returns:
            None
        """
        global info_fit, time_stamp

        # Create a new popup window
        popup_window = tk.Toplevel(main_window)
        popup_window.title("Fit data")
        popup_window.grab_set()  # Grab the focus for the popup window
        popup_window.protocol("WM_DELETE_WINDOW",
                              lambda: on_popup_close(popup_window))

        # Create a frame to hold the text widget and save button
        box_frame_pw = tk.Frame(popup_window, borderwidth=1, relief="groove")
        box_frame_pw.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

        # Create a text widget to display the fit information
        text_widget = tk.Text(box_frame_pw, padx=20, pady=10)
        text_widget.tag_configure("bold", font=("TkDefaultFont", 11, "bold"))
        text_widget.tag_configure("red", foreground="blue")

        # Display fitting time
        text_widget.insert("end", "Fitting time " + time_stamp + "\n", "bold")

        # Display fit information
        text_widget.insert("end", info_fit, "bold")

        text_widget.grid(row=0, column=0, sticky='nsew')

        # Create a new style for the button
        save_button_style = ttk.Style()
        save_button_style.configure(
            'Save.TButton', foreground='Blue', font=('Arial', 14))

        # Create the "Save" button
        save_button = ttk.Button(popup_window, text="Save", command=lambda: save_text(
            info_fit), style='Save.TButton')
        save_button.grid(row=1, column=0, sticky='nsew')

        popup_window.mainloop()

    def info_display():
        """
        Extracts and filters peak information from the global variable info_fit.

        Returns:
            dict: A dictionary containing filtered peak information, organized by peak number and parameter.
        """
        global info_fit, time_stamp

        # Extract peak information from info_fit using extract_peak_info function
        filtered_text = extract_peak_info(info_fit)

        return filtered_text

    def button_postpro_clicked():
        global formated_info

        if len(formated_info.items()) > 1:
            # get fitting parameters:
            # Get the keys from the dictionary
            keys = list(formated_info.keys())
            # Get the parameters from the first key
            param_names = list(formated_info[keys[0]].keys())

            # Create a window to show the tables:
            # Create a new popup window
            popup_window_pp = tk.Toplevel(main_window)
            popup_window_pp.title("Fit data")
            popup_window_pp.grab_set()  # Grab the focus for the popup window
            popup_window_pp.protocol(
                "WM_DELETE_WINDOW", lambda: on_popup_close(popup_window_pp))

            # Create a notebook widget to hold the tabs
            notebook = ttk.Notebook(popup_window_pp)
            notebook.pack(fill=tk.BOTH, expand=True)
            table_to_save = []
            iterator = 0
            for parameter in param_names:
                tab1 = ttk.Frame(notebook)
                notebook.add(tab1, text=parameter)

                # Create a frame to hold the text widget and save button
                box_frame_pw = tk.Frame(tab1, borderwidth=1, relief="groove")
                box_frame_pw.grid(row=0, column=0, padx=10,
                                  pady=10, sticky='nsew')
                table_to_save.append(parameter+":\n")
                labels = ["P"+str(item+1)
                          for item in range(len(formated_info.items()))]

                table_to_save.append(create_table(box_frame_pw, calculate_quotients(
                    formated_info, parameter)[iterator], labels))
                table_to_save.append("\n\n")

            box_frame_pwb = tk.Frame(
                popup_window_pp, borderwidth=1, relief="groove")
            box_frame_pwb.pack(fill=tk.BOTH, expand=True)
            # Create a new style for the button
            save_button_style = ttk.Style()
            save_button_style.configure(
                'Save.TButton', foreground='Blue', font=('Arial', 14))

            # Create the "Save" button
            save_button = ttk.Button(box_frame_pwb, text="Save",
                                     command=lambda: save_text(table_to_save), style='Save.TButton')
            save_button.pack(fill=tk.BOTH, expand=True)

            popup_window_pp.mainloop()
        else:
            error("Not enough peaks to compare between them")

    def save_fit_info_clicked():
        # Retrieve the text from the widget
        text_to_save = text_widget_2.get("1.0", "end-1c")
        save_text(text_to_save)

    #####################################################################
    ####                Create window for fitting:                    ###
    #####################################################################
    # Area to create the peak fittingwindow
    fit_window = tk.Toplevel(main_window)
    fit_window.title('Raman peak analizer')
    fit_window.geometry("755x900")
    fit_window.resizable(False, False)  # Disable resizing
    fit_window.attributes("-topmost", True)
    fit_window.protocol("WM_DELETE_WINDOW", on_close)

    # Grid layout configuration
    fit_window.grid_columnconfigure(0, weight=1)
    fit_window.grid_rowconfigure(0, weight=1)

   # Creation of main panel elements
    main_panel = tk.Frame(fit_window, bg='white')
    main_panel.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

    # Grid layout configuration
    main_panel.grid_columnconfigure(0, weight=1)
    main_panel.grid_rowconfigure(0, weight=1)

    # Create subpanel frames within the main panel
    subpanel1 = tk.Frame(main_panel, bg='lightblue')
    subpanel2 = tk.Frame(main_panel, bg='lightgrey')

    # Grid layout configuration
    main_panel.grid_columnconfigure(0, weight=1)
    main_panel.grid_rowconfigure(0, weight=1)
    main_panel.grid_rowconfigure(1, weight=1)

    # Grid placement of subpanel frames
    subpanel1.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')
    # Subpanel 1
    subpanel1.grid_rowconfigure(1, weight=1)  # Configure row weight
    subpanel1.grid_columnconfigure(0, weight=1)  # Configure column weight

    subpanel2.grid(row=1, column=0, padx=10, pady=10, sticky='nsew')
    # Subpanel 2
    subpanel2.grid_rowconfigure(1, weight=1)  # Configure row weight
    subpanel2.grid_columnconfigure(0, weight=1)  # Configure column weight

    # Create the decoration of the subpanels
    # Create a style for the labels
    style = ttk.Style()
    style.configure("Box.TLabel", background=main_window["background"])

    box_frame = ttk.Frame(subpanel1, borderwidth=1, relief="groove")
    box_frame.grid(row=1, column=0, padx=10, pady=10, sticky='nsew')
    # Configure box_frame to expand in both directions
    box_frame.grid_rowconfigure(0, weight=1)
    box_frame.grid_columnconfigure(0, weight=1)

    # Create a new frame as a container for frame1 and the scrollbar
    container_frame = ttk.Frame(box_frame, borderwidth=1, relief="groove")
    container_frame.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

    # Configure container_frame to expand in both directions
    container_frame.grid_rowconfigure(0, weight=1)
    container_frame.grid_columnconfigure(0, weight=1)

    # Create a canvas widget inside the container frame
    canvas_peaks = tk.Canvas(container_frame)
    canvas_peaks.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')
    # Configure canvas_peaks to expand in both directions
    canvas_peaks.grid_rowconfigure(0, weight=1)
    canvas_peaks.grid_columnconfigure(0, weight=1)

    # Create a scrollbar widget
    scrollbar = ttk.Scrollbar(
        container_frame, orient="vertical", command=canvas_peaks.yview)
    scrollbar.grid(row=0, column=1, sticky='ns')
    # Add this line to make the scrollbar appear shaded
    scrollbar.config(takefocus=0)

    # Configure the canvas to use the scrollbar
    canvas_peaks.configure(yscrollcommand=scrollbar.set)

    # Create a frame to hold the fields
    frame1 = ttk.Frame(canvas_peaks, borderwidth=2, relief="groove")
    # Configure canvas_peaks to expand in both directions
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=1)

    # Add the frame to the canvas using grid
    canvas_peaks.grid_rowconfigure(0, weight=1)
    canvas_peaks.grid_columnconfigure(0, weight=1)
    canvas_peaks.create_window(
        (0, 0), window=frame1, anchor="nw", tags="frame1")

    # Configure the canvas to adjust scroll region when the frame size changes
    frame1.bind("<Configure>", lambda event: canvas_peaks.configure(
        scrollregion=canvas_peaks.bbox("all")))

    # Add elements to frame1
    [field_creator(frame1, item+1, x_peak[item], peak_labels[item])
     for item in range(len(x_peak))]

    # Second child frame
    frame2 = ttk.Frame(box_frame, borderwidth=2, relief="groove")
    frame2.grid(row=0, column=1, padx=10, pady=10, sticky='nsew')

    button_fit = tk.Button(frame2, text='Fit data', command=button_fit_clicked)
    button_fit.grid(row=0, column=0, padx=10, pady=10)

    # Postprocesing all data:
    button_postpro = tk.Button(frame2, text='Postpro\nPeak data',
                               command=button_postpro_clicked,
                               state="disabled")
    button_postpro.grid(row=1, column=0, padx=10, pady=10)

    # Second panel configuration
    box_frame_2 = ttk.Frame(subpanel2, borderwidth=1, relief="groove")
    box_frame_2.grid(row=1, column=0, padx=10, pady=10, sticky='nsew')
    # Grid layout configuration
    box_frame_2.grid_columnconfigure(0, weight=1)
    box_frame_2.grid_rowconfigure(0, weight=1)
    # Create a new frame as a container for frame3 and the scrollbar
    container_frame_2 = ttk.Frame(box_frame_2, borderwidth=1, relief="groove")
    container_frame_2.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')
    # Grid layout configuration
    container_frame_2.grid_columnconfigure(0, weight=1)
    container_frame_2.grid_rowconfigure(0, weight=1)
    # Create a canvas widget inside the container frame
    canvas_peak_info = tk.Canvas(container_frame_2)
    canvas_peak_info.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')
    # Grid layout configuration
    canvas_peak_info.grid_columnconfigure(0, weight=1)
    canvas_peak_info.grid_rowconfigure(0, weight=1)
    # Create a scrollbar widget
    scrollbar_2 = ttk.Scrollbar(
        container_frame_2, orient="vertical", command=canvas_peak_info.yview)
    scrollbar_2.grid(row=0, column=1, sticky='ns')
    # Add this line to make the scrollbar appear shaded
    scrollbar_2.config(takefocus=0)

    # Configure the canvas to use the scrollbar
    canvas_peak_info.configure(yscrollcommand=scrollbar_2.set)

    # Create a frame to hold the fields
    frame3 = ttk.Frame(canvas_peak_info, borderwidth=2, relief="groove")
    # Grid layout configuration
    frame3.grid_columnconfigure(0, weight=1)
    frame3.grid_rowconfigure(0, weight=1)

    # Create a frame to hold the fields
    canvas_peak_info.grid_rowconfigure(0, weight=1)
    canvas_peak_info.grid_columnconfigure(0, weight=1)
    canvas_peak_info.create_window(
        (0, 0), window=frame3, anchor="nw", tags="frame1")

    # Configure the canvas to adjust scroll region when the frame size changes
    frame3.bind("<Configure>",
                lambda event: canvas_peak_info.configure(scrollregion=canvas_peak_info.bbox("all")))

    # First child frame
    text_widget_2 = tk.Text(frame3, padx=20, pady=10)
    text_widget_2.tag_configure("bold", font=("TkDefaultFont", 11, "bold"))
    text_widget_2.tag_configure("red", foreground="blue")
    text_widget_2.insert("end", "Fitting Info:\n", "bold")
    text_widget_2.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

    # Second child frame

    frame4 = ttk.Frame(box_frame_2, borderwidth=2, relief="groove")
    frame4.grid(row=0, column=1, padx=10, pady=10, sticky='nsew')

    button_fit_info_save = tk.Button(frame4, text='Save info',
                                     command=save_fit_info_clicked,
                                     state="disabled")
    button_fit_info_save.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

    button_fit_info = tk.Button(frame4, text='Fit details',
                                command=fit_info_clicked,
                                state="disabled")
    button_fit_info.grid(row=1, column=0, padx=10, pady=10)
