# -*- coding: utf-8 -*-
"""
Created on Wed May 24 16:31:53 2023

Module to create simple plots

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

# Raman_plot.py

import matplotlib.pyplot as plt
import numpy as np
import os
import csv
import re
import sys
# Import FigureCanvasTkAgg class and NavigationToolbar2Tk from matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk  # Import the Tkinter module for GUI
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox


from bin import Raman_dataloader

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

def replace_invalid_characters(input_string):
    # Define a regular expression pattern to match invalid characters
    invalid_char_pattern = re.compile(
        r'[^\x00-\x7F]+')  # Match non-ASCII characters

    # Replace invalid characters with a replacement character (e.g., '?')
    clean_string = invalid_char_pattern.sub('?', input_string)

    return clean_string

# Plotting


def unit_conv(x):
    """
    Parameters
    ----------
    x : float/double
        Number in cm

    Returns
    -------
    double/float
        DESCRIPTIONthe number in inches
    """
    return x/2.54
def error(text):
    """
    Display an error message based on the given type.

    Args:
        text (str): The message to show.


    """
    messagebox.showerror("Error", text)


def adjust_label_positions(labels, positions, threshold=5, fontsize=6):
    adjusted_positions = positions.copy()

    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            if np.abs(adjusted_positions[i] - adjusted_positions[j]) < threshold:
                offset = threshold - \
                    np.abs(adjusted_positions[i] - adjusted_positions[j])
                sign = np.sign(adjusted_positions[j] - adjusted_positions[i])
                adjusted_positions[i] -= sign * offset / 2
                adjusted_positions[j] += sign * offset / 2

    for i in range(len(labels)):
        if i > 0 and np.abs(adjusted_positions[i] - adjusted_positions[i-1]) < threshold:
            adjusted_positions[i] = max(
                adjusted_positions[i-1] + threshold, adjusted_positions[i])

    return adjusted_positions

def create_dropdown(master, dictionary, row, column):
    # Get the labels from the dictionary
    labels = list(dictionary.keys())

    # Create a Tkinter variable
    var = tk.StringVar(master)

    # Set the default option
    var.set(labels[0])

    # Create the dropdown menu using a Combobox
    dropdown = ttk.Combobox(master, textvariable=var, values=labels)
    dropdown.grid(row=row, column=column, padx=10, pady=10, sticky='nsew')  # Place the dropdown at the specified row and column
    return var  # Return the variable so you can get the selected option later

def plotter(data_list, labels, title, error=None, lines=True, leyends=None, size=None, res=None, arrow=None, text=None, leyend_frame=[False, False]):

    # Edit the style
    
    style_file =  resource_path("Resources\plotstyle.mplstyle")  # APS style
    plt.style.use(style_file)
    line_type = ['-', '--', ':']  # intial dashlines
    
    marker_type = ['o', '^', 'v', '<', '>','s', 'd', 'p', 'h', '*']

    data = np.array(data_list, dtype=object)
    # Get the default color cycle from matplotlib
    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Generate a list of 20 differentiated colors
    num_colors = len(data)
    differentiated_colors = default_colors[:num_colors]

    title = replace_invalid_characters(title)

    if size == "double_width":
        plt.rcParams['figure.figsize'] = [
            unit_conv(2*8.6), unit_conv(5)]  # APS style  two colums
    elif size == "double_height":
        plt.rcParams['figure.figsize'] = [
            unit_conv(8.6), unit_conv(2*5)]  # APS style  two rows
    elif size == "double_size":
        plt.rcParams['figure.figsize'] = [
            unit_conv(2*8.6), unit_conv(2*5)]  # APS style  two colums two rows
    elif size == "double_size_double_heigh":
        plt.rcParams['figure.figsize'] = [
            unit_conv(2*8.6), unit_conv(4*5)]  # APS style  two colums four rows
    elif size == "normal":
        plt.rcParams['figure.figsize'] = [
            unit_conv(8.6), unit_conv(5)]  # APS single colums
    else:
        plt.rcParams['figure.figsize'] = [
            unit_conv(16), unit_conv(10)]  # APS single colums

    if res != None:
        plt.rcParams['figure.dpi'] = res  # APS style  

    fig, ax = plt.subplots()

    if leyends == None:
        leyends = ["_" for element in data]

    for iterator in range(len(data)):
        # Plot data with error bars
        if iterator < 2:
            line_iter = iterator
            line_thick = 2.0
        else:
            line_iter = 2
            line_thick = 1.0
        if iterator < len(differentiated_colors):
            colour_iter = iterator
        else:
            colour_iter = len(differentiated_colors)-1
        if lines:
            ax.errorbar(data[iterator][0], data[iterator][1], fmt=line_type[line_iter],
                        color=differentiated_colors[colour_iter], label=leyends[iterator],
                        linewidth=line_thick)  # Increase line thickness to 2.0)
        else:
            ax.errorbar(data[iterator][0], data[iterator][1], fmt=marker_type[line_iter],
                        color=differentiated_colors[colour_iter], label=leyends[iterator],
                        linewidth=line_thick)  # Increase line thickness to 2.0)

    # Set plot properties
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])

    ax.set_title(title)
    ax.legend()

    # Add borders to the plot
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)

    # Add major grid lines
    ax.grid(color='gray', linestyle='--', linewidth=0.5)

    # Add minor grid lines
    ax.minorticks_on()
    ax.grid(visible=True, which='both', color='gray',
            linestyle=':', linewidth=0.2)

    # Adjust the tick parameters
    ax.tick_params(axis='both', which='both',
                   direction='in', length=4, width=1, pad=4)
    ax.tick_params(axis='both', which='minor', length=2, width=1)

    if text != None:
        for item in text:
            plt.text(item)
    if leyend_frame[0]:
        legend_positions = {
            'b': 'best',
            'ur': 'upper right',
            'ul': 'upper left',
            'll': 'lower left',
            'lr': 'lower right',
            'r': 'right',
            'cl': 'center left',
            'cr': 'center right',
            'lc': 'lower center',
            'uc': 'upper center',
            'c': 'center'
        }
        location = legend_positions.get(leyend_frame[1], 'best')
        legend = plt.legend(loc=location, facecolor='lightgray',
                            framealpha=0.5, frameon=True)
        legend.get_frame().set_facecolor('lightgray')

    if arrow is not None:
        ax.set_ylim([1.2*np.min([np.min(iter[1]) for iter in data]),
                    1.5*np.max([np.max(iter[1]) for iter in data])])
        peak_positions = [data[0][0][arrow_data[0]] for arrow_data in arrow]
        peak_values = [arrow_data[1] for arrow_data in arrow]

        adjusted_positions = adjust_label_positions(
            peak_values, peak_positions, 5, 6)

        for label, position in zip(peak_values, adjusted_positions):
            ax.annotate(f'{position:.2f} (1/cm)-{label:.2f} ',
                        xy=(position, label), xytext=(0, 10),
                        textcoords='offset points',
                        arrowprops=dict(arrowstyle='->'),
                        fontsize=6,
                        rotation=90)

    return fig, ax


def fig_clean():
    figures = plt.get_fignums()

    for fig_num in figures:
        fig = plt.figure(fig_num)
        plt.close(fig)

def set_figure_properties(fig, size, res=None):
    """
    Sets the figure size and resolution based on the given parameters.

    Parameters:
    fig (matplotlib.figure.Figure): The figure to modify.
    size (str): The size of the figure. Can be "double_width", "double_height", "double_size", "double_size_double_heigh", or "normal".
    res (int, optional): The resolution (DPI) of the figure. If not provided, the resolution is not changed.
    """
    # APS sizes 
    size_dict = {
        "normal": [unit_conv(8.6), unit_conv(5)],
        "double_width": [unit_conv(2*8.6), unit_conv(5)],
        "double_height": [unit_conv(8.6), unit_conv(2*5)],
        "double_size": [unit_conv(2*8.6), unit_conv(2*5)],
        "double_size_double_heigh": [unit_conv(2*8.6), unit_conv(4*5)]
    }

    # Set the figure size
    fig.set_size_inches(size_dict[size])

    # Set the figure resolution
    if res is not None:
        fig.set_dpi(res)

def update_plot(canvas, canvas_panel, fig, ax, data):
    global dat
    fig_clean()
    dat=data
    ax = fig.get_axes()[0]

    canvas.get_tk_widget().destroy()
    canvas = FigureCanvasTkAgg(fig, master=canvas_panel)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
    # Force an update of the window before updating the canvas
    canvas.get_tk_widget().update_idletasks()

    toolbar = NavigationToolbar2Tk(canvas, canvas_panel, pack_toolbar=False)
    toolbar.update()
    toolbar.grid(row=1, column=0, sticky="ew")

    def log_axes():
        
        ax.set_yscale("log")
        canvas.draw()
       

    def reset_axes():
        global dat
        ax.set_yscale("linear")
        xmin, xmax = ax.get_xlim()  # Get the current limits
        if xmin > xmax:  # If the x-axis is inverted
            ax.invert_xaxis()  # Invert it back
            dat=data            
            canvas.draw()
    
    def reverse_axes():
        global dat
        xmin, xmax = ax.get_xlim()  # Get the current limits
        if xmin < xmax:  # If the x-axis is inverted
            ax.invert_xaxis()  # This will reverse the x-axis
            canvas.draw()
            dat=Raman_dataloader.reorder_data_descending(data,multiple=True,transpose=True)
            
           
                
    def save_text():
        global dat
        file_path = tk.filedialog.asksaveasfile(defaultextension=".csv",
                                                filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])

        rows = []
        iter = 0
        for entry in dat:

            if iter == 0:
                rows.append(entry[0])
                rows.append(entry[1])
                iter = iter+1
            else:
                rows.append(entry[1])
                iter = iter+1

        headers = []
        if len(dat) < 2:
            headers.append(ax.get_xlabel())
            headers.append(ax.get_ylabel())
        elif len(dat) < 3:
            headers.append('Wavenumber (1/cm)')
            headers.append('Raman intensity (1/cm)')
            headers.append('Model Raman intensity (1/cm)')
        else:
            headers.append('Wavenumber (1/cm)')
            headers.append('Raman intensity (1/cm)')
            headers.append('Model Raman intensity (1/cm)')
            [headers.append(f"P_{iterator+1} intesity (1/cm)")
             for iterator in range(len(dat)-2)]
        if file_path is not None:
            with open(file_path.name, "w", newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(headers)  # Write the headers as the first row
                writer.writerows(np.transpose(rows))  # Write the rows of data

    def save_high_res_image():
        local_fig=fig
        # Create a new Tkinter window
        window = tk.Toplevel(canvas_panel.winfo_toplevel())

        #Label
        label = tk.Label(window, text='Select graph options:')
        label.config(font=("bold", 14))
        label.config(font=("underline"))
        label.grid(row=0, column=0, padx=10, pady=2, sticky='nsew')

        # Define the options for size and resolution
        size_options = ["normal", "double_width", "double_height", "double_size", "double_size_double_heigh"]
        res_options = [150, 300, 600, 1200]
        size_dict = {option: option for option in size_options}
        res_dict = {option: option for option in res_options}


        # Create the OptionMenu widgets
        
        size_var=create_dropdown(window, size_dict, 1, 0)
       

        res_var = create_dropdown(window, res_dict, 2, 0)

         # Create three entry fields in the panel
        title_entry = ttk.Entry(window)
        title_entry.insert(0, "Enter title (LaTex allowed)")
        title_entry.grid(row=3, column=0, sticky="nsew")

        x_label_entry = ttk.Entry(window)
        x_label_entry.insert(0, "Enter x-axis label (LaTex allowed)")
        x_label_entry.grid(row=4, column=0, sticky="nsew")

        y_label_entry = ttk.Entry(window)
        y_label_entry.insert(0, "Enter y-axis label (LaTex allowed)")
        y_label_entry.grid(row=5, column=0, sticky="nsew")
        

        # Create the confirmation button
        def confirm():
            # Get the user's selections
            size = size_var.get()
            res = int(res_var.get())

            #Set the title and labels:
            # Set plot properties
            try:
                ax.set_xlabel(r'$' + x_label_entry.get() + '$')
                ax.set_ylabel(r'$' + y_label_entry.get() + '$')

                ax.set_title(r'$' +title_entry.get()+ '$')

                # Set the figure properties
                set_figure_properties(local_fig, size, res)

                # Open a file dialog to save the figure
                file_path = filedialog.asksaveasfilename(defaultextension=".png",
                                                        filetypes=[("PNG files", "*.png"), 
                                                                    ("TIFF files", "*.tiff"), 
                                                                    ("EPS files", "*.eps"),
                                                                    ("PDF files", "*.pdf")])
                if file_path:
                    local_fig.savefig(file_path, dpi=res)
            except:
                error("Error in LaTex format, not saved")

            # Close the window
            window.destroy()

        confirm_button = tk.Button(window, text="Confirm", command=confirm)
        confirm_button.grid(row=6, column=0, padx=10, pady=10, sticky='nsew')
    

   # Define the button styles
           
    button_style = {
        'background': '#4caf50',
        'foreground': 'white',
        'activebackground': '#45a049',
        'activeforeground': 'white',
        'font': ('Arial', 12),
        'bd': 1,  # Border width
        'relief': tk.SOLID  # Border style
    }
    reverse_button = tk.Button(toolbar,
                           text='Reverse x', command=reverse_axes, **button_style)
    reverse_button.pack(side=tk.LEFT)
    log_button = tk.Button(toolbar,
                           text='y to log', command=log_axes, **button_style)
    log_button.pack(side=tk.LEFT)
    reset_button = tk.Button(toolbar,
                             text='Reset axis', command=reset_axes, **button_style)
    reset_button.pack(side=tk.LEFT)

    button_style_2 = {
        'background': "#FFC0CB",
        'foreground': 'white',
        'activebackground': '#45a049',
        'activeforeground': 'white',
        'font': ('Arial', 12),
        'bd': 1,  # Border width
        'relief': tk.SOLID  # Border style
    }
    save_image_button = tk.Button(toolbar,
                            text='Save HR image', command=save_high_res_image, **button_style_2)
    save_image_button.pack(side=tk.LEFT)
    save_button = tk.Button(toolbar,
                            text='Save current data', command=save_text, **button_style_2)
    save_button.pack(side=tk.LEFT)

    

    return fig


def plotterTest():
    # Generate sample data
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    std = 0.2*y  # Standard deviation of y

    return plotter([[x, y, std], [x, y, 0]], ["X", "Y"], "Foo graph", size="double_size")


def plotterTest2():
    # Generate sample data
    path = 'D:\\OneDrive - UVa\\Program_devs\\RamanSpectra\\Dev_0\\Test data\\Si ref\\OptronLab\\AutoCalibration - Laser 785nm_Edge - Grating 600 (500nm) - Offset -238  PASS.txt'
    x = Raman_dataloader.load_spectra_data(path)
    info = Raman_dataloader.load_spectra_info(path)
    print(x)
    return plotter([[x[:, 0], x[:, 1]]],
                   [info['AxisType[1]']+" ("+info['AxisUnit[1]']+")",
                    "Wavenumber"+"("+info['AxisUnit[0]']+")"],
                   info['Title'],
                   lines=True)


def test_3():
    # Sample data
    x = np.linspace(0, 10, 100)
    y = np.sin(x)

    # Create an arrow data list
    arrow_data = [
        (30, 0.5),  # Arrow at index 30 with value 0.5
        (30, -0.8)  # Arrow at index 70 with value -0.8
    ]

    # Call the plotter function with arrow data
    fig, ax = plotter([(x, y)], ['x', 'y'], 'Sample Plot',
                      arrow=arrow_data, leyend_frame=[False, False])

    # Show the figure
    plt.show()

# test_3()
