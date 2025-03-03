# -*- coding: utf-8 -*-
"""
This is a GUI to load and plot the data for a single peak analysis

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

# Import FigureCanvasTkAgg class and NavigationToolbar2Tk from matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk  # Import the Tkinter module for GUI
from tkinter import ttk, filedialog, font
import numpy as np#from numpy import arange
from tkinter import messagebox
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')  # Set the backend to TkAgg

plt.ioff()  # To spyder not showing two windows

# Raman_single_GUI
from bin import Raman_fit
from bin import Raman_lateralpanel
from bin import Raman_dataloader
from bin import Raman_plot  # Import the FDTR_plot module for data plotting


def check_key(dictionary, key):
    """
    Check if a key exists in a dictionary.

    Parameters:
        dictionary (dict): The dictionary to check.
        key (any): The key to check for existence.

    Returns:
        bool: True if the key exists in the dictionary, False otherwise.

    Example:
        my_dict = {'a': 1, 'b': 2, 'c': 3}
        result = check_key(my_dict, 'b')
        print(result)  # Output: True
    """
    return key in dictionary

def error(text):
    """
    Display an error message based on the given type.

    Args:
        text (str): The message to show.


    """
    messagebox.showerror("Error", text)

def main(window_parent):
    # Global variables:
    global model_data

    filepath = "None"
    get_fields_data = []
    variation_data = []
    graph_type = ["linear", "linear"]
    foo_x = np.arange(100)
    foo_y = Raman_fit.lorentz_f(foo_x, 10, 5, 1)+Raman_fit.lorentz_f(
        foo_x, 20, 10, 2)+Raman_fit.lorentz_f(foo_x, 50, 10, 10)
    fig, ax = Raman_plot.plotter([[foo_x, foo_y]],
                                 ["Wavenumber (1/cm)",
                                  "Intensity (A.U.)"],
                                 "Load Data",
                                 leyends=['Example data'],
                                 res=150)
    selection_prev = 1
    # Area to create functions for GUI

    def buttonCall(path,time):
        """
        Callback function for the button to load data.

        Parameters:
            path (str): The path to the data file.
        """
        global ax, fig, model_data
        data_type = model_data.get()
        try:
            x, key = Raman_dataloader.load_spectra_data(path, data_type)
            info, key = Raman_dataloader.load_spectra_info(path, data_type)
            if key:
                  dat = [[x[:, 0], x[:, 1]/time]]
            fig, ax = Raman_plot.plotter(dat,
                                        ["Wavenumber (1/cm)",
                                        "Intensity (A.U.)"],
                                        info['Title'],
                                        leyends=['Raw data'],
                                        res=150,
                                        lines=True)
            fig = Raman_plot.update_plot(canvas, canvas_panel, fig, ax, dat)
            
        except:
            "key missing"

      

    def metadata(path, window):
        """
        Displays the metadata of the loaded data file in a popup window.

        Parameters:
            path (str): The path to the data file.
            window (Tk): The parent window.
        """
        global model_data
        data_type = model_data.get()
        popup_window = tk.Toplevel(window)
        popup_window.title("Experiment metadata")

        info, key = Raman_dataloader.load_spectra_info(path, data_type)

        text_widget = tk.Text(popup_window, padx=20, pady=5)
        text_widget.tag_configure("bold", font=("TkDefaultFont", 11, "bold"))
        text_widget.tag_configure("red", foreground="red")

        for key, value in info.items():
            text_widget.insert("end", f"{key}: ", "bold")
            text_widget.insert("end", value + "\n", "red")

        text_widget.pack(fill="both", expand=True)

    def button_clicked():
        """
        Callback function for the button to retrieve experiment data.
        """
        global filepath
        # Call the metadata function when the button is clicked
        metadata(filepath, window)

    def open_file():
        """
        Opens a file dialog to select and load the data file.
        """
        global filepath, ax, fig, model_data
        data_type = model_data.get()

        filepath = filedialog.askopenfilename(
            filetypes=[("TXT Files", "*.txt"), ("Dat Files","*.dat")])

        if filepath != '':
            # Call the lateral panel:
            # default values for the global variables
            buttonCall(filepath,1.0)
            raw_dat, key = Raman_dataloader.load_spectra_data(filepath, data_type)
            info, key2 = Raman_dataloader.load_spectra_info(filepath, data_type)
            if key and key2:
                # Call the buttonCall() function or perform any desired action
              
                lateral_panel = Raman_lateralpanel.create_lateral_panel(
                    canvas, canvas_panel, window, filepath, fig, raw_dat, info, data_type)
                lateral_panel.grid(row=0, column=1, rowspan=1, sticky="nsew")
                button_meta.config(state="normal")
                apply_button.config(state="normal")
           
    def apply_time():
        global filepath,ax, fig, model_data
        data_type = model_data.get()

        try:
            time=float(time_var.get())
            if time <= 0:
                raise ValueError("The input must be a number greater than 0")
           
            # Call the lateral panel:
            # default values for the global variables
            buttonCall(filepath,time)
            raw_dat, key = Raman_dataloader.load_spectra_data(filepath, data_type)  
                
            raw_dat[:,1]=raw_dat[:,1]/time
            info, key2 = Raman_dataloader.load_spectra_info(filepath, data_type)
            if key and key2:
                # Call the buttonCall() function or perform any desired action
                
                lateral_panel = Raman_lateralpanel.create_lateral_panel(
                    canvas, canvas_panel, window, filepath, fig, raw_dat, info, data_type, time_norm=time)
                lateral_panel.grid(row=0, column=1, rowspan=1, sticky="nsew")
                button_meta.config(state="normal")
        except:
            error("Time must be positive and >0")

    def on_popup_close(popup):
        """
        Callback function for the window close event.
        """
        figures = plt.get_fignums()

        for fig_num in figures:
            fig = plt.figure(fig_num)
            plt.close(fig)
        popup.grab_release()  # Release the grab
        window_parent.deiconify()  # Minimize the first window

        popup.destroy()

    # Area to create the main window
    window = tk.Toplevel(window_parent)
    window.title('Single_Raman_Peak_Processing_V1.0_JAC')
    window.geometry("1600x900")

    # minimize the launcher:
    window_parent.iconify()  # Minimize the first window

    window.rowconfigure(0, weight=1)
    window.columnconfigure(0, weight=1)
    window.grab_set()  # Grab the focus for the popup window
    window.protocol("WM_DELETE_WINDOW", lambda: on_popup_close(window))

    style = ttk.Style(window)
    style.configure("TButton", font=("Arial", 12))

    canvas_panel = ttk.Frame(window)
    canvas_panel.grid(row=0, column=0, sticky="nsew")
    canvas_panel.rowconfigure(0, weight=1)
    canvas_panel.columnconfigure(0, weight=1)

    # Area for the figure:
    # Create the initial graph
    canvas = FigureCanvasTkAgg(fig, master=canvas_panel)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
    # toolbar = NavigationToolbar2Tk(canvas, canvas_panel,pack_toolbar=False)
    # toolbar.update()
    # toolbar.grid(row=1, column=0, sticky="ew")

    # Area to create the buttons
    button_panel = ttk.Frame(window, borderwidth=2, relief="groove")
    button_panel.grid(row=2, column=0, sticky="ew", padx=10, pady=5)

    # Load data button
    # First select the data format:
    # Add the combo box to select the fitting method
    label_model = ttk.Label(button_panel, text="Select the data format:")
    label_model .grid(row=0, column=0, padx=5, pady=5)

    options = ['Horiba', 'B&Wtek', 'Brukker IR']

    combobox = ttk.Combobox(button_panel, values=options)
    combobox.set(options[0])
    model_data = combobox
    combobox.grid(row=0, column=1, padx=5, pady=5)
    # Create style
    # Create a style
    style = ttk.Style(window)
    style.configure('TButton', 
                background='#ADD8E6', 
                foreground='black',  # Black text color
                font=('Helvetica', 18, 'bold'))  # Bold font 

    data_button = ttk.Button(
        master=button_panel,
        command=open_file,        
        text="Load data",
        style='TButton')
    data_button.grid(row=0, column=2, padx=10, pady=5)
# Create a style
    style = ttk.Style(window)
    style.configure('TButton2.TButton', 
                background='#ADD8E6', 
                foreground='black',  # Black text color
                font=('Helvetica', 14))  # Bold font 
    button_meta = ttk.Button(master=button_panel, 
                            text='Metadata',
                            command=button_clicked,                            
                            style='TButton2.TButton',
                            state="disabled")
    button_meta.grid(row=0, column=3, padx=10, pady=5)
    # Remove the buttons from the window before using grid

    # Create a StringVar to hold the time value
    time_var = tk.StringVar()
    time_var.set("1")  # Set the default time value to 1

    # Create an Entry widget for the time input
    time_label = tk.Label(button_panel, text="Enter adquisition time value (s):")
    time_label.grid(row=0, column=4, padx=10, pady=5)
    time_entry = tk.Entry(button_panel, textvariable=time_var)
    time_entry.grid(row=0, column=5, padx=10, pady=5)

    apply_button = tk.Button(button_panel, text="Apply", command=apply_time,state="disabled")
    apply_button.grid(row=0, column=6, padx=10, pady=5)

    window.mainloop()
