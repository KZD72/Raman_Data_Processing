# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 12:53:27 2023

This creates a smple GUI for the single peak fitting panel dashboard

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

# Raman_single_peak_fit_dashboard

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


def on_popup_close(popup):
    popup.grab_release()  # Release the grab
    popup.destroy()


def create_dropdown(master, dictionary, row, column):
    # Get the labels from the dictionary
    labels = list(dictionary.keys())

    # Create a Tkinter variable
    var = tk.StringVar(master)

    # Set the default option
    var.set(labels[0])

    # Create the dropdown menu using a Combobox
    # Set a fixed width for the combobox
    dropdown = ttk.Combobox(master, textvariable=var, values=labels, width=20)
    dropdown.grid(row=row, column=column, padx=10, pady=10, sticky='nsew')  # Place the dropdown at the specified row and column
    return var  # Return the variable so you can get the selected option later


class FilteredListbox:
    def __init__(self, master, dictionary, row, column):
        self.master = master
        self.row = row
        self.column = column

        # All of the items for the listbox.
        self.items = list(dictionary.keys())

        # The current filter. Setting it to None initially forces the first update.
        self.curr_filter = None

        # Create the filter label and entry box.
        tk.Label(master, text='Filter Search:').grid(row=row, column=column)
        self.filter_box = tk.Entry(master)
        self.filter_box.grid(row=row+1, column=column)

        # A listbox with scrollbars.
        tk.Label(master, text='Select the file:').grid(row=row+2, column=column)

        yscrollbar = tk.Scrollbar(master, orient='vertical')
        yscrollbar.grid(row=row+3, column=column+1, sticky='ns')

        xscrollbar = tk.Scrollbar(master, orient='horizontal')
        xscrollbar.grid(row=row+4, column=column, sticky='we')

        self.listbox = tk.Listbox(master)
        self.listbox.grid(row=row+3, column=column, sticky='nswe')

        yscrollbar.config(command=self.listbox.yview)
        xscrollbar.config(command=self.listbox.xview)

        # The initial update.
        self.on_tick()

    def on_tick(self):
        if self.filter_box.get() != self.curr_filter:
            # The contents of the filter box has changed.
            self.curr_filter = self.filter_box.get()

            # Refresh the listbox.
            self.listbox.delete(0, 'end')

            for item in self.items:
                if self.curr_filter in item:
                    self.listbox.insert('end', item)

        self.master.after(200, self.on_tick)

    def get(self):
        return self.listbox.get(tk.ACTIVE)
    
    

def extract_peaks(dictionary, window):
    """
    Extracts and groups peak data from a nested dictionary structure.
    
    The function first extracts 'Center' values from the 'fit_results' of each entry in the input dictionary.
    It then groups these 'Center' values into groups based on a given window.
    Each group is represented by the average of its values.
    
    Parameters:
    dictionary (dict): The input dictionary containing the 'fit_results'.
    window (int): The window size for grouping 'Center' values.
    
    Returns:
    dict: A dictionary where each key is a group label and each value is the average of the 'Center' values in that group.
    """
    
    new_dictionary = {}

    for key, value in dictionary.items():
        fit_results = value.get('fit_results', {})
        center_values = {k: round(v.get('Center', 0)) for k, v in fit_results.items()}
        new_dictionary[key] = center_values

    grouped_dictionary = {}
    group_counter = 1

    for key, value_dict in new_dictionary.items():
        for sub_key, center_value in value_dict.items():
            if not any(abs(center_value - val) <= window for group in grouped_dictionary.values() for val in group):
                grouped_dictionary[f'P{group_counter}'] = [center_value]
                group_counter += 1
            else:
                for group_key, group_values in grouped_dictionary.items():
                    if any(abs(center_value - val) <= window for val in group_values):
                        grouped_dictionary[group_key].append(center_value)
                        break

    # Calculate the average of each group
    averaged_dictionary = {key+"_"+f"{round(sum(values) / len(values))}": round(sum(values) / len(values)) for key, values in grouped_dictionary.items()}
    #print(averaged_dictionary)
    return averaged_dictionary

def extract_keys(dictionary):
    """
    Extracts the keys of the nested dictionaries inside 'fit_results' from the first entry of the input dictionary.
    
    Parameters:
    dictionary (dict): The input dictionary containing the 'fit_results'.
    
    Returns:
    dict: A dictionary where each key is a key in the nested dictionaries inside 'fit_results' of the first entry in the input dictionary, and the value is the same as the key.
    """
    
    # Get the first entry of the input dictionary
    first_entry = next(iter(dictionary.values()))
    
    # Get the 'fit_results' dictionary of the first entry
    fit_results = first_entry.get('fit_results', {})
    
    # Get the first nested dictionary inside 'fit_results'
    first_nested_dict = next(iter(fit_results.values()), {})
    
    # Create a dictionary where the keys are the keys of the first nested dictionary and the values are the same as the keys
    keys_dict = {key: key for key in first_nested_dict.keys()}

    return keys_dict

def recover_values(dict1_value, dict2, window, key_to_recover):
    """
    Recovers values from dict2 based on a given value and a given window.
    
    Parameters:
    dict1_value (float): The value to look for in the 'Center' values of dict2.
    dict2 (dict): The dictionary containing the 'fit_results' to recover values from.
    window (int): The window size for comparing 'Center' values.
    key_to_recover (str): The key whose value should be recovered from the 'fit_results' dictionary.
    
    Returns:
    dict: A dictionary where each key is a key from the 'fit_results' dictionary in dict2 that matches dict1_value within the given window, and each value is the value of key_to_recover.
    """
    
    # Create a new dictionary to store the recovered values
    recovered_values = {}

    # Iterate over each key-value pair in dict2
    #print(dict2.items())
    for key, value in dict2.items():
        print(key)
        # Iterate over each key-value pair in the 'fit_results' dictionary
        for fit_key, fit_value in value['fit_results'].items():
            # Check if the 'Center' value is within the window of dict1_value
            if abs(fit_value['Center'] - dict1_value) <= window:
                # If it is, add the value of key_to_recover to the recovered values
                recovered_values[key] = fit_value.get(key_to_recover)

    return recovered_values

###############################################################################
def create_dashboard(main_window, canvas, canvas_panel, dictionary):

    global selected_option, peaks,peaks_compare, peak_1_sel,peak_2_sel,peak_params,peak_params_sel
    
    

    def on_closing():
        """
        Handles the closing event of the main window.
        Performs necessary cleanup actions before closing the application.
        """
        main_window.quit()  # Quit the main window event loop
        main_window.destroy()  # Destroy the main window

   

    def fit_display():

        global selected_option

        plots_to_show=dictionary[selected_option.get()]["plots_to_show"]

        fig, ax = Raman_plot.plotter(plots_to_show,
                                         ["Wavenumber (1/cm)",
                                          "Intensity (A.U.)"],
                                         dictionary[selected_option.get()]["title"],
                                         leyends=dictionary[selected_option.get()]["leyends"],
                                         lines=True,
                                         res=150,
                                         leyend_frame=[True, 'b']
                                         )
        Raman_plot.update_plot(canvas, canvas_panel,
                                   fig, ax, plots_to_show)
        
        text_widget.delete("1.0", tk.END)

        for peak_number, parameters in dictionary[selected_option.get()]["fit_results"].items():
                text_widget.insert(
                    "end", f"Peak {peak_number}:" + "\n", "bold")
                for parameter, value in parameters.items():
                    text_widget.insert(
                        "end", f"{parameter}: {value}" + "\n", "bold")
                text_widget.insert("end", "\n", "bold")
        button_full_fit_info.config(state="normal")
        
        
    def fit_info_clicked():
            """
            Handles the event when the fit info button is clicked.
            Displays fit information in a popup window and provides an option to save the information.

            Returns:
                None
            """
            global selected_option

            # Create a new popup window
            popup_window = tk.Toplevel(main_window)
            popup_window.title("Full Fit info")
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

            # Display fit information
            text_widget.insert("end",dictionary[selected_option.get()]["full_fit_info"], "bold")

            text_widget.grid(row=0, column=0, sticky='nsew')

            # Create a new style for the button
            save_button_style = ttk.Style()
            save_button_style.configure(
                'Save.TButton', foreground='Blue', font=('Arial', 14))

            # Create the "Save" button
            save_button = ttk.Button(popup_window, text="Save", command=lambda: save_text(
                dictionary[selected_option.get()]["full_fit_info"]), style='Save.TButton')
            save_button.grid(row=1, column=0, sticky='nsew')

            popup_window.mainloop()
    
    def update_window():
        global peaks,peaks_compare, peak_params, peak_1_sel,peak_2_sel,peak_params_sel
        try:
            if float(window_width.get())>=1 and float(window_width.get())<100:
                peaks=extract_peaks(dictionary,float(window_width.get()))        
                peaks_compare={"None": 1}
                peaks_compare.update(peaks)
                peak_1_sel = create_dropdown(box_frame2, peaks,1,0)
                peak_params_sel = create_dropdown(box_frame2, peak_params,1,1)    
                peak_2_sel = create_dropdown(box_frame2, peaks_compare,1,2)
            else:
                error("Field must be a valid number [1,100]")
        except:
            error("Field must be a valid number [1,100]")

    def show_value():
        global peaks,peaks_compare,peak_1_sel,peak_2_sel,peak_params,peak_params_sel
    
        print(recover_values(peaks[peak_1_sel.get()], dictionary, float(window_width.get()), peak_params_sel.get()))


    ###############################################################
    ###############################################################
    ### Graphic crreation                                       ###
    ###############################################################
    ###############################################################       
    # Creation of main panel elements
    main_panel = tk.Frame(main_window, bg='white')
    main_panel.grid(row=0, column=0, padx=1, pady=1, sticky='nsew')
    
    # Grid layout configuration
    main_panel.grid_columnconfigure(0, weight=1)
    main_panel.grid_rowconfigure(0, weight=1)
    main_panel.grid_rowconfigure(1, weight=1)


    # Create subpanel frames within the main panel
    subpanel1 = tk.Frame(main_panel, bg='lightblue')
    subpanel2 = tk.Frame(main_panel, bg='lightgrey')   

    # Grid placement of subpanel frames
    subpanel1.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')
    # Subpanel 1
    subpanel1.grid_rowconfigure(0, weight=1)  # Configure row weight
    subpanel1.grid_columnconfigure(0, weight=1)  # Configure column weight
    #subpanel1.grid_columnconfigure(1, weight=1)  # Configure column weight

    subpanel2.grid(row=1, column=0, padx=10, pady=10, sticky='nsew')
    # Subpanel 2
    subpanel2.grid_rowconfigure(0, weight=1)  # Configure row weight
    subpanel2.grid_columnconfigure(0, weight=1)  # Configure column weight
    
    ###############################################################
    ### Subpanel 1                                              ###
    ###############################################################
    # Create the decoration of the subpanels
    # Create a style for the labels
    style = ttk.Style()
    style.configure("Box.TLabel", background=main_window["background"])

    box_frame = ttk.Frame(subpanel1, borderwidth=1)
    box_frame.grid(row=1, column=0, padx=10, pady=10, sticky='nsew')
    # Configure box_frame to expand in both directions
    box_frame.grid_rowconfigure(0, weight=1)
    box_frame.grid_columnconfigure(0, weight=1)
    box_frame.grid_columnconfigure(1, weight=1)

    # Create a new frame as a container for frame1 and the scrollbar
    container_frame = ttk.Frame(box_frame, borderwidth=1)
    container_frame.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

    # Configure container_frame to expand in both directions
    container_frame.grid_rowconfigure(0, weight=1)
    container_frame.grid_rowconfigure(1, weight=1)
    container_frame.grid_columnconfigure(0, weight=1)

    # Create a canvas widget inside the container frame
    canvas_1 = tk.Canvas(container_frame)
    canvas_1.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')
    # Configure canvas_peaks to expand in both directions
    canvas_1.grid_rowconfigure(0, weight=1)
    canvas_1.grid_columnconfigure(0, weight=1)
   
    # Create a scrollbar widget
    scrollbar = ttk.Scrollbar(
        container_frame, orient="vertical", command=canvas_1.yview)
    scrollbar.grid(row=0, column=1, sticky='ns')
    # Add this line to make the scrollbar appear shaded
    scrollbar.config(takefocus=0)

    # Configure the canvas to use the scrollbar
    canvas_1.configure(yscrollcommand=scrollbar.set)

    # Create a frame to hold the fields
    frame1 = ttk.Frame(canvas_1, borderwidth=2, relief="groove")
    # Configure canvas_peaks to expand in both directions
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=1)
    # Configure the canvas to adjust scroll region when the frame size changes
    frame1.bind("<Configure>", lambda event: canvas_1.configure(
        scrollregion=canvas_1.bbox("all")))
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_columnconfigure(0, weight=1)

    canvas_1.create_window(
        (0, 0), window=frame1, anchor="nw", tags="frame1",width=450)
    #Add text widge to show fit info:
    text_widget = tk.Text(frame1, padx=2, pady=2, wrap="none")
    text_widget.tag_configure("bold", font=("TkDefaultFont", 11, "bold"))
    text_widget.tag_configure("red", foreground="blue")
    text_widget.insert("end", "Fitting Info:\n", "bold")
    text_widget.grid(row=0, column=0, padx=2, pady=2, sticky='nsew')

    button_full_fit_info= tk.Button(container_frame, text='Show full fit info', command= fit_info_clicked, state="disabled")
    button_full_fit_info.grid(row=1, column=0, padx=10, pady=10,sticky='w')


    # Second child frame
 
    frame2 = ttk.Frame(box_frame, borderwidth=2, relief="groove")    
    frame2.grid(row=0, column=1, padx=10, pady=10, sticky='nsew')   
    frame2.grid_columnconfigure(0, weight=1)
    # dropdown to select the file:
    selected_option = FilteredListbox(frame2, dictionary, 0,0)
    button_fit = tk.Button(frame2, text='Show fit', command=fit_display)
    button_fit.grid(row=5, column=0, padx=10, pady=10)
    ###############################################################
    ### Subpanel 2                                              ###
    ###############################################################
   
    box_frame2 = ttk.Frame(subpanel2, borderwidth=1, relief="groove")
    box_frame2.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')
    # Configure box_frame to expand in both directions
    box_frame2.grid_rowconfigure(0, weight=1)
    box_frame2.grid_columnconfigure(0, weight=1)

    box_frame3 = ttk.Frame(subpanel2, borderwidth=1, relief="groove")
    box_frame3.grid(row=0, column=1, padx=10, pady=10, sticky='nsew')
    # Configure box_frame to expand in both directions
    box_frame3.grid_rowconfigure(0, weight=1)
    box_frame3.grid_columnconfigure(0, weight=1)

    
    box_frame4 = ttk.Frame(subpanel2, borderwidth=1, relief="groove")
    box_frame4.grid(row=1, column=0, padx=10, pady=10, sticky='nsew')
    # Configure box_frame to expand in both directions
    box_frame4.grid_rowconfigure(0, weight=1)
    box_frame4.grid_columnconfigure(0, weight=1)

    ### Box frame
    tk.Label(box_frame2, text="Spectral peak window").grid(row=0)
    window_width = tk.Entry(box_frame2)
    window_width.insert(0, "5.0")
    window_width.grid(row=1, column=0)

    peaks=extract_peaks(dictionary,float(window_width.get()))
    peaks_compare={"None": 1}
    peaks_compare.update(peaks)

    peak_params=extract_keys(dictionary)

    window_update = tk.Button(box_frame2, text='Update', command=update_window)
    window_update.grid(row=1, column=1, padx=10, pady=10)

     ### Box frame2
    tk.Label(box_frame3, text="Peak:").grid(row=0,column=0)
    peak_1_sel = create_dropdown(box_frame3, peaks,1,0)
    tk.Label(box_frame3, text="Parameter:").grid(row=0,column=1)
    peak_params_sel = create_dropdown(box_frame3, peak_params,1,1)    
    tk.Label(box_frame3, text="Normalise by:").grid(row=0,column=2)
    peak_2_sel = create_dropdown(box_frame3, peaks_compare,1,2)

    ### Box frame3
    postpro_tab = ttk.Notebook(box_frame4)
    postpro_tab.grid(row=0, column=0, padx=10, pady=2, sticky='nsew')
    tab1 = ttk.Frame(postpro_tab)
    tab2 = ttk.Frame(postpro_tab)
    postpro_tab.add(tab1, text="Inner data")
    postpro_tab.add(tab2, text="External data")
    #postpro_tab.bind("<<NotebookTabChanged>>", clean_external)
    tab1.grid_columnconfigure(0, weight=4)
    tab1.grid_columnconfigure(1, weight=2)
    tab2.grid_columnconfigure(3, weight=6)
    tab2.grid_columnconfigure(0, weight=4)

    ## Tab_1
    tk.Label(tab1, text="Inner Parameter").grid(row=0,column=0)
    peak_inner = create_dropdown(tab1, peaks,1,0)
    tk.Label(box_frame3, text="Parameter:").grid(row=0,column=1)
    peak_inner_param = create_dropdown(tab1, peak_params,1,1)    
    plot_inner = tk.Button(tab1, text='Update', command=show_value)
    plot_inner.grid(row=2, column=1, padx=10, pady=10)
   





 