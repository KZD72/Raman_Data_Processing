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
from tkinter import ttk, filedialog
from tkinter import messagebox
import os
import scipy.integrate as spi
import csv
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

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
    keys_dict = {key: key for key in first_nested_dict.keys() if key != 'Peak Model'}

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
    for key, value in dict2.items():
        # Iterate over each key-value pair in the 'fit_results' dictionary
        # Flag to indicate if a peak has been found
        peak_found = False
        
        for fit_key, fit_value in value['fit_results'].items():
            # Check if the 'Center' value is within the window of dict1_value
          
            if abs(float(fit_value['Center']) - float(dict1_value)) <= window:
                if not peak_found:
                    # If a peak has not been found yet, add the value of key_to_recover to the recovered values
                    recovered_values[key] = fit_value.get(key_to_recover)
                    peak_found = True
                    
                else:
                    # If a peak has already been found, issue a warning
                    error(f"Warning: More than one peak detected in the window for key {key}.")
                    recovered_values[key] = np.nan
                    
            
            if not key in recovered_values: 
                recovered_values[key] = np.nan
                

    return recovered_values, peak_found

def create_range(range_type, start, end, n):
    """
    Creates a range of values based on the specified type, start value, step size, and number of steps.

    Parameters:
    range_type (str): The type of range to create. Options are 'Linear', 'Logarithmic (base 10)', 'Logarithmic (base e)', and 'Exponential'.
    start (float): The start value of the range.
    step (float): The step size between values in the range.
    n (int): The number of steps to generate in the range.

    Returns:
    numpy.ndarray: A numpy array containing the generated range of values.
    """
    key=True
   
    if end<=start:
        key=False
    
    if range_type == 'Linear':
        try:
    
            values = np.linspace(start, end, num=n)
        except:
            key=False
            values=np.array([])
            error("Add a valid range")

    elif range_type == 'Log (base 10)':
        try:            
            values = np.logspace(np.log10(start), np.log10(end), num=n, base=10)
            print("log 10 val")
            print(values)
        except:
            key=False
            values=np.array([])
            error("Add a valid range")
        
    return values, key

def get_nan_indices(y_list):
    """
    Get the indices of NaN values in a list.

    Parameters:
    y_list (list): The list to check for NaN values.

    Returns:
    list: A list of indices where y_list is NaN.
    """
    # Convert the list to a numpy array
    y_array = np.array(y_list)

    # Get a boolean array that is True where y_array is NaN
    is_nan = np.isnan(y_array)

    # Get the indices where y_array is NaN
    nan_indices = np.where(is_nan)[0]

    return nan_indices.tolist()

###############################################################################
def create_dashboard(main_window, canvas, canvas_panel, dictionary):

    global selected_option, peaks,peaks_compare
    global peak_1_sel,peak_2_sel,peak_params,peak_params_sel
    global peak_1_sel_inner,peak_3_sel,window_width
    global custom_list_type
    global custom_list_type, start_entry, step_entry
    
    

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
                peak_1_sel = create_dropdown(box_frame3, peaks,1,0)
                peak_params_sel = create_dropdown(box_frame3, peak_params,1,1)    
                peak_2_sel = create_dropdown(box_frame3, peaks_compare,1,2)

            else:
                error("Field must be a valid number [1,100]")
        except:
            error("Field must be a valid number [1,100]")

    def show_inner_plot():
        global peaks,peaks_compare,peak_1_sel,peak_2_sel,peak_params,peak_params_sel
        global peak_1_sel_inner,peak_3_sel,peak_2_sel

        def update_inner_plot():
            if isinstance(x_label_entry.get(), str) and isinstance(y_label_entry.get(), str) and isinstance(title_entry.get(), str) :      
            
                fig, ax = Raman_plot.plotter(plots_to_show,
                                                [x_label_entry.get(),
                                                y_label_entry.get()],
                                                title_entry.get(),                                            
                                                lines=False,
                                                res=150,
                                                size="",
                                                leyends=peak_params_sel.get(),
                                                leyend_frame=[True, 'b']
                                                )
                Raman_plot.update_plot(new_canvas, new_canvas_panel, fig, ax, plots_to_show)
            else:
                error("Labels must be strings")
        
        
        inner_x, flag2=recover_values(peaks[peak_3_sel.get()], dictionary, float(window_width.get()), peak_1_sel_inner.get())
        inner_y, flag=recover_values(peaks[peak_1_sel.get()], dictionary, float(window_width.get()), peak_params_sel.get())
      
        if peak_2_sel.get()=="None":
            inner_y_norm={key: 1 for key in inner_y.keys()}
        else:
            inner_y_norm,flag3=recover_values(peaks[peak_2_sel.get()], dictionary, float(window_width.get()), peak_params_sel.get())
          
        
        inner_x_list=list(inner_x.values())
        nan_indices_x = get_nan_indices(inner_x_list)
        inner_y_list=list(inner_y.values())
        nan_indices_y = get_nan_indices(inner_y_list)
        inner_y_norm_list=list(inner_y_norm.values())
        nan_indices_y_norm = get_nan_indices(inner_y_norm_list)
       
        
        filtered_y_list = [inner_y_list[i]/inner_y_norm_list[i]
                            for i in range(len(inner_y_list))
                              if i not in nan_indices_y
                                and i not in nan_indices_x
                                  and i not in nan_indices_y_norm]
        filtered_x_list = [inner_x_list[i] 
                           for i in range(len(inner_x_list))
                               if i not in nan_indices_y
                                and i not in nan_indices_x
                                  and i not in nan_indices_y_norm]

        plots_to_show=[[filtered_x_list,filtered_y_list]]

        fig, ax = Raman_plot.plotter(plots_to_show,
                                            [peak_params_sel.get(),
                                             peak_1_sel_inner.get()],
                                            peak_params_sel.get(),                                            
                                            lines=False,
                                            res=150,
                                            leyends=peak_params_sel.get(),
                                            leyend_frame=[True, 'b']
                                            )
        # Create a new Toplevel window
        window = tk.Toplevel(main_window)
        window.title('Dahsboard graph')
                
        window.resizable(False, False)  # Disable resizing

        # Create a new canvas and add it to the new window
        new_canvas_panel = ttk.Frame(window)
        new_canvas_panel.grid(row=0, column=0, sticky="nsew")
        new_canvas_panel.rowconfigure(0, weight=1)
        new_canvas_panel.columnconfigure(0, weight=1)

        
        # Area for the figure:
        # Create the initial graph
        new_canvas = FigureCanvasTkAgg(fig, master=new_canvas_panel)
        new_canvas.draw()
        new_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        
        # Update the plot on the new canvas
        Raman_plot.update_plot(new_canvas, new_canvas_panel, fig, ax, plots_to_show)

            ## Area to introduce the custom inputs:
        # Create a panel on the right side of the canvas
        panel = ttk.Frame(window)
        panel.grid(row=0, column=1, sticky="nsew")

        # Create three entry fields in the panel
        title_entry = ttk.Entry(panel)
        title_entry.insert(0, "Enter title (LaTex allowed)")
        title_entry.grid(row=0, column=0, sticky="nsew")

        x_label_entry = ttk.Entry(panel)
        x_label_entry.insert(0, "Enter x-axis label (LaTex allowed)")
        x_label_entry.grid(row=1, column=0, sticky="nsew")

        y_label_entry = ttk.Entry(panel)
        y_label_entry.insert(0, "Enter y-axis label (LaTex allowed)")
        y_label_entry.grid(row=2, column=0, sticky="nsew")

        update_plot = tk.Button(panel, text="Update", command=update_inner_plot)
        update_plot.grid(row=3,column=0,sticky='ne')
    
    def show_outer_plot(option=1):
     
        global custom_list_type, start_entry, step_entry
        global peaks,peak_1_sel,window_width,peak_params_sel,peak_2_sel

        def update_inner_plot():
                if isinstance(x_label_entry.get(), str) and isinstance(y_label_entry.get(), str) and isinstance(title_entry.get(), str) :     
                    
                    try:
                
                        fig, ax = Raman_plot.plotter(plots_to_show,
                                                        [r'$' + x_label_entry.get() + '$',
                                                        r'$' + y_label_entry.get() + '$'],
                                                        r'$' +title_entry.get()+ '$',                                            
                                                        lines=False,
                                                        res=150,
                                                        size="",
                                                        leyends=peak_params_sel.get(),
                                                        leyend_frame=[True, 'b']
                                                        )
                        Raman_plot.update_plot(new_canvas, new_canvas_panel, fig, ax, plots_to_show)
                    except:
                        error("Error in LaTex format")    
                else:
                    error("Labels must be strings")

        key=False
        inner_y, flag=recover_values(peaks[peak_1_sel.get()], dictionary, float(window_width.get()), peak_params_sel.get())
        if option==1:
            x_list,key=create_range(custom_list_type.get(), float(start_entry.get()), float(step_entry.get()), len(inner_y))
        else:
            filepath = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv"), (
            "TXT Files", "*.txt")])
            if filepath != '':
                data = np.genfromtxt(filepath, dtype=float, delimiter=' ')
                print(data.shape)
            if data.ndim == 1 and len(data) == len(inner_y) and not np.any(np.isnan(data)):
                key=True
                x_list=data
            else:
                error("Imported file with wrong dimensions or non numerical entries")
    
       
        if key:
            
            if peak_2_sel.get()=="None":
                inner_y_norm={key: 1 for key in inner_y.keys()}
            else:
                inner_y_norm,flag3=recover_values(peaks[peak_2_sel.get()], dictionary, float(window_width.get()), peak_params_sel.get())
            
            
            inner_y_list=list(inner_y.values())
            nan_indices_y = get_nan_indices(inner_y_list)
            inner_y_norm_list=list(inner_y_norm.values())
            nan_indices_y_norm = get_nan_indices(inner_y_norm_list)
        
            
        
            
            filtered_y_list = [inner_y_list[i]/inner_y_norm_list[i]
                                for i in range(len(inner_y_list))
                                if i not in nan_indices_y                               
                                    and i not in nan_indices_y_norm]
            
            filtered_x_list = [x_list[i]
                                for i in range(len(x_list))
                                if i not in nan_indices_y                               
                                    and i not in nan_indices_y_norm]
            

            plots_to_show=[[filtered_x_list,filtered_y_list]]

            fig, ax = Raman_plot.plotter(plots_to_show,
                                                [peak_params_sel.get(),
                                                peak_1_sel_inner.get()],
                                                peak_params_sel.get(),                                            
                                                lines=False,
                                                res=150,
                                                size="",
                                                leyends=peak_params_sel.get(),
                                                leyend_frame=[True, 'b']
                                                )
            
             # Create a new Toplevel window
            window = tk.Toplevel(main_window)
            window.title('Dahsboard graph')
                 
            window.resizable(False, False)  # Disable resizing

            # Create a new canvas and add it to the new window
            new_canvas_panel = ttk.Frame(window)
            new_canvas_panel.grid(row=0, column=0, sticky="nsew")
            new_canvas_panel.rowconfigure(0, weight=1)
            new_canvas_panel.columnconfigure(0, weight=1)

           
            # Area for the figure:
            # Create the initial graph
            new_canvas = FigureCanvasTkAgg(fig, master=new_canvas_panel)
            new_canvas.draw()
            new_canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
           
            # Update the plot on the new canvas
            Raman_plot.update_plot(new_canvas, new_canvas_panel, fig, ax, plots_to_show)

             ## Area to introduce the custom inputs:
            # Create a panel on the right side of the canvas
            panel = ttk.Frame(window)
            panel.grid(row=0, column=1, sticky="nsew")

            # Create three entry fields in the panel
            title_entry = ttk.Entry(panel)
            title_entry.insert(0, "Enter title")
            title_entry.grid(row=0, column=0, sticky="nsew")

            x_label_entry = ttk.Entry(panel)
            x_label_entry.insert(0, "Enter x-axis label")
            x_label_entry.grid(row=1, column=0, sticky="nsew")

            y_label_entry = ttk.Entry(panel)
            y_label_entry.insert(0, "Enter y-axis label")
            y_label_entry.grid(row=2, column=0, sticky="nsew")

            update_plot = tk.Button(panel, text="Update", command=update_inner_plot)
            update_plot.grid(row=3,column=0,sticky='ne')
            
     
    def show_outer_plot_1():
        show_outer_plot(option=1)      

    def show_outer_plot_2():
        show_outer_plot(option=2)


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

    box_frame2_m = ttk.Frame(subpanel2, borderwidth=1, relief="groove")
    box_frame2_m.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')
    tk.Label(box_frame2_m, text="Select the magnitude for the y axis:").grid(row=0)
    box_frame4_m = ttk.Frame(subpanel2, borderwidth=1, relief="groove")
    box_frame4_m.grid(row=1, column=0, padx=10, pady=10, sticky='nsew')
    box_frame4_m.grid_rowconfigure(0, weight=1)
    box_frame4_m.grid_rowconfigure(1, weight=1)
    box_frame4_m.grid_columnconfigure(0, weight=1)
    box_frame4_m.grid_columnconfigure(1, weight=1)
    tk.Label(box_frame4_m, text="Select the magnitude for the x axis:").grid(row=0,column=0)
       
    box_frame2 = ttk.Frame(box_frame2_m, borderwidth=1, relief="groove")
    box_frame2.grid(row=1, column=0, padx=10, pady=10, sticky='nsew')
    # Configure box_frame to expand in both directions
    box_frame2.grid_rowconfigure(0, weight=1)
    box_frame2.grid_columnconfigure(0, weight=1)

    box_frame3 = ttk.Frame(box_frame2_m, borderwidth=1, relief="groove")
    box_frame3.grid(row=1, column=1, padx=10, pady=10, sticky='nsew')
    # Configure box_frame to expand in both directions
    box_frame3.grid_rowconfigure(0, weight=1)
    box_frame3.grid_columnconfigure(0, weight=1)

      
    box_frame4 = ttk.Frame(box_frame4_m, borderwidth=1, relief="groove")
    box_frame4.grid(row=1, column=1, padx=10, pady=10, sticky='nsew')
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
    postpro_tab.add(tab1, text="External data")
    postpro_tab.add(tab2,text="Inner data")
    #postpro_tab.bind("<<NotebookTabChanged>>", clean_external)
    tab1.grid_columnconfigure(0, weight=1)
    tab1.grid_columnconfigure(1, weight=1)
    tab2.grid_columnconfigure(0, weight=1)
    tab2.grid_columnconfigure(1, weight=1)

    ## Tab_1
    tk.Label(tab1, text="Select type of data for x values:").grid(row=0,column=0)
    ### SubTab1    
    subtab1_notebook=ttk.Notebook(tab1)
    subtab1_notebook.grid(row=1,column=0)
    subtab1 = ttk.Frame(subtab1_notebook)
    subtab1_notebook.add(subtab1, text='Predefined x')
    box_frame_subtab1 = ttk.Frame(subtab1, borderwidth=1, relief="groove")
    box_frame_subtab1 .grid(row=0, column=0, padx=1, pady=1, sticky='nsew')
    box_frame_subtab2 = ttk.Frame(subtab1, borderwidth=1, relief="groove")
    box_frame_subtab2 .grid(row=0, column=1, padx=1, pady=1, sticky='nsew')
    custom_list_type= tk.StringVar(value='Linear')
    item_custom_list=0
    for text in ['Linear', 'Log (base 10)']:
        custom_list_type_var = tk.Radiobutton(box_frame_subtab1, text=text, variable=custom_list_type, value=text)
        custom_list_type_var.grid(row=item_custom_list,column=0,sticky='nw')
        item_custom_list+=1

    start_label = tk.Label(box_frame_subtab2, text="Start:")
    start_label.grid(row=0,column=0,sticky='nw')
    start_entry = tk.Entry(box_frame_subtab2)
    start_entry.grid(row=1,column=0,sticky='ne')
    start_entry.insert(0, '1')  # Set default start value to 1

    step_label = tk.Label(box_frame_subtab2, text="End:")
    step_label.grid(row=2,column=0,sticky='nw')
    step_entry = tk.Entry(box_frame_subtab2)
    step_entry.grid(row=3,column=0,sticky='ne')
    step_entry.insert(0, '10')  # Set default start value to 1

    outer_dat_button_1 = tk.Button(box_frame_subtab2, text="Create range plot", command=show_outer_plot_1)
    outer_dat_button_1.grid(row=4,column=0,sticky='ne')

    ### SubTab2
    subtab2 = ttk.Frame(subtab1_notebook)
    subtab1_notebook.add(subtab2, text='Upload data for x')
    outer_dat_button_2 = tk.Button(subtab2, text="Load data and plot", command=show_outer_plot_2)
    outer_dat_button_2.grid(row=0,column=1,sticky='e')

     ## Tab_2
    tk.Label(tab2, text="Inner Parameter:").grid(row=0,column=0)
    peak_3_sel = create_dropdown(tab2, peaks,1,0)
    tk.Label(tab2, text="Parameter:").grid(row=0,column=1)
    peak_1_sel_inner= create_dropdown(tab2, peak_params,1,1)    
    plot_inner = tk.Button(tab2, text='Update', command=show_inner_plot)
    plot_inner.grid(row=2, column=1, padx=10, pady=10)
   





 