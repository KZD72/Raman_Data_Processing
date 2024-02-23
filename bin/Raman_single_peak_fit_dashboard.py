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
    dropdown = ttk.Combobox(master, textvariable=var, values=labels)
    dropdown.grid(row=row, column=column, padx=10, pady=10, sticky='nsew')  # Place the dropdown at the specified row and column
    return var  # Return the variable so you can get the selected option later

###############################################################################
def create_dashboard(main_window, canvas, canvas_panel, dictionary):

    global selected_option
    
   

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

    # Creation of main panel elements
    main_panel = tk.Frame(main_window, bg='white')
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
    subpanel1.grid_rowconfigure(0, weight=1)  # Configure row weight
    subpanel1.grid_columnconfigure(0, weight=1)  # Configure column weight

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
        (0, 0), window=frame1, anchor="nw", tags="frame1")
    #Add text widge to show fit info:
    text_widget = tk.Text(frame1, padx=10, pady=10)
    text_widget.tag_configure("bold", font=("TkDefaultFont", 11, "bold"))
    text_widget.tag_configure("red", foreground="blue")
    text_widget.insert("end", "Fitting Info:\n", "bold")
    text_widget.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

    button_full_fit_info= tk.Button(container_frame, text='Show full fit info', command= fit_info_clicked, state="disabled")
    button_full_fit_info.grid(row=1, column=0, padx=10, pady=10,sticky='w')


    # Second child frame
 
    frame2 = ttk.Frame(box_frame, borderwidth=2, relief="groove")    
    frame2.grid(row=0, column=1, padx=10, pady=10, sticky='nsew')
    # label dropdown 
    label_drop = ttk.Label(frame2, text="Select the file:")
    label_drop.grid(row=0, column=0, padx=5, pady=5)
    # dropdown to select the file:
    selected_option = create_dropdown(frame2, dictionary,1,0)

    button_fit = tk.Button(frame2, text='Show fit', command=fit_display)
    button_fit.grid(row=2, column=0, padx=10, pady=10)
    ###############################################################
    ### Subpanel 2                                              ###
    ###############################################################
   
