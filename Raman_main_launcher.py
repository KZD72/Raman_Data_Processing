# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 10:05:19 2023

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

@author:  J.Anaya
"""

# Raman processing program launcher

import tkinter as tk

from tkinter import ttk
from PIL import ImageTk, Image
import os
from tendo import singleton
import multiprocessing
import sys

from bin import Raman_single_GUI
############################################################


def main_window():
    # Area to create the main window
    window = tk.Tk()
    window.title('Raman_Data_Processing_V6 JAC')
    window.geometry("500x500")
    window.resizable(False, False)  # Disable resizing
    window.attributes("-topmost", False)

    # Apply Material Design style
    style = ttk.Style()
    style.configure("TButton",
                    font=("Roboto", 14),
                    relief="flat",
                    background="#2196f3",
                    foreground="white",
                    width=10,
                    height=5,
                    anchor="center")
    style.map("TButton",
              background=[("active", "#1976d2")])

    def on_closing():

        window.quit()
        window.destroy()

    # Define the button click handlers
    def launch_program1():
        Raman_single_GUI.main(window)

    def launch_program2():
        print("disabled")

    def launch_program3():
        print("disabled")

    def launch_program4():
        print("disabled")

    def load_images():
        global button_bg_1, button_bg_2, button_bg_3, button_bg_4

        image1 = Image.open(image1_path)
        button_bg_1 = ImageTk.PhotoImage(image1)

        image2 = Image.open(image2_path)
        button_bg_2 = ImageTk.PhotoImage(image2)

        image3 = Image.open(image3_path)
        button_bg_3 = ImageTk.PhotoImage(image3)

        image4 = Image.open(image4_path)
        button_bg_4 = ImageTk.PhotoImage(image4)

    # Get the path to the image files
    current_dir = os.path.dirname(os.path.abspath(__file__))
    image1_path = os.path.join(current_dir, "Resources", "Icon2.PNG")
    image2_path = os.path.join(current_dir, "Resources", "Icon3.PNG")
    image3_path = os.path.join(current_dir, "Resources", "Icon4.PNG")
    image4_path = os.path.join(current_dir, "Resources", "Icon5.PNG")

    # Load the images
    button_bg_1 = None
    button_bg_2 = None
    button_bg_3 = None
    button_bg_4 = None

    # Create a frame to hold the buttons
    button_frame = ttk.Frame(window)
    button_frame.pack(expand=True)

    # Load the images
    image1 = Image.open(image1_path)
    button_bg_1 = ImageTk.PhotoImage(image1)
    image2 = Image.open(image2_path)
    button_bg_2 = ImageTk.PhotoImage(image2)
    image3 = Image.open(image3_path)
    button_bg_3 = ImageTk.PhotoImage(image3)
    image4 = Image.open(image4_path)
    button_bg_4 = ImageTk.PhotoImage(image4)

    # Create the buttons using grid layout
    button1 = ttk.Button(button_frame, text="",
                         command=launch_program1, compound="center")
    button1.image = button_bg_1  # Store the image object as an attribute
    button1.config(image=button1.image)
    button1.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

    button2 = ttk.Button(
        button_frame, text="", command=launch_program2, compound="center", state="disabled")
    button2.image = button_bg_2  # Store the image object as an attribute
    button2.config(image=button2.image)
    button2.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")

    button3 = ttk.Button(
        button_frame, text="", command=launch_program3, compound="center", state="disabled")
    button3.image = button_bg_3  # Store the image object as an attribute
    button3.config(image=button3.image)
    button3.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")

    button4 = ttk.Button(
        button_frame, text="", command=launch_program4, compound="center", state="disabled")
    button4.image = button_bg_4  # Store the image object as an attribute
    button4.config(image=button4.image)
    button4.grid(row=1, column=1, padx=10, pady=10, sticky="nsew")

    # Configure row and column weights
    button_frame.grid_rowconfigure(0, weight=1)
    button_frame.grid_rowconfigure(1, weight=1)
    button_frame.grid_columnconfigure(0, weight=1)
    button_frame.grid_columnconfigure(1, weight=1)

    # Configure row and column weights
    button_frame.grid_rowconfigure(0, weight=1)
    button_frame.grid_rowconfigure(1, weight=1)
    button_frame.grid_columnconfigure(0, weight=1)
    button_frame.grid_columnconfigure(1, weight=1)

    # Calculate the center position of the window
    window.update_idletasks()
    window_width = window.winfo_width()
    window_height = window.winfo_height()
    screen_width = window.winfo_screenwidth()
    screen_height = window.winfo_screenheight()
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2
    window.geometry(f"+{x}+{y}")  # Center the window on the screen

    window.protocol("WM_DELETE_WINDOW", on_closing)
    window.mainloop()


###############################################################################


if __name__ == "__main__":
    multiprocessing.freeze_support()
    # Allow a single instance running.
    try:
        me = singleton.SingleInstance()
    except singleton.SingleInstanceException:
        print("Another instance of this application is already running, kill the terminal please.")
        sys.exit(1)
    app = main_window()
    sys.exit(app)

###############################################################################
