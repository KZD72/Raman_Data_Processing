# -*- coding: utf-8 -*-
"""
Created on Thu May 25 10:37:26 2023

Library to import data

@author:J.Anaya

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
# Raman_dataloader.py


from numpy import array, empty_like, transpose, argsort, expand_dims
import tkinter as tk
import numexpr as ne

from tkinter import messagebox

def check_file_open(file_path):
    try:
        with open(file_path, 'r') as file:
            return True
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except IOError:
        print(f"Error opening file '{file_path}'.")

def error(type):

    if type==0:
        messagebox.showerror("Error", "The CSV must contain 2 or 3 columns")
    elif type==1:
        messagebox.showerror("Error", "The CSV contains non-numerical values")
def gen_error(e):
    messagebox.showerror("Error", e)

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

def robust_float(value):
    """
    Convert a string value to a floating-point number, handling different decimal separators and scientific notation.

    Parameters:
        value (str): The string value to convert.

    Returns:
        float: The converted floating-point number.

    Example:
        result = robust_float('1,234.56')
        print(result)  # Output: 1234.56
    """
    try:
        # First, try to convert the string directly to a float
        value = value.replace(',', '.')
        return float(value)
    except ValueError:
        # If that fails, try replacing commas with periods and evaluating
        result = ne.evaluate(value.replace(',', '.'))
        return result.item()

def reorder_data(arr,multiple=False, transpose=False, squeeze=False):
    """
    Reorders a 2D array based on the values in the first column in descending order.

    Parameters:
    arr (numpy.ndarray or list): The input 2D array to be reordered.

    Returns:
    numpy.ndarray: The reordered array where rows are sorted based on the values in the first column.
    """
    # Convert arr to a NumPy array if it's not already
    arr_inner = array(arr)
  
    if squeeze:
        # Remove single-dimensional entries from the shape of the array
        arr_inner = array(arr[0])
    sorted_array=empty_like(arr_inner)
    if multiple:
        for i in range(arr_inner.shape[0]):
            # Get the (2, m) slice
            slice_arr = array(arr[i])
            if transpose:
                slice_arr=transpose(slice_arr)

            # Get the indices that would sort the first column
            sorted_indices = argsort(slice_arr[:, 0])  # Note the [::-1] to reverse the order

            if transpose:
                out=transpose(slice_arr[sorted_indices])
            else:
                out=slice_arr[sorted_indices]

            # Use the sorted indices to reorder the slice
            sorted_array[i] = out
    else:
        if transpose:
            slice_arr=transpose(arr_inner)
        # Get the indices that would sort the first column
        sorted_indices = argsort(arr_inner[:, 0]) # Note the [::-1] to reverse the order
        # Use the sorted indices to reorder the entire array
        sorted_array = arr_inner[sorted_indices]
        if transpose:
            sorted_array=transpose(sorted_array)
   
    if squeeze:
        sorted_array=expand_dims(sorted_array, axis=0)

    return sorted_array

def reorder_data_descending(arr,multiple=False, transpose=False, squeeze=False):
    """
    Reorders a 2D array based on the values in the first column in descending order.

    Parameters:
    arr (numpy.ndarray or list): The input 2D array to be reordered.

    Returns:
    numpy.ndarray: The reordered array where rows are sorted based on the values in the first column.
    """
    # Convert arr to a NumPy array if it's not already
    arr_inner = array(arr)

    if squeeze:
        # Remove single-dimensional entries from the shape of the array
        arr_inner = array(arr[0])
    sorted_array=empty_like(arr_inner)
    if multiple:
        for i in range(arr_inner.shape[0]):
            # Get the (2, m) slice
            slice_arr = array(arr[i])
            if transpose:
                slice_arr=transpose(slice_arr)

            # Get the indices that would sort the first column
            sorted_indices = argsort(slice_arr[:, 0])[::-1]  # Note the [::-1] to reverse the order

            if transpose:
                out=transpose(slice_arr[sorted_indices])
            else:
                out=slice_arr[sorted_indices]

            # Use the sorted indices to reorder the slice
            sorted_array[i] = out
    else:
        if transpose:
            slice_arr=transpose(arr_inner)
        # Get the indices that would sort the first column
        sorted_indices = argsort(arr_inner[:, 0])[::-1]  # Note the [::-1] to reverse the order
        # Use the sorted indices to reorder the entire array
        sorted_array = arr_inner[sorted_indices]
        if transpose:
            sorted_array=transpose(sorted_array)
   
    if squeeze:
        sorted_array=expand_dims(sorted_array, axis=0)
    
    return sorted_array

def load_spectra_info(path, data_type):
    """
    Load spectral information from a file.

    Parameters:
        path (str): The path to the file.
        data_type (str): The type of spectral data.

    Returns:
        dict: The dictionary containing the spectral information.

    Example:
        info = load_spectra_info('spectra.txt', 'Horiba')
        print(info)
    """
    key=False
    if data_type == "Horiba":
        try:
           
            return load_spectra_info_horiba(path), True
           
        except Exception as e:
            messagebox.showerror("Error", "Metadata not in a known format")
        return [], key
    elif data_type == "B&Wtek":
        try:
            return load_spectra_info_bwtech(path), True           
        except Exception as e:
            messagebox.showerror("Error", "Metadata not in a known format")
            return [], key
    elif data_type == "Brukker IR":
        try:
            return load_spectra_info_brukkerIR(path), True
            key=True
        except Exception as e:
            messagebox.showerror("Error", "Metadata not in a known format")
            return [], key
   


def load_spectra_data(path, data_type, silent=True):
    """
    Load spectral data from a file.

    Parameters:
        path (str): The path to the file.
        data_type (str): The type of spectral data.

    Returns:
        ndarray: The array containing the spectral data.

    Example:
        data = load_spectra_data('spectra.txt', 'Horiba')
        print(data)
    """
    key=False
    if data_type == "Horiba":
        try:
            return load_spectra_data_horiba(path),True
        except Exception as e:
            if  silent:
                messagebox.showerror("Error", "Data not in a known format")
            return [], key

    elif data_type == "B&Wtek":
        try:
            return load_spectra_data_bwtech(path),True
        
        except Exception as e:
            if  silent:
                messagebox.showerror("Error", "Data not in a known format")
            return [], key
    elif data_type == "Brukker IR":
        try:
            return load_spectra_data_brukkerIR(path),True
            
        except Exception as e:
            messagebox.showerror("Error", "Data not in a known format")
        return [], key



def load_spectra_info_horiba(path):
    """
    Load Horiba spectral information from a file.

    Parameters:
        path (str): The path to the file.

    Returns:
        dict: The dictionary containing the spectral information.

    Example:
        info = load_spectra_info_horiba('spectra.txt')
        print(info)
    """
    keys = []
    entries = []

    with open(path, 'r', encoding='utf-8', errors='replace')  as file:
        for line in file:
            if line.startswith('#'):
                line = line.strip('#').strip()  # Remove '#' and whitespace
                columns = line.split("=")  # Split line into columns
                # Extract text and numeric data
                keys.append(columns[0])
                value = columns[1].split('\t')
                if len(value) > 1:
                    entries.append(value[1])
                else:
                    entries.append(value[0])

    dictionary = {key: value for key, value in zip(keys, entries)}
    if len(dictionary) == 0:
        dictionary["Title"] = path
        dictionary["Warning"] = "No metadata"
    return dictionary

def load_spectra_data_horiba(path):
    """
    Load Horiba spectral data from a file.

    Parameters:
        path (str): The path to the file.

    Returns:
        ndarray: The array containing the spectral data.

    Example:
        data = load_spectra_data_horiba('spectra.txt')
        print(data)
    """
    wavenumber = []
    counts = []
    check_file_open(path)

    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                columns = line.split("	")  # Split line into columns
                dat = columns[1].split('\n')
                # Extract text and numeric data
                wavenumber.append(robust_float(columns[0]))
                #Check if direction

                counts.append(robust_float(dat[0]))

    data = transpose(array((wavenumber, counts), dtype='float64'))
    data=reorder_data(data)
    return data

def load_spectra_data_bwtech(filename):
    """
    Read numerical data from a text file starting from the line started by "Raw data".

    Parameters:
        filename (str): The name of the text file.

    Returns:
        data (list): A list of lists containing the numerical data.

    Example:
        data = read_data_from_file("data.txt")
        print(data)
    """
    data = []


    with open(filename, 'r') as file:
        iterator=0
        headers=[]
        wavenumber = []
        counts = []

        for line in file:
            line = line.strip()

            c1= any("raw data" in s.lower() for s in line.split(';'))

            if c1 :
                iterator=iterator+1
                headers.append(line.split(';'))

                index_1=headers[0].index("Raman Shift")
                index_2=headers[0].index("Dark Subtracted #1")

            if iterator>0:
                values = line.split(';')
                x=values[index_1]
                y=values[index_2]
                c11=check_number_f(x)
                c12=check_number_f(x)
                if c11 and c12:
                    wavenumber.append(robust_float(x))
                    counts.append(robust_float(y))

    data =  transpose(array((wavenumber, counts), dtype='float64'))
    data=reorder_data(data)
    return data

def load_spectra_info_bwtech(filename):
    """
    Read the remaining file info from a text file and return it as a dictionary.

    Parameters:
        filename (str): The name of the text file.

    Returns:
        file_info (dict): A dictionary containing the file info.

    Example:
        file_info = load_spectra_info_bwtech("data.txt")
        print(file_info)
    """
    file_info = {}

    with open(filename, 'r', encoding='utf-8', errors='replace')  as file:
        start_reading = False
        for line in file:
            line = line.strip()
            if line.startswith("Pixel"):
                start_reading = True
                continue
            if not start_reading:
                key_value = line.split(';')
                if len(key_value) > 1:
                    key, value = key_value[0].strip(), key_value[1].strip()
                    capitalized_key = key.capitalize()  # Capitalize the first letter of the key
                    file_info[capitalized_key] = value
    if len(file_info) == 0:
        file_info["Title"] = filename
        file_info["Warning"] = "No metadata"
    return file_info


def load_spectra_data_brukkerIR(filename):
    """
    Read numerical data from a text file starting from the line started by "Raw data".

    Parameters:
        filename (str): The name of the text file.

    Returns:
        data (list): A list of lists containing the numerical data.

    Example:
        data = read_data_from_file("data.txt")
        print(data)
    """
    data = []


    with open(filename, 'r') as file:
       
        headers=[]
        x_values = []
        y_values = []

        for line in file:
            line = line.strip()

            # If the line starts with '#c', it's the X values
            if line.startswith('#c'):
                # Split the line into words, convert each word to a float, and store the result in x_values
                x_values = [float(word) for word in line[2:].split()]
            # If the line doesn't start with '#', it's a Y value
            elif not line.startswith('#'):
                # Convert the line to a float and store the result in y_values
                y_values.append(float(line))

    data =  transpose(array((x_values, y_values), dtype='float64'))
    data=reorder_data(data)
    return data

def load_spectra_info_brukkerIR(filename):
    """
    This function reads a .dat file and extracts the metadata stored in lines starting with '#'.

    Parameters:
    filename: The path to the .dat file.

    Returns:
    file_info: A dictionary containing the title of the file and other metadata. If no metadata is found, a warning is included.
    """
    # Initialize an empty dictionary to store the file information
    file_info = {}

    # Open the file
    with open(filename, 'r', encoding='utf-8', errors='replace')  as file:
        # Initialize an empty list to store the metadata
        metadata = []
        
        # Iterate over each line in the file
        for line in file:
            # Remove leading and trailing whitespace from the line
            line = line.strip()

            # If the line starts with '#', it's metadata
            if line.startswith('#'):
                # Add the line to the metadata list (excluding the '#' character)
                metadata.append(line[1:])
    
    # If no metadata was found
    if len(metadata) == 0:
        # Set the title to the filename and include a warning
        file_info["Title"] = filename
        file_info["Warning"] = "No metadata"
    else:
        # Set the title and other information from the metadata
        file_info["Title"] = metadata[1]
        file_info["Other"] = metadata[0]

    # Return the file information
    return file_info

# Define the function to save the text widget contents to a file
def save_text(text_to_save):
    file_path = tk.filedialog.asksaveasfile(defaultextension=".txt",
                                          filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
    if file_path is not None:
        with open(file_path.name, "w") as f:
            f.write(text_to_save)





def test():
    path=r"D:\OneDrive - UVa\Program_devs\RamanSpectra\Dev_0\Test data\PCLS_22_1_20 s_785nm_Edge_600 (500nm)_100x_200 µm_2% (0_78mW)_01.txt"
    # print(load_spectra_info_horiba(path))
    # print(load_spectra_data_horiba(path))
    foo=load_spectra_data_horiba(path)
    print(foo[:,0])

def test2():
    path=r"E:\OneDrive - UVa\Program_devs\RamanSpectra\Dev_0\Test data\1_2_1.txt"
    #print(load_spectra_info_bwtech(path))
    foo=load_spectra_data_bwtech(path)
    print(foo[:,0])
def test3():
    filename=r"C:\Users\Admin\Documents\GitHub\Raman_Data_Processing\Test data\Other\IR\Arbocel_BC_200.0.dat"
    foo=load_spectra_data_brukkerIR(filename)
    print(foo)
#test()

#test3()