#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plot
import numpy as np
import os
import zipfile
import re
from datetime import datetime

import streamlit as st

# 0410 version. 0408 is gone. 0407 is on Lab Archives.


# In[8]:


def run_MaxInstPower(folder_path, start_cutoff=50, end_cutoff=215, baseline_cutoff=45):

    def MaxInstPower(file_path):
        df = pd.read_csv(file_path, delimiter='\t',skiprows=32, usecols=[0, 1, 2])
        data_array = df.to_numpy()

        # truncate these based on "int_sg":
        x = data_array[:,0][start_cutoff:end_cutoff]*0.001 # time, currently in seconds
        y = data_array[:,1][start_cutoff:end_cutoff]*0.5 # position, currently in millimeters
        z = data_array[:,2][start_cutoff:end_cutoff]*99.06 # force, currently in mN

        baseline_start_y = y[0:baseline_cutoff] # edit this based on "int_sg"
        mean_start_y = np.mean(baseline_start_y)
        y = -data_array[:,1][start_cutoff:end_cutoff]*0.5 + mean_start_y

        baseline_start_z = z[0:baseline_cutoff] # edit this based on "int_sg"
        mean_start_z = np.mean(baseline_start_z)
        z = data_array[:,2][start_cutoff:end_cutoff]*99.06 - mean_start_z

        deriv = np.gradient(y,x)
        inst_p = 10**(-6)* z * deriv # currently in Watts

        max_power = np.max(inst_p)

        return (max_power)

    # Custom sort function to sort files naturally (alphabet first, then numerically)
    def natural_sort_key(s):
        return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

    # sorted_files = sorted(files, key=natural_sort_key)

    num_folders = len(os.listdir(folder_path))-1

    outputs = []
    for _ in range(num_folders):
        outputs.append([])

    x_coord = np.linspace(0,149,num=150)
    q=0
    excel_files = []

    fig, ax = plot.subplots(figsize=(11, 8))

    sorted_names = sorted(os.listdir(folder_path), key=natural_sort_key)

    # Loop through the contents of the folder
    for name in sorted_names:
        if name != ".DS_Store":
            print("processing "+name+" ...")
            excel_files = []
            mouse_files = os.path.join(folder_path, name)

            # Only look inside folders (skip .zip files or __MACOSX)
            if os.path.isdir(mouse_files):
                for f in os.listdir(mouse_files):
                    if f.lower().endswith(".xlsx") or f.lower().endswith(".xls"):
                        print("Excel file "+f+" found in "+mouse_files)
                        excel_files.append(os.path.join(mouse_files, f))
                        excel_path = pd.read_excel(excel_files[0], sheet_name=0, header=None)
                        e=excel_path.iloc[6,1]*0.001 # mass(kg)

            sorted_files = sorted(os.listdir(mouse_files), key=natural_sort_key)
            for f in sorted_files:
                if f.lower().endswith(".xlsx") or f.lower().endswith(".xls"):
                    print("skipped Excel file "+f+" in "+mouse_files)
                else: 
                    ddf_files = os.path.join(mouse_files, f)
                    act = MaxInstPower(ddf_files)/e
                    outputs[q].append(act)             
            q=q+1
        else: 
            print("skipped "+name)

    for index, result in enumerate(outputs):
        ax.plot(x_coord, np.array(result), label=f"Animal {index+1}")
    ax.set_xlabel('Contraction Index')
    ax.set_ylabel('Normalized Power (W/kg)')
    ax.set_title('Peak Power')
    ax.legend()
    ax.grid()

    return fig,outputs

st.title("Peak Power Analysis")

uploaded_zip = st.file_uploader("Upload a .zip file", type="zip")

folder_structure = """
The folder you upload should be in the format: 
    name_of_folder.zip
      |
       -> M2034
           |
            -> ___.ddf
            -> ...
            -> ___.xlsx
       -> M2035
           |
            -> ...
       -> M2036
           |
            -> ...
       -> ...

Every 'Mouse ID' folder should consist of '.ddf' contraction files 
    AND an excel datasheet with tissue information.
"""

st.code(folder_structure, language='text')

if uploaded_zip:
    # Define the folder where we will save the uploaded file
    upload_folder = "uploaded_files"

    # Create the folder if it doesn't exist
    if not os.path.exists(upload_folder):
        os.makedirs(upload_folder)

    # Save the uploaded zip file
    zip_file_path = os.path.join(upload_folder, uploaded_zip.name)

    with open(zip_file_path, "wb") as f:
        f.write(uploaded_zip.getbuffer())

    st.success(f"File {uploaded_zip.name} uploaded successfully!")

    # Unzip the file to a folder
    unzip_folder = os.path.join(upload_folder, uploaded_zip.name.split('.')[0])

    # with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        # zip_ref.extractall(unzip_folder)

    with zipfile.ZipFile(uploaded_zip, 'r') as zip_ref:
        for member in zip_ref.namelist():
            if '__MACOSX' not in member:
                zip_ref.extract(member, unzip_folder)

    unzip_folder = os.path.join(unzip_folder, uploaded_zip.name.split('.')[0])

    st.write(f"Contents of {unzip_folder}:")
    st.write(os.listdir(unzip_folder))  # Display contents

    # You can now loop through the contents of the folder
    # for filename in os.listdir(unzip_folder):
        # st.write(filename)

start_cutoff = st.number_input("Start Cutoff:", min_value=0, value=50, step=1)
end_cutoff = st.number_input("End Cutoff:", min_value=start_cutoff+1, value=215, step=1)
baseline_cutoff = st.number_input("Baseline Cutoff:", min_value=0, value=45, step=1)

parameter_instructions = """
(Generally, there's no need to adjust these parameters!)
"""

st.code(parameter_instructions, language='text')

if st.button("Run Analysis"):
    st.write("Calculating...")
    fig,csv_output = run_MaxInstPower(unzip_folder,start_cutoff,end_cutoff,baseline_cutoff)
    st.write("Graphing...")
    st.pyplot(fig)

now = datetime.now()
timestamp = now.strftime("%Y_%m_%d_%H_%M_%S")
csv_name=f"peak_power_{timestamp}.csv"

csv = csv_output.to_csv(index=False).encode('utf-8')
st.download_button(label="Download CSV",
        data=csv,
        file_name=csv_name,
        mime='text/csv')

