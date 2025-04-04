#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plot
import numpy as np
import os

import streamlit as st


# In[6]:


def run_IsotonicWork(multiple_animals_folder_path,drop,mouseid_min,mouseid_max,start_cutoff=50, end_cutoff_set=[240,210], baseline_cutoff=45,jump_threshold_set=[0.025,0.004]):

    def IsotonicWork(file_path,end_cutoff_set,jump_threshold_set):
        df = pd.read_csv(file_path, delimiter='\t',skiprows=31, usecols=[0, 1, 2])
        data_array = df.to_numpy()

        clicker = 0

        while clicker<2:
            end_cutoff = np.array(end_cutoff_set)[int(clicker)]
            jump_threshold = np.array(jump_threshold_set)[int(clicker)]

            x = data_array[:,0][start_cutoff:end_cutoff]*0.001 # time, currently in seconds
            y = data_array[:,1][start_cutoff:end_cutoff]*0.5 # position, currently in millimeters
            z = data_array[:,2][start_cutoff:end_cutoff]*99.06 # force, currently in mN

            baseline_start_y = y[0:baseline_cutoff]
            mean_start_y = np.mean(baseline_start_y)
            y = -data_array[:,1][start_cutoff:end_cutoff]*0.5 + mean_start_y

            baseline_start_z = z[0:baseline_cutoff]
            mean_start_z = np.mean(baseline_start_z)
            z = data_array[:,2][start_cutoff:end_cutoff]*99.06 - mean_start_z

            deriv = np.gradient(y,x)
            inst_p = 10**(-6)* z * deriv # currently in Watts

            zero_crossings = np.where(np.diff(np.sign(inst_p)) != 0)[0]
            x_mid = (x[:-1] + x[1:]) / 2

            zero_diffs = np.diff(x_mid[zero_crossings])
            significant_jumps = np.where(zero_diffs > jump_threshold)[0]

            if clicker==1 and len(significant_jumps)==0: 
                clicker = clicker+ 1
            if clicker==0 and len(significant_jumps)==0: 
                clicker = clicker+ 1
            if clicker==1 and len(significant_jumps)>0: 
                clicker = clicker+ 1
            if clicker==0 and len(significant_jumps)>0: 
                clicker = clicker+ 2

        in_start = x_mid[zero_crossings][np.min(significant_jumps)]
        in_end = x_mid[zero_crossings][np.max(significant_jumps) + 1]

        index_start = int(1000*in_start-start_cutoff)
        index_end = int(1000*in_end-start_cutoff)

        total_work = np.trapezoid(inst_p[index_start:index_end],x[index_start:index_end])

        return (float(total_work))

    num_folders = sum(os.path.isdir(os.path.join(multiple_animals_folder_path, item)) for item in os.listdir(multiple_animals_folder_path))

    outputs = []
    for _ in range(num_folders):
        outputs.append([])

    x_coord = np.linspace(0,149,num=150)
    p = mouseid_min

    fig, ax = plot.subplots(figsize=(11, 8))

    for j in range(mouseid_min,mouseid_max+1): 
        number = str(j)
        if os.path.exists(multiple_animals_folder_path+"M"+number+drop): 
            q = np.copy(j-p)

            folder_path = multiple_animals_folder_path+"M"+number+drop
            excel_files = [f for f in os.listdir(folder_path) if f.endswith(".xlsx") or f.endswith(".xls")]
            if excel_files:
                file_path = os.path.join(folder_path, excel_files[0])
                excel_path = pd.read_excel(file_path, sheet_name=0, header=None)
                e=excel_path.iloc[6,1]*0.001 # mass(kg)
            else:
                print("No Excel files found in "+folder_path)

            for filename in os.listdir(multiple_animals_folder_path+"M"+number+drop):
                if len(filename) >= 22 and filename[18:21] == "149":
                    before_149 = filename.split("149")[0]
                    after_149 = filename.split("149")[1]

            for i in range(0,150): 
                act = IsotonicWork(multiple_animals_folder_path+"M"+number+drop+before_149+str(i)+after_149,end_cutoff_set,jump_threshold_set)/e
                outputs[q].append(act)

            ax.plot(x_coord, np.array(outputs[q]), label=f"Animal {number}")
            ax.set_xlabel('Contraction Index')
            ax.set_ylabel('Normalized Work (J/kg)')
            ax.set_title('Total Isotonic Work')
            ax.legend()
            ax.grid()
        else: 
            p = p+1
    # np.savetxt("name.csv", outputs, delimiter=",", fmt='%s')

    return fig

# e.g. drop = "_Fatigue_rate_analysis/"

st.title("Isotonic Work Analysis")

multiple_animals_folder_path = st.text_input("Folder Path:")
drop = st.text_input("Folder Suffix:", value="_Fatigue/")

mouseid_min = st.number_input("Min Mouse ID:", min_value=0, step=1)
mouseid_max = st.number_input("Max Mouse ID:", min_value=0, step=1)
start_cutoff = st.number_input("Start Cutoff:", min_value=0, value=50, step=1)
end_cutoff_1 = st.number_input("End Cutoff (Primary):", min_value=start_cutoff+1, value=230, step=1)
end_cutoff_2 = st.number_input("End Cutoff (Secondary):", min_value=start_cutoff+1, value=210, step=1)
baseline_cutoff = st.number_input("Baseline Cutoff:", min_value=0, value=45, step=1)
jump_threshold_1 = st.number_input("Jump Threshold (Primary):", value=0.025, step=0.001)
jump_threshold_2 = st.number_input("Jump Threshold (Secondary):", value=0.005, step=0.001)

if st.button("Run Analysis"):
    if os.path.exists(multiple_animals_folder_path):
        fig = run_IsotonicWork(multiple_animals_folder_path,drop,mouseid_min,mouseid_max,
                               start_cutoff,[end_cutoff_1,end_cutoff_2],baseline_cutoff,[jump_threshold_1,jump_threshold_2])
        st.pyplot(fig)
    else:
        st.error("Invalid folder path. Please check and try again.")

