"""
power.py

Description: Process all animals in a folder, compute normalized maximum instantaneous 
power for each contraction, and return a figure plus the raw results.
 [y] fixed contraction detector
 [y] add debug function
 [] write start & end detection into max_inst function, not parse
 [] make x axis of graph the file names
 [] add the option to save mass data
    add the option to upload mass data (from previous generation)
"""

__author__ = "Chaoran Yang"
__version__ = "2.1"
__email__ = "cy197@duke.edu"
__date__ = "2026-01-15"

import os
import re
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import streamlit as st
import zipfile
from datetime import datetime

import random

def natural_sort_key(s: str):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r"([0-9]+)", s)]

def parse_dmc_file(file_path: str) -> dict[str, np.ndarray]:
    """
    Parse a DMCv5.x-style data file.

    Extracts:
      - sample frequency
      - calibration info per channel (scale, offset, units)
      - data table (Sample, AI0, AI1, ...)

    Returns
    -------
    dict with keys:
      'time'        : time axis in seconds
      'length_mm'   : muscle length in mm (using calibrated AI channel)
      'force_mN'    : force in mN (using calibrated AI channel)
      'raw_df'      : full pandas DataFrame of the test data
    """
    
    with open(file_path, "r", encoding="latin-1") as f:
        lines = f.readlines()

    data_marker_idx = None
    sample_freq_hz = None

    # Parse header-level info, including calibration block
    channel_names: List[str] = []
    units: List[str] = []
    scales: List[float] = []
    offsets: List[float] = []

    i = 0
    prot_skiprows = None  # row index where "Protocol Array" line occurs

    while i < len(lines):
        line = lines[i].strip()

        # Sample frequency
        if line.startswith("Sample Frequency"):
            try:
                sample_freq_hz = float(line.split(":")[1].strip())
            except (IndexError, ValueError):
                raise ValueError(f"Could not parse sample frequency in {file_path}")
            i += 1
            continue

        # Protocol Array start
        if line.startswith("Protocol Array"):
            prot_skiprows = i
            i += 1
            continue

        # Locate data marker (end of protocol section)
        if line.startswith("Test Data in Volts"):
            data_marker_idx = i

            if prot_skiprows is None:
                raise ValueError(f'Found "Test Data in Volts" before "Protocol Array" in {file_path}')

            # Assume protocol table begins on the next line after "Protocol Array"
            protocol_start = prot_skiprows + 1
            protocol_nrows = data_marker_idx - protocol_start
            if protocol_nrows <= 0:
                raise ValueError(f"No protocol rows found in {file_path}")

            prot_df = pd.read_csv(
                file_path,
                delimiter="\t",
                skiprows=protocol_start,
                nrows=protocol_nrows,
                engine="python",
                header=None,
            )

            initial_baseline_end = None
            for j in range(len(prot_df)):
                if str(prot_df.iloc[j, 1]).strip() == "Stimulus-Tetanus":
                    initial_baseline_end = 0.8*(float(prot_df.iloc[j, 0]) + 
                                                float(prot_df.iloc[j, 3].split(",")[0].strip()))*sample_freq_hz
                    print("(cy) initial baseline end")
                    print(initial_baseline_end)

                    final_baseline_start = 1.2*(float(prot_df.iloc[j, 0]) + 
                                                float(prot_df.iloc[j, 3].split(",")[3].strip()))*sample_freq_hz
                    print("(cy) final baseline start")
                    print(final_baseline_start)
                    break

            break

        # Calibration block
        if line.startswith("Channel"):
            channel_names = lines[i].strip().split("\t")[1:]
            units = lines[i + 1].strip().split("\t")[1:]
            scales = [float(x) for x in lines[i + 2].strip().split("\t")[1:]]
            offsets = [float(x) for x in lines[i + 3].strip().split("\t")[1:]]
            i += 4
            continue

        i += 1

    if data_marker_idx is None:
        raise ValueError(f"'Test Data in Volts:' not found in {file_path}")

    if sample_freq_hz is None:
        raise ValueError(f"Sample frequency not found in {file_path}")

    # Build calibration dictionary
    calib: dict[str, dict[str, float]] = {}
    for name, unit, scale, offset in zip(channel_names, units, scales, offsets):
        calib[name] = {
            "units": unit,
            "scale": scale,
            "offset": offset,
        }

    # Read test data using pandas:
    # skip all lines up to and including 'Test Data in Volts:'
    # so the next line ('Sample\tAI0\tAI1...') becomes the header row.
    skiprows = data_marker_idx + 1
    df = pd.read_csv(
        file_path,
        delimiter="\t",
        skiprows=skiprows,
        engine="python"
    )

    # Identify channels:
    #   - Length: first AI channel with units 'mm'
    #   - Force : first AI channel with units 'Ref' and positive scale
    length_chan = None
    force_chan = None

    for name in calib:
        meta = calib[name]
        if name.startswith("AI") and meta["units"] == "mm" and length_chan is None:
            length_chan = name
        if name.startswith("AI") and meta["units"] == "Ref" and meta["scale"] > 0 and force_chan is None:
            force_chan = name

    if length_chan is None or force_chan is None:
        raise ValueError(
            f"Could not identify length/force channels in calibration table for {file_path}"
        )

    # Raw volts
    length_volts = df[length_chan].to_numpy(dtype=float)
    force_volts = df[force_chan].to_numpy(dtype=float)

    # Convert to physical units: units = (volts - offset) * scale
    length_mm = (length_volts - calib[length_chan]["offset"]) * calib[length_chan]["scale"]
    force_ref = (force_volts - calib[force_chan]["offset"]) * calib[force_chan]["scale"]

    # Treat 'Ref' for force as mN (matches your previous scaling)
    force_mN = force_ref

    # Build time axis
    # Sample column is an integer index; but we can compute time directly from index
    num_samples = len(df)
    sample_indices = np.arange(num_samples, dtype=float)
    time_s = sample_indices / sample_freq_hz

    return {
        "time": time_s,
        "length_mm": length_mm,
        "force_mN": force_mN,
        "raw_df": df,
        "start_idx": int(0.8*initial_baseline_end),
        "end_idx": int(1.2*final_baseline_start)
    }

def max_instantaneous_power_from_file(file_path: str) -> float:
    """
    Compute the maximum instantaneous power for a single file.

    Steps:
      - Parse DMC file (header + data).
      - Detect contraction window automatically from force trace.
      - Use a pre-contraction middle segment as baseline.
      - Compute instantaneous power: P(t) = F(t) * v(t), converted to Watts.
    """
    parsed = parse_dmc_file(file_path)
    time = parsed["time"]
    length_mm = parsed["length_mm"]
    force_mN = parsed["force_mN"]
    start_idx = parsed["start_idx"]
    end_idx = parsed["end_idx"]

    # Baseline correction using the pre-contraction middle segment
    length_baseline = float(np.mean(length_mm[0:start_idx]))
    force_baseline = float(np.mean(force_mN[0:start_idx]))

    length_seg = length_mm[start_idx:end_idx] - length_baseline
    force_seg = force_mN[start_idx:end_idx] - force_baseline
    time_seg = time[start_idx:end_idx]

    # Velocity (mm/s)
    velocity_mm_s = np.gradient(length_seg, time_seg)

    # Power: force (mN) * velocity (mm/s) -> W via 1e-6 factor
    inst_power_W = - 1e-6 * force_seg * velocity_mm_s
    min_idx = np.argmin(inst_power_W)
    inst_power_sliced = inst_power_W[:min_idx + 1]
    max_power_W = float(np.max(inst_power_sliced))

    return max_power_W

def run_max_inst_power(folder_path: str,
                       mass_kg: float) -> List[List[float]]:
    """
    Process one animal folder under `folder_path` and compute normalized
    peak instantaneous power (W/kg) for each contraction.

    Parameters
    ----------
    folder_path : str
        Root directory of the one animal folder.
    mass_kg : float
        Animal mass in kilograms used to normalize peak power.

    Returns
    -------
    outputs : list[list[float]]
        Normalized peak power values (W/kg) for each contraction and animal.
    """
    animal_results: List[float] = []

    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"Folder does not exist: {folder_path}")

    files = [
        f for f in os.listdir(folder_path)
        if f.lower().endswith(".ddf") and not f.startswith(".")
    ]

    for f in sorted(files, key=natural_sort_key):
        data_file = os.path.join(folder_path, f)

        print(f"Processing {f} ...")

        peak_power_W = max_instantaneous_power_from_file(data_file)
        normalized_power = peak_power_W / mass_kg  # W/kg
        animal_results.append(normalized_power)

    return animal_results

def val_max_inst_power(file_path: str,
                       i: float) -> plt:
    """
    Selects specific .ddf files and shows the power vs. time plot
    to validate run_max_inst_power results 

    Steps:
      - Parse DMC file (header + data).
      - Detect contraction window automatically from force trace.
      - Use a pre-contraction middle segment as baseline.
      - Compute instantaneous power: P(t) = F(t) * v(t), converted to Watts.
    """
    parsed = parse_dmc_file(file_path)
    time = parsed["time"]
    length_mm = parsed["length_mm"]
    force_mN = parsed["force_mN"]
    start_idx = parsed["start_idx"]
    end_idx = parsed["end_idx"]

    # Baseline correction using the pre-contraction middle segment
    length_baseline = float(np.mean(length_mm[0:start_idx]))
    force_baseline = float(np.mean(force_mN[0:start_idx]))

    length_seg = length_mm[start_idx:end_idx] - length_baseline
    force_seg = force_mN[start_idx:end_idx] - force_baseline
    time_seg = time[start_idx:end_idx]

    # Velocity (mm/s)
    velocity_mm_s = np.gradient(length_seg, time_seg)

    # Power: force (mN) * velocity (mm/s) -> W via 1e-6 factor
    inst_power_W = - 1e-6 * force_seg * velocity_mm_s
    min_idx = np.argmin(inst_power_W)
    inst_power_sliced = inst_power_W[:min_idx + 1]
    max_power_W = float(np.max(inst_power_sliced))
    max_idx = np.argmax(inst_power_sliced)

    fig_l, ax_l = plt.subplots(figsize=(11, 8))
    ax_l.plot(time_seg[:min_idx + 1], length_seg[:min_idx + 1])
    ax_l.set_xlabel("Time (s)")
    ax_l.set_ylabel("Length (mm))")
    ax_l.set_title(f"{file_path} (index {i}) Length")
    ax_l.grid(True)

    fig_f, ax_f = plt.subplots(figsize=(11, 8))
    ax_f.plot(time_seg[:min_idx + 1], force_seg[:min_idx + 1])
    ax_f.set_xlabel("Time (s)")
    ax_f.set_ylabel("Force (mN)")
    ax_f.set_title(f"{file_path} (index {i}) Force")
    ax_f.grid(True)

    fig_p, ax_p = plt.subplots(figsize=(11, 8))
    ax_p.plot(time_seg[:min_idx + 1], inst_power_W[:min_idx + 1], label="Power")
    ax_p.plot(time_seg[max_idx], max_power_W, marker="o", label="Max Power")
    ax_p.set_xlabel("Time (s)")
    ax_p.set_ylabel("Power (W)")
    ax_p.legend()
    ax_p.set_title(f"{file_path} (index {i}) Power")
    ax_p.grid(True)

    return fig_l, fig_f, fig_p

st.title("Peak Power Analysis")

uploaded_zip = st.file_uploader("Upload a .zip file", type="zip")

folder_structure = """
The folder you upload should be in the format: 
    name_of_folder.zip
      |
       -> animal_1
           |
            -> ___.ddf
            -> ...
       -> animal_2
           |
            -> ...
       -> animal_3
           |
            -> ...
       -> ...
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
    sorted_folder_names = sorted(os.listdir(unzip_folder), key=natural_sort_key)
    st.write(sorted_folder_names)  # Display contents

    mass_g: List[float] = []

    for i, filename in enumerate(os.listdir(unzip_folder)):
        value = st.number_input(f"{i} Mass (g):", min_value=0.0, value=1.0, step=0.00001) 
        mass_g.append(value)

if "analysis_done" not in st.session_state:
    st.session_state.analysis_done = False
if "csv_output" not in st.session_state:
    st.session_state.csv_output = None
if "animal_folders" not in st.session_state:
    st.session_state.animal_folders = None
if "csv_name" not in st.session_state:
    st.session_state.csv_name = None

def list_subfolders(parent_dir: str):
    return sorted(
        [os.path.join(parent_dir, d) for d in os.listdir(parent_dir)
         if os.path.isdir(os.path.join(parent_dir, d))],
        key=natural_sort_key
    )

def list_files(folder: str):
    return sorted(
        [os.path.join(folder, f) for f in os.listdir(folder)
         if os.path.isfile(os.path.join(folder, f))],
        key=natural_sort_key
    )

if st.button("Run Analysis"):
    st.write("Calculating...")

    csv_output = []
    for i, foldername in enumerate(os.listdir(unzip_folder)):
        run_path = os.path.join(unzip_folder, foldername)
        value = run_max_inst_power(run_path, mass_kg=mass_g[i] * 0.001)
        csv_output.append(value)

    # store results so they survive reruns
    st.session_state.csv_output = csv_output
    st.session_state.animal_folders = list_subfolders(unzip_folder)

    now = datetime.now()
    st.session_state.csv_name = f"peak_power_{now.strftime('%Y_%m_%d_%H_%M_%S')}.csv"
    st.session_state.analysis_done = True

if st.session_state.analysis_done:
    csv_output = st.session_state.csv_output
    animal_folders = st.session_state.animal_folders or []

    # main plot
    st.write("Graphing...")
    fig, ax = plt.subplots(figsize=(11, 8))
    for idx, result in enumerate(csv_output):
        x_coord = np.arange(len(result))
        ax.plot(x_coord, np.array(result), label=f"Folder {idx + 1}")
    ax.set_xlabel("Contraction Index")
    ax.set_ylabel("Normalized Power (W/kg)")
    ax.set_title("Peak Power")
    ax.legend()
    ax.grid(True)
    st.pyplot(fig)

    # download
    df = pd.DataFrame(csv_output)
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="Download CSV",
        data=csv_bytes,
        file_name=st.session_state.csv_name,
        mime="text/csv",
    )

    st.subheader("Validate Results")

    if len(animal_folders) == 0:
        st.warning("No animal subfolders found.")
    else:
        # pick which animal folder to inspect
        animal_choice = st.selectbox(
            "Select animal folder",
            options=animal_folders,
            format_func=lambda p: os.path.basename(p),
        )

        open_animal_folder = [f for f in list_files(animal_choice)
                              if f.lower().endswith(".ddf")]

        if len(open_animal_folder) == 0:
            st.warning("No files found in the selected animal folder.")
        else:
            mode = st.radio(
                "Validation mode",
                ["Random Sample", "Specific Inspection"],
                horizontal=True,
            )

            if mode == "Random Sample":
                pct = st.slider("Sample fraction", 0.0, 1.0, 0.05, 0.01)
                n_samples = max(1, int(pct * len(open_animal_folder)))

                if st.button("Run Random Sample"):
                    random_indices = sorted(random.sample(
                        range(len(open_animal_folder)),
                        k=min(n_samples, len(open_animal_folder)),
                    ))
                    for i in random_indices:
                        val_path = open_animal_folder[i]
                        fig_l, fig_f, fig_p = val_max_inst_power(val_path, i)
                        st.pyplot(fig_l)
                        st.pyplot(fig_f)
                        st.pyplot(fig_p)

            else:  # Specific Inspection
                start = st.number_input(
                    "Start index (inclusive)",
                    min_value=0,
                    max_value=len(open_animal_folder) - 1,
                    value=0,
                    step=1,
                )
                end = st.number_input(
                    "End index (inclusive)",
                    min_value=0,
                    max_value=len(open_animal_folder) - 1,
                    value=len(open_animal_folder) - 1,
                    step=1,
                )

                if start > end:
                    st.error("Start index must be less than or equal to end index.")
                else:
                    if st.button("Run Specific Inspection"):
                        for i in range(int(start), int(end) + 1):
                            val_path = open_animal_folder[i]
                            fig_l, fig_f, fig_p = val_max_inst_power(val_path, i)
                            st.pyplot(fig_l)
                            st.pyplot(fig_f)
                            st.pyplot(fig_p)
