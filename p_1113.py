"""
power.py

Description: Process all animals in a folder, compute normalized maximum instantaneous 
power for each contraction, and return a figure plus the raw results.
"""

__author__ = "Chaoran Yang"
__version__ = "2.0"
__email__ = "cy197@duke.edu"
__date__ = "2025-11-13"

import os
import re
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import streamlit as st
import zipfile
import datetime

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
    
    with open(file_path, "r") as f:
        lines = f.readlines()

    data_marker_idx = None
    sample_freq_hz = None

    # Parse header-level info, including calibration block
    channel_names: List[str] = []
    units: List[str] = []
    scales: List[float] = []
    offsets: List[float] = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Locate data marker
        if line.startswith("Test Data in Volts"):
            data_marker_idx = i
            break

        # Sample frequency
        if line.startswith("Sample Frequency"):
            # e.g. "Sample Frequency (Hz): 1000"
            try:
                sample_freq_hz = float(line.split(":")[1].strip())
            except (IndexError, ValueError):
                raise ValueError(f"Could not parse sample frequency in {file_path}")

        # Calibration block
        if line.startswith("Channel"):
            # Format:
            # Channel <tab> AI0 AI1 ... AO0 AO1
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
    }

def detect_contraction_window(
    force_mN: np.ndarray,
    fs_hz: float,
    threshold_std: float = 5.0,
    min_duration_ms: float = 20.0
) -> tuple[int, int, slice]:
    """
    Automatically detect contraction start/end and define a baseline segment.

    Strategy:
      1. Estimate an initial baseline from the first ~20% of the trace.
      2. Find where |force - baseline_mean| exceeds baseline_std * threshold_std.
      3. Define the contraction as the contiguous region where this condition holds.
      4. Define the baseline segment as the middle 60% of the pre-contraction region.

    Parameters
    ----------
    force_mN : np.ndarray
        Force trace in mN.
    fs_hz : float
        Sample frequency in Hz.
    threshold_std : float
        Multiplier on baseline standard deviation to define onset.
    min_duration_ms : float
        Minimal contraction duration to avoid false positives.

    Returns
    -------
    start_idx : int
        Index where contraction begins.
    end_idx : int
        Index where contraction ends.
    baseline_slice : slice
        Slice object selecting the baseline region (pre-contraction).
    """
    n = len(force_mN)
    if n == 0:
        raise ValueError("Empty force trace.")

    # Take the first 20% (at least 100 samples) as an initial baseline region
    initial_baseline_end = max(int(0.2 * n), 100)
    initial_baseline_end = min(initial_baseline_end, n)
    baseline_region = force_mN[:initial_baseline_end]

    baseline_mean = float(np.mean(baseline_region))
    baseline_std = float(np.std(baseline_region))

    if baseline_std == 0:
        # Flat signal; cannot detect contraction meaningfully
        raise ValueError("Baseline standard deviation is zero; cannot detect contraction.")

    deviation = np.abs(force_mN - baseline_mean)
    threshold = threshold_std * baseline_std
    active_mask = deviation > threshold

    min_samples = int((min_duration_ms / 1000.0) * fs_hz)
    if min_samples < 1:
        min_samples = 1

    # Find the first index where we have a contiguous run of 'min_samples' active points
    start_idx = None
    i = 0
    while i < n - min_samples:
        if np.all(active_mask[i:i + min_samples]):
            start_idx = i
            break
        i += 1

    if start_idx is None:
        raise ValueError("No contraction detected in force trace.")

    # End index: last point where deviation is above threshold after the start
    active_indices = np.where(active_mask)[0]
    end_idx = active_indices[active_indices >= start_idx].max()

    # Define a baseline slice: middle 60% of the pre-contraction region
    if start_idx < 5:
        raise ValueError("Contraction starts too early; not enough baseline region.")

    pre_region_len = start_idx
    baseline_start = int(0.2 * pre_region_len)
    baseline_end = int(0.8 * pre_region_len)
    if baseline_end <= baseline_start:
        baseline_start = 0
        baseline_end = pre_region_len

    baseline_slice = slice(baseline_start, baseline_end)

    return start_idx, end_idx, baseline_slice

def max_instantaneous_power_from_file(
    file_path: str,
    threshold_std: float = 5.0,
    min_duration_ms: float = 20.0,
    flip_length_sign: bool = True
) -> float:
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

    # Optionally flip length so shortening is positive
    if flip_length_sign:
        length_mm = -length_mm

    fs_hz = 1.0 / np.mean(np.diff(time))

    start_idx, end_idx, baseline_slice = detect_contraction_window(
        force_mN, fs_hz, threshold_std=threshold_std, min_duration_ms=min_duration_ms
    )

    # Baseline correction using the pre-contraction middle segment
    length_baseline = float(np.mean(length_mm[baseline_slice]))
    force_baseline = float(np.mean(force_mN[baseline_slice]))

    length_seg = length_mm[start_idx:end_idx] - length_baseline
    force_seg = force_mN[start_idx:end_idx] - force_baseline
    time_seg = time[start_idx:end_idx]

    # Velocity (mm/s)
    velocity_mm_s = np.gradient(length_seg, time_seg)

    # Power: force (mN) * velocity (mm/s) -> W via 1e-6 factor
    inst_power_W = 1e-6 * force_seg * velocity_mm_s
    max_power_W = float(np.max(inst_power_W))

    return max_power_W

def run_max_inst_power(
    folder_path: str,
    mass_kg: float,
    threshold_std: float = 5.0,
    min_duration_ms: float = 20.0
) -> List[List[float]]:
    """
    Process one animal folder under `folder_path` and compute normalized
    peak instantaneous power (W/kg) for each contraction.

    Parameters
    ----------
    folder_path : str
        Root directory of the one animal folder.
    mass_kg : float
        Animal mass in kilograms used to normalize peak power.
    threshold_std : float
        Threshold in units of baseline standard deviation for contraction detection.
    min_duration_ms : float
        Minimum contraction duration in milliseconds to avoid false positives.

    Returns
    -------
    outputs : list[list[float]]
        Normalized peak power values (W/kg) for each contraction and animal.
    """
    animal_results: List[float] = []

    for f in sorted(os.listdir(folder_path), key=natural_sort_key):
        data_file = os.path.join(folder_path, f)

        # Skip directories and hidden files
        if os.path.isdir(data_file) or f.startswith("."):
            print(f"  Skipped {f} in {folder_path}")
            continue

        print(f"Processing {f} ...")

        peak_power_W = max_instantaneous_power_from_file(
            data_file,
            threshold_std=threshold_std,
            min_duration_ms=min_duration_ms,
        )
        normalized_power = peak_power_W / mass_kg  # W/kg
        animal_results.append(normalized_power)

    return animal_results

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

    mass_kg = []

    # You can now loop through the contents of the folder
    for i, filename in os.listdir(unzip_folder):
        mass_kg[i] = st.number_input("Mass (kg):", min_value=0.0, value=1.0, step=0.001) 
        # check if mass is actually in kg

if st.button("Run Analysis"):
    st.write("Calculating...")
    csv_output = []
    for i, filename in enumerate(os.listdir(unzip_folder)): 
        csv_output[i] = run_max_inst_power(filename, mass_kg=mass_kg[i])
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

    now = datetime.now()
    timestamp = now.strftime("%Y_%m_%d_%H_%M_%S")
    csv_name=f"peak_power_{timestamp}.csv"

    df = pd.DataFrame(csv_output)
    csv = df.to_csv(index=False).encode('utf-8')    
    st.download_button(label="Download CSV",
            data=csv,
            file_name=csv_name,
            mime='text/csv')

# add the option to save mass data
# add the option to upload mass data (from previous generation)