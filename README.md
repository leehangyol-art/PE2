# 📂 Wafer Data Analysis & Integrated Reporting System

![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python&logoColor=white)
![Pandas](https://img.shields.io/badge/Pandas-Data%20Analysis-150458?logo=pandas&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-Visualization-black?logo=python)
---
# 📝Contents
1. Introduction
2. Project information
3. Directory Structure
4. Install and Run
---
> **This is an automated system that parses XML-based wafer measurement data and performs precise analysis of IL, ER, VpiL, and Phase Shift to generate visualizations and integrated Excel reports.**

This project is an end-to-end pipeline that automatically analyzes wafer-level test data from optical devices such as Mach-Zehnder Modulators and outputs the results as a dashboard-style report for at-a-glance visualization.

---
# 1. Introduction
-PE2: WOP project

Our goal is to build a system that automatically analyzes optical characteristics using Python and XML data.
## 👥 Contributors

|     Name      |       E-mail       | 
|:-------------:|:------------------:|
|  Lee Hangyol  | 0000@hanyang.ac.kr |
| Jeong Jae-min | 0000@hanyang.ac.kr |
| Lee HyoSeong  | 0000@hanyang.ac.kr |
|  Kim HanSeo   | 0000@hanyang.ac.kr |
---
# 2. Project information
### 📖 About This Project

Our goal is to build a system that automatically analyzes
optical characteristics using Python and XML data.

The system receives raw wafer data from the customer along
with specific analysis requests, including:

- Target Wafer
- Die Row & Column
- Analysis Options

Based on these inputs, the pipeline acts as a **black box** —
automatically processing the requested data and delivering
the final results without manual intervention.


### ▶️ How It Works

1. Extract files named `HY202103` from the customer-provided dataset
2. Process raw XML data through a series of analysis modules
3. Export the final results as `.csv` and `.xlsx` reports
---
# 3. ⚙️ Key Features


### 1) Data Extraction & Visualization Preparation

- **`ref_poly.py`** — Removes noise from REF signals and establishes a stable baseline for analysis
- **`data_parser.py`** — Parses target band (LMZC / LMZO) spectrum data from raw XML files
- **`plot.py`** — Generates baseline wavelength-transmission spectrum plots from parsed raw data


### 2) Signal Correction & Target Region Extraction

- **`flatting.py`** — Flattens the signal baseline by correcting offset errors between Reference and MZM devices
- **`zoom.py`** — Zooms into key analysis wavelength ranges per band (LMZC: 1550 nm / LMZO: 1310 nm)
- **`Fitting.py`** — Removes high-frequency noise (ripple) and applies polynomial fitting for smooth data refinement


### 3) Device Performance Metric Calculation

- **`Phase shift - V.py`** — Tracks wavelength shifts according to applied bias voltage based on dip positions in fitted graphs
- **`VpiL.py`** — Calculates VπL (electro-optic modulation efficiency) by converting phase shifts into half-wave voltage (Vπ) and multiplying by device length (L)


### 4) Wafer Map & Box Plot Auto-Generation

- **`ER_Analysis.py`** — Generates Wafer Map and Box Graph for Extinction Ratio (ER)
- **`IL_Analysis.py`** — Generates Wafer Map and Box Graph for Insertion Loss (IL)
- **`VpiL_Analysis.py`** — Generates Wafer Map and Box Graph for VpiL


### 5) Visualization Merging & Report Auto-Generation

- **`combine_plot.py`** — Merges individual analysis graphs into a single summary dashboard image per wafer and measurement date
- **`export_summary.py`** — Exports final IL / ER / VπL data as `.csv` and `.xlsx` files, with hyperlinks in Excel mapped to merged summary images for intuitive one-click data verification
------
# 4. 📁 Directory Structure

The project is structured to separate raw inputs, 
processed outputs, and source logic clearly

```text
📁 pycharm-project-root/
│
├── 📄 run.py                       # Integrated automation entry point (One-click Execution)
│
├── 📁 dat/
│   └── 📄 data                     # Raw wafer XML compressed data (ex: HY202103.zip)
│
├── 📁 src/                         # Core analysis module directory
│   ├── 📄 data_parser.py           
│   ├── 📄 plot.py                  
│   ├── 📄 flatting.py              
│   ├── 📄 zoom.py                  
│   ├── 📄 Fitting.py               
│   ├── 📄 Phase shift - V.py       
│   ├── 📄 VpiL.py                  
│   ├── 📄 ER_Analysis.py           
│   ├── 📄 IL_Analysis.py           
│   ├── 📄 VpiL_Analysis.py         
│   ├── 📄 combine_plot.py          
│   └── 📄 export_summary.py        
│
└── 📁 res/                         # Output directory automatically generated at runtime
    ├── 📁 csv/                     # Collection of summary and consolidated data CSVs
    │   ├── 📄 Analysis.csv             
    │   ├── 📄 Total_Process_result.csv 
    │   └── 📄 {Wafer_ID}_Process_result.csv 
    │
    ├── 📁 xlsx/                    # Collection of consolidated Excel reports with hyperlinks
    │   ├── 📄 Analysis.xlsm            
    │   ├── 📄 Total_Process_result.xlsx 
    │   └── 📄 {Wafer_ID}_Process_result.xlsx 
    │
    └── 📁 png/                     # Visualization image storage
        ├── 📁 WaferMap/                # Per-wafer ER, IL, and VpiL heatmaps
        ├── 📁 BoxPlot/                 # Per-wafer Center vs Edge box plots
        └── 📁 {Wafer_ID}/              # Individual die data
            └── 📁 {Date_YYYYMMDD}/     # Per-measurement-date folders 
                └── 📄 HY202103_{Wafer}_({C},{R})_LION1_DCM_{Band}.png  # Merged summary images
```
---
# 5. Install and Run

### Installation

```bash
pip install -r requirements.txt
```

### How to Run

1. Place the raw XML data folder you want to analyze into the `dat/` directory
2. Execute the pipeline

```bash
python run.py
```

3. Select analysis options
   - Target **Wafer** and **Die** (row & column)
   - Figure options : `Show Figure` / `Save Figure` / `Save CSV`
4. Results are automatically saved in the `res/` directory
---
# 6. Input & Output Examples

여기에 그래프와 csv 사진 첨부
---

# 7. ⚙️ Data Analysis Pipeline

When `run.py` is executed, a total of 9 core modules operate 
sequentially to process the data.

### 1. Data Extraction & Visualization Preparation

- **`data_parser.py`** (Data Parsing)
  - Loads target band (LMZC, LMZO) spectrum data required 
    for analysis from the raw XML data files.
- **`plot.py`** (Raw Data Plot)
  - Generates basic wavelength-transmission spectrum plots 
    based on the parsed raw data.


<img width="989" height="590" alt="D07_C0_R0_LMZC_Raw" src="https://github.com/user-attachments/assets/96f24f7d-ff9e-44ef-946a-5471b433d626" />

### 2. Signal Correction & Target Region Extraction

- **`flattening.py`** (Signal Flattening)
  - Corrects offset errors between the Reference device and 
    MZM device to flatten the signal baseline.
<img width="989" height="590" alt="D07_C0_R0_LMZC_Flat" src="https://github.com/user-attachments/assets/27be193e-6cf5-4301-8ad5-cbfc1fffd25f" />
- **`zoom.py`** (Target Wavelength Zoom-in)
  - Zooms into key analysis wavelength ranges according to 
    each band's characteristics. 
    (LMZC: 1550 nm / LMZO: 1310 nm)
<img width="990" height="590" alt="D07_C0_R0_LMZC_Zoom" src="https://github.com/user-attachments/assets/e3782f8c-73f1-4885-9f2d-93f716d63e9d" />
- **`Fitting.py`** (Noise Filtering & Fitting)
  - Removes high-frequency noise such as ripple from the 
    measured signal and applies polynomial fitting to 
    smooth and refine the data.
<img width="987" height="590" alt="D07_C0_R0_LMZC_Fitting" src="https://github.com/user-attachments/assets/344548af-9787-42b6-9549-4a50cc5da6ec" />

### 3. Device Performance Metric Calculation

- **`Phase shift - V.py`** (Phase Shift Calculation)
  - Tracks how much the wavelength shifts according to the 
    applied bias voltage based on the dip positions in the 
    fitted graph, and calculates the phase shift.
<img width="989" height="590" alt="D07_C0_R0_LMZC_Phase" src="https://github.com/user-attachments/assets/74cbb1c8-9cdb-4ef1-98b4-fc615ac42cac" />

- **`VpiL.py`** (VπL Extraction)
  - Converts the bias-dependent phase shift into half-wave 
    voltage (Vπ) and multiplies by the device length (L) to 
    calculate the final electro-optic modulation efficiency 
    index, VπL.
<img width="989" height="590" alt="D07_C0_R0_LMZC_VpiL" src="https://github.com/user-attachments/assets/5f2617c1-50b5-404b-869d-b2d1ef6dbbb6" />
---

### 4. Wafer Map & Box Graph Auto-Generation

- **`ER_Analysis.py`** (Wafer Map & Box Graph)
  - Generates `Wafer Map` and `Box Graph` files from the 
    extracted ER data.
<img width="698" height="654" alt="Map_D07_LMZC_20190715_ER" src="https://github.com/user-attachments/assets/bcfdc059-7823-4d33-adef-03bc0929266e" />

<img width="790" height="790" alt="Box_D07_LMZC_20190715_ER_Flattened" src="https://github.com/user-attachments/assets/9b1665b5-2bec-47fb-9969-4c5249594398" />

- **`IL_Analysis.py`** (Wafer Map & Box Graph)
  - Generates `Wafer Map` and `Box Graph` files from the 
    extracted IL data.
<img width="729" height="654" alt="Map_D07_LMZC_20190715_IL" src="https://github.com/user-attachments/assets/e326e389-0a04-4acd-8ef7-09a30112c107" />

<img width="790" height="790" alt="Box_D07_LMZC_20190715_IL" src="https://github.com/user-attachments/assets/9a246003-ef34-4431-a691-c1fecb7f8c77" />

- **`VpiL_Analysis.py`** (Wafer Map & Box Graph)
  - Generates `Wafer Map` and `Box Graph` files from the 
    extracted VpiL data.

<img width="704" height="645" alt="Map_D07_LMZC_20190715_VpiL_0V" src="https://github.com/user-attachments/assets/d3f860ed-76dd-44fb-bd93-2afafd4bcfdf" />

<img width="790" height="790" alt="Box_D07_LMZC_20190715_VpiL_0V" src="https://github.com/user-attachments/assets/1481388c-0f16-419b-8881-049a174d9c78" />

### 5. Visualization Merging & Report Auto-Generation

- **`combine_plot.py`** (Dashboard Image Merging)
  - Merges multiple graphs generated throughout the analysis 
    into a single summary dashboard image, grouped by 
    wafer and measurement date.
 <img width="2187" height="654" alt="Summary_WaferMap_D07_LMZC_20190715" src="https://github.com/user-attachments/assets/4051eff6-9de0-4185-b99d-5940eab6fc8e" />
 
<img width="2370" height="790" alt="Summary_BoxPlot_D07_LMZC_20190715" src="https://github.com/user-attachments/assets/9de893ca-1882-49ac-aaf7-7b962bc48d3f" />

- **`export_summary.py`** (Consolidated Report Export)
  - Saves the extracted key metric data (IL, ER, VπL) as 
    `.csv` and `.xlsx` files. The Excel file includes 
    hyperlinks mapped to merged summary images (PNG), 
    enabling intuitive one-click data verification.
---
# ⚠️ Precautions