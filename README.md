# 📂 **Wafer Data Analysis & Automated Reporting System**

![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-Data%20Computation-013243?logo=numpy&logoColor=white)
![Pandas](https://img.shields.io/badge/Pandas-Data%20Analysis-150458?logo=pandas&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-Advanced%20Math-8CAAE6?logo=scipy&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-Visualization-11557c?logo=python&logoColor=white)
![Jupyter/iPyWidgets](https://img.shields.io/badge/iPyWidgets-Interactive-F37626?logo=jupyter&logoColor=white)

---
# 📝Contents
1. Introduction
2. Quick Run
3. Project information
4. Directory Structure
5. Data Analysis Pipeline
6. Jupyter Notebook
---
# 1. Introduction
-PE2: team2 project

- We have automated a 3-stage pipeline, encompassing everything from raw data processing to visualization and analysis.

### 1️⃣ Data Extraction & Preprocessing
* **Execution**: The system automatically searches for and loads raw `XML` data files, followed by a preprocessing stage to cleanse missing values (nulls) and data anomalies.
* **Deliverables**: The cleansed data is stored in `CSV` and `Xlsx` formats for seamless subsequent analysis and archiving, and initial `Png` plots are automatically generated for primary data validation.

### 2️⃣ Statistical Analysis
* **Execution**: Conducts large-scale statistical analysis based on the preprocessed structured data.
* **Deliverables**: Automatically renders and reports **Box plots (Box graphs)** to grasp data distribution and outliers at a glance, along with **Wafer Maps** that visualize defect locations and spatial uniformity across the wafer.

### 3️⃣ Die-Level Interactive Analysis
* **Execution**: Moves beyond static reporting to provide a customized analytical environment where clients can track specific die-level sensor data in-depth.
* **Deliverables**: Delivers a **Jupyter Notebook-based interactive dashboard UI** that enables real-time user-driven data filtering.

```markdown
```text
[Raw Data (XML)] 
        │
        ▼ 
 ┌──────────────────────────────────────────┐
 │  Data Extraction & Preprocessing         │ ───▶ CSV, Xlsx, Png
 └──────────────────────────────────────────┘
        │
        ▼ 
 ┌──────────────────────────────────────────┐
 │  Statistical Analysis                    │ ───▶ Wafer Map, Box graph
 └──────────────────────────────────────────┘
        │
        ▼ 
 ┌──────────────────────────────────────────┐
 │  Die-Level Interactive Analysis          │ ───▶ Jupyter Notebook
 └──────────────────────────────────────────┘
```
---

## 👥 Contributors

|     Name     |               E-mail               | 
|:------------:|:----------------------------------:|
| Lee Hangyol  |       LeeHanGyol@hanyang.kr        |
| Jung Jae min |      jjmin0118@hanyang.ac.kr       |
| Lee HyoSeong |      lhs012100@hanyang.ac.kr       |
|  Kim HanSeo  |     lincoln1111@hanyang.ac.kr      |
---
# 2. Quick Run

### Installation

```bash
pip install numpy matplotlib scipy pandas ipywidgets
```

```bash
pip install -r requirements.txt
```

### ▶️ How to Run

1. Place the raw XML data folder you want to analyze into the `dat/` directory
2. Execute the pipeline

```bash
python run.py
```

3. Results are automatically saved in the `res/` directory
---
# 3. Project information
### 📖 About This Project

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

### 6) Jupyter Notebook
- **`wafer_analysis`** — The purpose is to visually extract the necessary data at the die level.

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
│   ├── 📄 ref_poly.py                   
│   ├── 📄 flatting.py              
│   ├── 📄 zoom.py                  
│   ├── 📄 Fitting.py               
│   ├── 📄 Phase shift - V.py       
│   ├── 📄 VpiL.py                  
│   ├── 📄 ER_Analysis.py           
│   ├── 📄 IL_Analysis.py           
│   ├── 📄 VpiL_Analysis.py         
│   ├── 📄 combine_plot.py
│   ├── 📄 wafer_analyzer.py          
│   └── 📄 export_summary.py        
│
└── 📁 res/                         # Output directory automatically generated at runtime
    ├── 📁 csv/                     # Collection of summary and consolidated data CSVs
    │   ├── 📄 Analysis.csv             
    │   ├── 📄 Total_Process_result.csv 
    │   └── 📄 {Wafer_ID}_Process_result.csv 
    │
    ├── 📁 xlsx/                    # Collection of consolidated Excel reports with hyperlinks
    │   ├── 📄 Analysis.xlsx            
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
# 5. ⚙️ Data Analysis Pipeline

When `run.py` is executed, a total of 9 core modules operate 
sequentially to process the data.

### 1. Data Extraction & Visualization Preparation

- **`ref_poly.py`** (Remove Ref)
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
<img width="989" height="590" alt="D07_C0_R0_LMZC_Phase" src="https://github.com/user-attachments/assets/718d1b59-97df-4144-933f-b793e9defc17" />

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
### CSV
<img width="1122" height="307" alt="image" src="https://github.com/user-attachments/assets/c0a1a98c-fcf3-42f5-8ae2-bfa8dcac8c6a" />

### Xlsx
<img width="1807" height="290" alt="image" src="https://github.com/user-attachments/assets/fef3c999-d39f-46d5-b748-cfd04a0c5e71" />

# 6. Jupyter Notebook

- **`wafer_analyzer.py`** (Die-Level Data Analysis)
