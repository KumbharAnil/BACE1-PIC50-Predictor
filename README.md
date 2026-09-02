
# 🧬 BACE1 pIC50 Predictor


## 📌 Overview

**BACE1 pIC50 Predictor** is a QSAR-based computational tool designed to predict the inhibitory activity of molecules against **BACE1 (Beta-secretase 1)**.

The tool accepts molecular **SMILES strings** as input, automatically calculates the required molecular descriptors using **PaDEL-Descriptor**, applies a pre-developed QSAR equation, and generates predicted:

- 🧪 **pIC50**
- 🔬 **IC50 (nM)**
- 📊 Molecular descriptor values
- ✅ Prediction status

The tool supports **batch processing**, automatic error handling, and CSV-based input/output.

---

## 🚀 Key Features

| Feature | Description |
|---|---|
| 🧪 BACE1 Prediction | Predict BACE1 inhibitory activity |
| 📂 CSV Input | Read molecules directly from CSV files |
| 🧬 SMILES Processing | Calculate molecular descriptors from SMILES |
| 📊 QSAR Prediction | Predict pIC50 using a linear QSAR equation |
| 🔬 IC50 Conversion | Convert predicted pIC50 to IC50 in nM |
| ⚡ Batch Processing | Process molecules in configurable chunks |
| 🛡️ Error Handling | Isolate problematic molecules automatically |
| ❌ Failure Tracking | Record molecules that cannot be processed |
| 📁 CSV Output | Automatically generate prediction results |

---

## 🔬 How It Works

```text
                    Input CSV
                        │
                        ▼
                Read Molecule IDs
                        │
                        ▼
                  Read SMILES
                        │
                        ▼
             PaDEL-Descriptor
                        │
                        ▼
        ┌───────────────────────────┐
        │ Molecular Descriptors     │
        │                           │
        │ • SIC1                    │
        │ • SpMin4_Bhm              │
        │ • SC-5                    │
        └───────────────────────────┘
                        │
                        ▼
                  QSAR Equation
                        │
                        ▼
                     pIC50
                        │
                        ▼
                  IC50 (nM)
                        │
                        ▼
                Prediction CSV
````

---

## 🧮 QSAR Model

The predictor uses three molecular descriptors:

* **SIC1**
* **SpMin4_Bhm**
* **SC-5**

The QSAR equation is:

$$
pIC50 =
-15.13463
+ 13.84656(SIC1)
+ 7.33758(SpMin4\_Bhm)
+ 1.90156(SC-5)
$$

The corresponding IC50 value is calculated as:

$$
IC50(M) = 10^{-pIC50}
$$

and:

$$
IC50(nM) = IC50(M) \times 10^9
$$

### Python implementation

```python
def compute_pic50(sic1, spmin4_bhm, sc5):

    pIC50 = (
        -15.13463
        + (13.84656 * sic1)
        + (7.33758 * spmin4_bhm)
        + (1.90156 * sc5)
    )

    IC50_nM = (10 ** (-pIC50)) * 1000000000

    return pIC50, IC50_nM
```

---

## 🧬 Molecular Descriptors

| Descriptor     | Description                                    |
| -------------- | ---------------------------------------------- |
| **SIC1**       | Information-content based molecular descriptor |
| **SpMin4_Bhm** | Burden/graph-based molecular descriptor        |
| **SC-5**       | Kier/graph-based molecular descriptor          |

These descriptors are automatically calculated from the molecular SMILES using **PaDEL-Descriptor** through `padelpy`.

---

## 📁 Project Structure

```text
BACE1-PIC50-Predictor/
│
├── app.py
│
├── assets/
│   ├── Inactive Molecule.csv
│   └── BACE1_Active_predictions.csv
│
├── README.md
│
└── .gitignore
```

---

# 🛠️ Installation

## 1️⃣ Requirements

Before running the project, make sure you have:

* Python **3.9+**
* Java Runtime Environment (JRE) or JDK
* pip
* Git *(optional)*

Python dependency:

```text
padelpy
```

---

## 2️⃣ Clone the Repository

```bash
git clone https://github.com/KumbharAnil/BACE1-PIC50-Predictor.git
```

Move into the project directory:

```bash
cd BACE1-PIC50-Predictor
```

---

## 3️⃣ Create a Virtual Environment

```bash
python -m venv .venv
```

### Windows PowerShell

```powershell
.venv\Scripts\Activate.ps1
```

You should see:

```text
(.venv)
```

at the beginning of your terminal prompt.

---

## 4️⃣ Install Dependencies

```bash
pip install padelpy
```

Verify the installation:

```bash
pip show padelpy
```

---

# ☕ Java Requirement

PaDEL-Descriptor requires Java.

Check whether Java is installed:

```bash
java -version
```

Example:

```text
java version "..."
```

If Java is not recognized, install a JRE/JDK and make sure Java is available in your system `PATH`.

---

# ▶️ Running the Predictor

Activate your virtual environment:

```powershell
.venv\Scripts\Activate.ps1
```

Then run:

```bash
python app.py
```

The program will display something similar to:

```text
==========================================
     BACE1 QSAR BATCH PREDICTION TOOL
==========================================

Loaded ... molecules from ./assets/Inactive Molecule.csv

Processing 1-50 of ...
Processing 51-100 of ...
Processing 101-150 of ...

==========================================
Done. ... succeeded, ... failed.

Results written to:
./assets/BACE1_Active_predictions.csv
==========================================
```

---

# 📥 Input Dataset

The default input file is:

```text
assets/Inactive Molecule.csv
```

The CSV should contain the required columns:

```text
pub chem SID
Smile
```

### Example

```csv
pub chem SID,Smile
Molecule_001,CCO
Molecule_002,CC(C)O
Molecule_003,c1ccccc1
```

> ⚠️ **Important:** The current configuration uses the exact column name `Smile ` with a trailing space.

If your CSV contains:

```text
Smile
```

instead of:

```text
Smile 
```

change:

```python
SMILES_COLUMN = "Smile "
```

to:

```python
SMILES_COLUMN = "Smile"
```

Column names must match the configuration exactly.

---

# 📤 Output

The prediction file is generated at:

```text
assets/BACE1_Active_predictions.csv
```

The output contains:

| Column       | Description                      |
| ------------ | -------------------------------- |
| `ID`         | Molecule identifier              |
| `SMILES`     | Molecular SMILES                 |
| `SIC1`       | Calculated SIC1 descriptor       |
| `SpMin4_Bhm` | Calculated SpMin4_Bhm descriptor |
| `SC-5`       | Calculated SC-5 descriptor       |
| `pIC50`      | Predicted pIC50                  |
| `IC50_nM`    | Predicted IC50 in nanomolar      |
| `status`     | Prediction status                |

### Example

```csv
ID,SMILES,SIC1,SpMin4_Bhm,SC-5,pIC50,IC50_nM,status
Molecule_001,CCO,...,...,...,...,...,OK
```

---

# ⚡ Batch Processing

The predictor processes molecules in configurable chunks.

Default:

```python
CHUNK_SIZE = 50
```

This means 50 molecules are processed at a time.

You can increase or decrease the value:

```python
CHUNK_SIZE = 100
```

Larger datasets can therefore be processed without attempting to calculate all descriptors in a single operation.

---

# 🛡️ Error Handling

The predictor uses a **two-stage processing strategy**.

### Stage 1 — Batch Processing

The program first attempts to process an entire chunk:

```python
desc_list = from_smiles(
    smiles_list,
    fingerprints=False
)
```

### Stage 2 — Individual Processing

If the complete batch fails, the program automatically attempts to process molecules individually.

This helps isolate problematic SMILES strings.

Instead of stopping the complete prediction process, failed molecules are recorded with:

```text
status = FAILED: ...
```

Therefore, a problematic molecule does not stop the entire dataset from being processed.

---

# 🔧 Configuration

The main configuration variables are located near the beginning of `app.py`.

```python
INPUT_CSV = "./assets/Inactive Molecule.csv"

ID_COLUMN = "pub chem SID"

SMILES_COLUMN = "Smile "

OUTPUT_CSV = "./assets/BACE1_Active_predictions.csv"

CHUNK_SIZE = 50
```

You can modify these values according to your dataset.

### Example

```python
INPUT_CSV = "./assets/my_molecules.csv"

OUTPUT_CSV = "./assets/my_predictions.csv"
```

---

# 🧪 Complete Example

### Windows PowerShell

```powershell
cd "C:\Users\Anil\OneDrive\Desktop\ProjectsDraft\Project_Data"
```

Activate the environment:

```powershell
.venv\Scripts\Activate.ps1
```

Install the dependency:

```powershell
pip install padelpy
```

Run the predictor:

```powershell
python app.py
```

View the generated output:

```powershell
Get-Content ".\assets\BACE1_Active_predictions.csv"
```

---

# 📊 Prediction Workflow

```text
SMILES
  │
  ▼
PaDEL-Descriptor
  │
  ├── SIC1
  ├── SpMin4_Bhm
  └── SC-5
        │
        ▼
   QSAR Equation
        │
        ▼
      pIC50
        │
        ▼
     IC50 (nM)
```

---

# 📌 Scientific Interpretation

The generated pIC50 and IC50 values are **computational QSAR predictions** and should not be interpreted as experimental measurements.

Prediction reliability depends on factors including:

* Training dataset quality
* Descriptor calculation
* Model validation
* Applicability domain
* Molecular similarity to the training compounds
* Chemical space represented by the model

Predictions outside the model's applicability domain should therefore be interpreted cautiously.

---

# ⚠️ Important Notes

### SMILES Column

Make sure the configured SMILES column exactly matches your CSV.

```python
SMILES_COLUMN = "Smile "
```

and:

```python
SMILES_COLUMN = "Smile"
```

are treated as different column names.

---

### Virtual Environment

Do **not** upload `.venv/` to GitHub.

Add it to `.gitignore`:

```text
.venv/
__pycache__/
*.pyc
```

---

### Large Datasets

Avoid committing extremely large datasets or generated prediction files directly to GitHub.

For large datasets, consider:

* Git LFS
* External data storage
* Releasing datasets separately

---

# 🔮 Future Improvements

The project can be extended with:

* 🖥️ Streamlit web interface
* 📊 Interactive prediction dashboard
* 📈 Activity visualization
* 🧬 Molecular structure visualization
* 📁 Drag-and-drop CSV upload
* 📊 Prediction statistics
* 🔍 Molecule filtering
* 📥 Downloadable prediction reports
* 🧪 Additional QSAR models
* 🤖 Machine-learning model comparison
* 📈 Model validation metrics
* 🧠 Applicability-domain visualization

---

# 🤝 Contributing

Contributions, suggestions, and improvements are welcome.

To contribute:

```bash
git fork
```

Create a new branch:

```bash
git checkout -b feature/new-feature
```

Commit your changes:

```bash
git add .
git commit -m "Add new feature"
```

Push the branch:

```bash
git push origin feature/new-feature
```

Then open a Pull Request.

---

# 📚 Technologies Used

| Technology          | Purpose                          |
| ------------------- | -------------------------------- |
| 🐍 Python           | Core programming                 |
| ☕ Java              | Runtime for PaDEL                |
| 🧬 PaDEL-Descriptor | Molecular descriptor calculation |
| 📦 padelpy          | Python interface for PaDEL       |
| 📄 CSV              | Dataset and prediction storage   |
| 🧪 QSAR             | Molecular activity prediction    |

---

# 👨‍🔬 Project

## BACE1 pIC50 Predictor

A computational drug-discovery project implementing a **QSAR-based batch prediction workflow for BACE1 inhibitory activity**.

The project demonstrates the integration of:

```text
SMILES
   ↓
Molecular Descriptors
   ↓
QSAR Model
   ↓
pIC50 Prediction
   ↓
IC50 Prediction
```

---

# ⭐ Support

If you find this project useful, consider giving the repository a ⭐ on GitHub.

**Repository:**

[KumbharAnil/BACE1-PIC50-Predictor](https://github.com/KumbharAnil/BACE1-PIC50-Predictor)

---

## 📜 Disclaimer

This software is intended for **research and educational purposes**.

The predictions generated by this tool are computational estimates from a QSAR model and should not be considered experimental BACE1 activity measurements or used as a substitute for experimental validation.

---

<p align="center">

🧬 <b>BACE1 pIC50 Predictor</b>

<br>

<sub>QSAR • Molecular Descriptors • Computational Drug Discovery</sub>

</p>
```

