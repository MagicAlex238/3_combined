# combined
combined repo for physicochemical and bacteriological data

| **1_visual_protein.ipynb** | A notebook to visualise the computational results of the picrust2 prediction|     
| **1_dashboard.py** | A dashboard file to visualise the computational results of the picrust2 prediction|    
| **2_cluster_prepare.ipynb** | A notebook to cluster the data and combine it|   
| **3_Nobel_MIC.ipynb** | A notebook to discuss nobel Microorganism associated with corrosion pro and contra|   
| **4_Neuronetwork.ipynb** | A notebook to do the prediction of corrosion |  

#=================================================================================
# 1. Install Python 3.11 and create a virtual environment

This project works on Linux and Windows. Choose the instructions that match your machine.

## Windows (recommended for quick setup)
- Official installer: download the Windows installer for Python 3.11 from https://www.python.org/downloads/windows/ and run it. Select **Add Python to PATH**.
- Microsoft Store: search for **Python 3.11** in the Microsoft Store and install.
- winget (may be blocked by Group Policy):
	```powershell
	winget install --id Python.Python.3.11 -e
	```
	If `winget` is blocked, use the official installer or Miniconda.

## Miniconda / Conda (alternative, useful on restricted systems)
- Download Miniconda from https://docs.conda.io/en/latest/miniconda.html and install for your user.
- Create an environment:
	```powershell
	conda create -n combi python=3.11 -y
	conda activate combi
	```

## Linux (Ubuntu / Debian) or WSL
Inside a Linux shell (or WSL) use apt:
```bash

# 1. Install Python 3.11 and venv
sudo apt update
sudo apt install python3.11 python3.11-venv

# 2. Install pip if needed
sudo apt install python3-pip

# 3. Create and activate virtual environment
python3.11 -m venv .venv_combi
source .venv_combi/bin/activate
pip install --upgrade pip setuptools wheel
pip install pip-tools pipdeptree
# 4. Upgrade pip and install requirements
pip install --upgrade pip
pip install -r requirements.txt
# NOTE: Google Chrome must be installed for full functionality.
wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb && sudo apt install ./google-chrome-stable_current_amd64.deb

# 5. Register the kernel with Jupyter
python -m ipykernel install --user --name .venv_combi --display-name ".venv_combi"

```

## Create & activate virtual environment (PowerShell)
From the repository root (`C:\Users\hhgwtemp39\Documents\GitHub\3_combined`):
```powershell
# Using the installed Python
python -m venv .venv_combi
.\.venv_combi\Scripts\Activate.ps1
```
If PowerShell blocks script execution, allow local scripts for the current user:
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
# then activate
.\.venv_combi\Scripts\Activate.ps1
```
Alternative activation without changing ExecutionPolicy:
```powershell
cmd /c .\.venv_combi\\Scripts\\activate.bat
```

## Install dependencies and register the Jupyter kernel
With the venv active:
```powershell
pip install --upgrade pip
pip install -r requirements.txt
python -m ipykernel install --user --name python3.11_combi --display-name "python3.11_combi"
```

## VS Code tips
- Install the `ms-python.python` and `ms-toolsai.jupyter` extensions.
- `Ctrl+Shift+P` → `Python: Select Interpreter` → choose `./.venv_combi/Scripts/python.exe`.
- For notebooks, select the kernel named `python3.11_combi`.

