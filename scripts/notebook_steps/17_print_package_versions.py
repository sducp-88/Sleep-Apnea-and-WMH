# Extracted from WMH_Analysis.ipynb, code cell 16.
# Notebook heading: #Print Versions of all packages used in the analysis for reproducibility purposes
# Run this file from the repository root unless a local CONFIG section is edited.

#Print Versions of all packages used in the analysis for reproducibility purposes
import importlib

packages = [
    "sys", "pandas", "numpy", "scipy", "statsmodels",
    "matplotlib", "seaborn", "sklearn",
    "lifelines", "gseapy", "networkx", "graphviz",
    "openpyxl", "docx", "tqdm"
]

for pkg in packages:
    try:
        if pkg == "sys":
            import sys
            print(f"Python: {sys.version}")
        else:
            module = importlib.import_module(pkg)
            print(f"{pkg}: {module.__version__}")
    except ImportError:
        print(f"{pkg}: not installed")
    except AttributeError:
        print(f"{pkg}: version attribute not found")
