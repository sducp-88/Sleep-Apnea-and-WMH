# Extracted from WMH_Analysis.ipynb, code cell 0.
# Notebook heading: # Graphviz Flowchart for submission (Option 1: inline display + file export)
# Run this file from the repository root unless a local CONFIG section is edited.

# Graphviz Flowchart for submission (Option 1: inline display + file export)
# If not installed, run:  !pip install graphviz
# (Optional for LZW compression of TIFF:  !pip install pillow)

import os
from datetime import datetime
from graphviz import Digraph

# Optional: Pillow for LZW-compressed TIFF post-processing
try:
    from PIL import Image
except Exception:
    Image = None  # If Pillow is not available, TIFF is still exported (without LZW)

# ------------------------------------------------------
# 1. Output directory
# ------------------------------------------------------
outdir = "Flow_Chart"
os.makedirs(outdir, exist_ok=True)

# Auto-generate filenames with timestamp
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
basename = "flowchart_neurology"

# ------------------------------------------------------
# 2. Initialize graph
# ------------------------------------------------------
# Note: 'dpi' controls raster resolution for PNG/TIFF outputs.
dot = Digraph("Study_Flowchart", format="png")
dot.attr(rankdir="TB", size="8,10", dpi="600")   # TB = top-to-bottom layout; 600 dpi for high-res raster

# Global font settings (sans-serif fonts, e.g., Arial)
dot.attr("graph", fontname="Arial")
dot.attr("node", shape="rectangle", style="rounded",
         fontsize="10", fontname="Arial", color="black", fillcolor="white")
dot.attr("edge", fontname="Arial", fontsize="9")


# ------------------------------------------------------
# 3. Define nodes
# ------------------------------------------------------
dot.node("A", "UK Biobank participants\n ≈ 502,000")
# dot.node("B", "I2 brain MRI (WMH IDPs) available\n≈ 61,000")
# dot.node("C", "SA identified (Self-report / Hospital)\n≈ 13,000")

dot.node("BC", "Brain MRI available ≈ 46,000\nSA identified ≈ 13,000")

dot.node("D", "SA diagnosed before the V1 MRI\nN = 598")

# Branches
dot.node("E1", "After excluding SA cases with pre-index neuro dx\nN = 578")
dot.node("E2", "No pre-index exclusion")

# Matching steps
dot.node("F1", "Propensity score matching (1:10)\nControl index assigned\nPre-index neuro dx exclusion applied")
dot.node("F2", "Propensity score matching (1:10)")

# Final cohorts
dot.node("G1", "Primary cohort\nSA: N = 578\nControls: N = 5,672")
dot.node("G2", "Main sensitivity cohort\nSA: N = 598\nControls: N = 5,980")

# Align F1 and F2 horizontally
with dot.subgraph() as s:
    s.attr(rank='same')
    s.node("F1")
    s.node("F2")

# Align G1 and G2 horizontally
with dot.subgraph() as s:
    s.attr(rank='same')
    s.node("G1")
    s.node("G2")

# Exploratory longitudinal subcohort (derived from Primary cohort only)
dot.node("H1", "Exploratory longitudinal subcohort (V1→V2)\nFrom Primary cohort\nSA: N = 47\nControls: N = 528",
         style="rounded,dashed")


# ------------------------------------------------------
# 4. Define edges
# ------------------------------------------------------
# Branches from A
dot.edge("A", "BC")
# dot.edge("A", "C")

# Intersection: D = B ∩ C
dot.edge("BC", "D")
# dot.edge("C", "D")

# Continue paths
# dot.edge("D", "E1", xlabel ="Primary path", labeldistance = "2.5")
# dot.edge("D", "E2", label ="Main sensitivity path", labeldistance="2.5")
dot.edge("D", "E1")
dot.edge("D", "E2")

dot.edge("E1", "F1")
dot.edge("F1", "G1")

dot.edge("E2", "F2")
dot.edge("F2", "G2")

# Dashed edge to show derivation
dot.edge("G1", "H1", style="dashed")


# #Invisible node for alignment (optional)
# dot.node("AlignLeft", "", style="invis", width="0", height="0")

# ------------------------------------------------------
# 5. Save outputs (PNG, PDF, SVG, DOT, and high-res TIFF)
# ------------------------------------------------------
# Base paths
base_path = os.path.join(outdir, basename)
png_path = base_path + ".png"
pdf_path = base_path + ".pdf"
svg_path = base_path + ".svg"
tif_path = base_path + ".tiff"
dot_path = base_path + ".dot"

# Render outputs
dot.render(filename=base_path, format="png", cleanup=True)
dot.render(filename=base_path, format="pdf", cleanup=True)
dot.render(filename=base_path, format="svg", cleanup=True)
dot.render(filename=base_path, format="tiff", cleanup=True)

# Save DOT source once
with open(dot_path, "w", encoding="utf-8") as f:
    f.write(dot.source)

# Optional TIFF LZW compression
if Image is not None and os.path.exists(tif_path):
    try:
        im = Image.open(tif_path)
        if im.mode in ("P", "L"):
            im = im.convert("RGB")
        im.save(tif_path, format="TIFF", compression="tiff_lzw")
    except Exception as e:
        print(f"[Warning] Pillow LZW compression failed: {e}")

print("✅ Flowchart generated. Files saved in Flow_Chart folder:")
print(" -", png_path)
print(" -", pdf_path)
print(" -", svg_path)
print(" -", tif_path)
print(" -", dot_path)


# ------------------------------------------------------
# 6. Display inline in Jupyter
# ------------------------------------------------------
dot
