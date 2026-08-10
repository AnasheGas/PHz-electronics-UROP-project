import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path


file_path_orig = Path(__file__).parent / "two_coupled_antennas_with_currents_ReadoutR"
data_original = pd.read_csv(file_path_orig, sep ="\t")

file_path_norm = Path(__file__).parent / "2_norm_antennasTest_ReadR"
data_normalized = pd.read_csv(file_path_norm, sep ="\t")



# 1. Scaled down figure size with high DPI for crisp printing on small posters
plt.figure(figsize=(10, 6), dpi=300)

# 2. Base line: Very thick, dark, solid midnight blue
plt.plot(
    data_original['time'], 
    data_original['I(readout)'], 
    label="Original Model",
    color="#0B2545",       # Deep navy blue
    linewidth=4.5,         # Extra thick to peek out from underneath
    zorder=1               # Places it in the background
)

# 3. Overlap line: Thinner, dashed, ultra-bright neon orange
plt.plot(
    data_normalized['time'], 
    data_normalized['I(readout1)'], 
    label="Updated Model",
    color="#FF6B35",       # High-contrast vibrant orange
    linewidth=2.0,         # Thinner so the navy borders it perfectly
    linestyle=(0, (3, 1.5)), # Custom tight dash pattern (3pt line, 1.5pt space)
    zorder=2               # Forces it on top
)

# 4. Large font sizes optimized for small posters
plt.xlabel("Time", fontsize=12, fontweight='bold')
plt.ylabel("Current", fontsize=12, fontweight='bold')
plt.title("Current vs Time Comparison", fontsize=14, fontweight='bold', pad=15)

# 5. Clean, high-contrast legend and grid
plt.legend(fontsize=10, loc="upper right", framealpha=0.9, facecolor="white", edgecolor="none")
plt.grid(True, linestyle=":", alpha=0.5, color="#888888")

# 6. Ensure labels don't get cut off
plt.tight_layout()
plt.savefig("antenna_comparisonFinal.png", bbox_inches='tight', dpi=300)
plt.show()