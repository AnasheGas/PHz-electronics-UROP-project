import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path


file_path_orig = Path(__file__).parent / "two_coupled_antennas_with_currents_ReadoutR"
data_original = pd.read_csv(file_path_orig, sep ="\t")

file_path_norm = Path(__file__).parent / "2_norm_antennasTest_ReadR"
data_normalized = pd.read_csv(file_path_norm, sep ="\t")



plt.figure(figsize=(10, 6))

plt.plot(
    data_original['time'], 
    data_original['I(readout)'], 
    label="Original Model",
)

plt.plot(
    data_normalized['time'], 
    data_normalized['I(readout1)'], 
    label="Normalized Model",
    linestyle="--"
)

plt.xlabel("Time")
plt.ylabel("Current")
plt.title("Current vs Time Comparison")
plt.legend()
plt.grid(True)

plt.show()
