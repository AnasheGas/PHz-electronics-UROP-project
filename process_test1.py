# imports 
import os
import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from scipy import interpolate 

from scipy import optimize
import sympy 
from scipy import signal 
from scipy.signal import savgol_filter
import matplotlib.colors as mcolors
from scipy import integrate

from scipy.signal import butter, sosfiltfilt

# Useful constants 
consts = {
    "c": 2.99792458e8, 
    "e": 1.60217663e-19,
    "eps0": 8.854187817e-12
}

def TD_and_freq_info_from_two_antenna_spice_sweep(path, resample_t, two_source=True):
    with open(path, 'r') as file:
        data = file.readlines()

    table = {}
    head = data[0]
    for line in data[1:]:
        if "Step Information" in line:
            print(line)
            cp = line.split("=")[1]
            if "M" in cp:
                cp = str(float(cp.split('M')[0])*1e6)
            elif 'm' in cp:
                cp = str(float(cp.split('m')[0])*0.001)
            if 'K' in cp:
                cp = str(float(cp.split('K')[0])*1e3)
            elif 'n' in cp:
                cp = str(float(cp.split('n')[0])*1e-9)
            elif 'f' in cp:
                cp = str(float(cp.split('f')[0])*1e-15)
            elif 'p' in cp:
                cp = cp.split(' ')[0]
                if 'p' in cp:
                    cp = str(float(cp.split('p')[0])*1e-12)
            table[cp] = [head.split('\n')[0].split('\t')]
        else:
            table[cp].append(line.split('\n')[0].split('\t'))

    del data

    if "V(tpEmtrY,clctrY)" in list(table.values())[0][0]:
        inpY = "V(tpEmtrY,clctrY)"
    elif "V(clctrY,tpEmtrY)" in list(table.values())[0][0]:
        inpY = "V(clctrY,tpEmtrY)"
    elif "V(tpemtry)" in list(table.values())[0][0]:
        inpY = "V(tpemtry)"
    else:
        inpY = None

    
    if "V(tpEmtrX,clctrX)" in list(table.values())[0][0]:
        inpX = "V(tpEmtrX,clctrX)"
    elif "V(clctrX,tpEmtrX)" in list(table.values())[0][0]:
        inpX = "V(clctrX,tpEmtrX)"
    elif "V(tpemtrx)" in list(table.values())[0][0]:
        inpX = "V(tpemtrx)"
    else: 
        inpX = None
        
    if "V(btmEmtrY,clctrY)" in list(table.values())[0][0]:
        outpY = "V(btmEmtrY,clctrY)" 
    elif "V(clctrY,btmEmtrY)" in list(table.values())[0][0]:
        outpY = "V(clctrY,btmEmtrY)"
    elif "V(btmemtry)" in list(table.values())[0][0]:
        outpY = "V(btmemtry)"
    else:
        outpY = "V(btmemtry)"

    if "V(btmEmtrX,clctrX)" in list(table.values())[0][0]:
        outpX = "V(btmEmtrX,clctrX)"
    elif "V(clctrX,btmEmtrX)" in list(table.values())[0][0]:
        outpX = "V(clctrX,btmEmtrX)"
    else:
        outpX = "V(btmemtrx)"

    if "I(B1)" in list(table.values())[0][0]:
        b1 = "I(B1)"
        b1_currents = []
        b1C_list = []
    elif "I(YEmissionSrc)" in list(table.values())[0][0]: 
        b1 = "I(YEmissionSrc)"
        b1_currents = []
        b1C_list = []
    else:
        b1 = None

    if "I(B2)" in list(table.values())[0][0]:
        b2 = "I(B2)"
        b2_currents = []
        b2C_list = []
    if "I(XEmissionSrc)" in list(table.values())[0][0]:
        b2 = "I(XEmissionSrc)"
        b2_currents = []
        b2C_list = []
    else: 
        b2 = None

    if "I(R5)" in list(table.values())[0][0]:
        readoutR = "I(R5)"
        readout_currents = []
        readoutC_list = []
    elif "I(readout)" in list(table.values())[0][0]:
        readoutR = "I(readout)"
        readout_currents = []
        readoutC_list = []
    else: 
        readoutR = None

    if "I(C11)" in list(table.values())[0][0]:
        end_cap = "I(C11)"
        end_capC= []
    else:
        end_cap = None

    cp_list = []
    vin_list = []
    gapY_list = []
    gapX_list = []

    if two_source:
        vinY_list = []
        vinX_list = []

    for key in table:
        ## extract the phase by stripping white space
        cp = key.split(" ")[0]
        df = pd.DataFrame(table[key][1:], columns=table[key][0])
        #print(df.columns)
        df["time"] = df["time"].astype(float)
        if inpY is not None:
            df[inpY] = df[inpY].astype(float)
        if inpX is not None:
            df[inpX] = df[inpX].astype(float)
        cp_list.append(cp)
        #print(len(df["time"].to_numpy()))
        if inpY is not None:
            finY = interpolate.interp1d(df["time"].to_numpy()*1e-15, df[inpY].to_numpy().astype(float)*1e12)
        else: 
            finY = None
        if inpX is not None:
            finX = interpolate.interp1d(df["time"].to_numpy()*1e-15, df[inpX].to_numpy().astype(float)*1e12)
        else: 
            finX = None
        
        fgapY = interpolate.interp1d(df["time"].to_numpy()*1e-15, df[outpY].to_numpy().astype(float)*1e12)
        fgapX = interpolate.interp1d(df["time"].to_numpy()*1e-15, df[outpX].to_numpy().astype(float)*1e12)
        # TODO: now we select either x or y source, but we should allow both in the future 

        if b1 is not None:
            fb1_cur = interpolate.interp1d(df["time"].to_numpy()*1e-15, df[b1].to_numpy().astype(float)*1e15)
            fb1_c = integrate.simpson(df[b1].to_numpy().astype(float)*1e15, df["time"].to_numpy()*1e-15)
            b1C_list.append(fb1_c)
            Ib1_n = fb1_cur(resample_t)
            b1_currents.append(Ib1_n)

        if b2 is not None:
            fb2_cur = interpolate.interp1d(df["time"].to_numpy()*1e-15, df[b2].to_numpy().astype(float)*1e15)
            fb2_c = integrate.simpson(df[b2].to_numpy().astype(float)*1e15, df["time"].to_numpy()*1e-15)
            b2C_list.append(fb2_c)
            Ib2_n = fb2_cur(resample_t)
            b2_currents.append(Ib2_n)

        if readoutR is not None:
            readout_cur = interpolate.interp1d(df["time"].to_numpy()*1e-15, df[readoutR].to_numpy().astype(float)*1e15)
            readout_c = integrate.simpson(df[readoutR].to_numpy().astype(float)*1e15, df["time"].to_numpy()*1e-15)
            readoutC_list.append(readout_c)
            readout_n = readout_cur(resample_t)
            readout_currents.append(readout_n)

        if end_cap is not None:
            cap_c = integrate.simpson(df[end_cap].to_numpy().astype(float)*1e15, df["time"].to_numpy()*1e-15)
            end_capC.append(cap_c)
            
        if two_source:
            if inpY is not None:
                VinY_n = finY(resample_t)
                vinY_list.append(VinY_n)
            if inpX is not None:
                VinX_n = finX(resample_t)
                vinX_list.append(VinX_n)
        else:
            if inpY is not None:
                Vin_n = finY(resample_t)
            elif inpX is not None:
                Vin_n = finX(resample_t)
            vin_list.append(Vin_n)
            
        VgapY_n= fgapY(resample_t)
        VgapX_n= fgapX(resample_t)
        gapY_list.append(VgapY_n)
        gapX_list.append(VgapX_n)

    del table
    vin_fr_list = []
    if two_source:
        vinY_fr_list = []
        vinX_fr_list = []
    vgapY_fr_list = []
    vgapX_fr_list = []
   
    for i in range(len(cp_list)):
        if two_source:
            if inpY is not None:
                vinY_fr = np.fft.fftshift(np.fft.fft(vinY_list[i], n=1*len(resample_t)))
                vinY_fr_list.append(vinY_fr)
            if inpX is not None:
                vinX_fr = np.fft.fftshift(np.fft.fft(vinX_list[i], n=1*len(resample_t)))
                vinX_fr_list.append(vinX_fr)
        else:
            vin_fr = np.fft.fftshift(np.fft.fft(vin_list[i], n=1*len(resample_t)))
            vin_fr_list.append(vin_fr)
        vgapY_fr = np.fft.fftshift(np.fft.fft(gapY_list[i], n=1*len(resample_t)))
        vgapY_fr_list.append(vgapY_fr)
        vgapX_fr = np.fft.fftshift(np.fft.fft(gapX_list[i], n=1*len(resample_t)))
        vgapX_fr_list.append(vgapX_fr)
    freq = np.fft.fftshift(np.fft.fftfreq(4*len(resample_t), d=resample_t[1]-resample_t[0]))
    del df

    if two_source:
        if  b1 is not None and  b2 is not None:
            if readoutR is not None:
                if end_cap is not None:
                     return cp_list, vinY_list, vinX_list, gapY_list, gapX_list, vinY_fr_list, vinX_fr_list, vgapY_fr_list, vgapX_fr_list, freq, b1_currents, b1C_list, b2_currents, b2C_list, readout_currents, readoutC_list, end_capC
                return cp_list, vinY_list, vinX_list, gapY_list, gapX_list, vinY_fr_list, vinX_fr_list, vgapY_fr_list, vgapX_fr_list, freq, b1_currents, b1C_list, b2_currents, b2C_list, readout_currents, readoutC_list
            return cp_list, vinY_list, vinX_list, gapY_list, gapX_list, vinY_fr_list, vinX_fr_list, vgapY_fr_list, vgapX_fr_list, freq, b1_currents, b1C_list, b2_currents, b2C_list
        elif b1 is not None:
            return cp_list, vinY_list, vinX_list, gapY_list, gapX_list, vinY_fr_list, vinX_fr_list, vgapY_fr_list, vgapX_fr_list, freq, b1_currents, b1C_list
        elif b2 is not None:
            return cp_list, vinY_list, vinX_list, gapY_list, gapX_list, vinY_fr_list, vinX_fr_list, vgapY_fr_list, vgapX_fr_list, freq, b2_currents, b2C_list
        return cp_list, vinY_list, vinX_list, gapY_list, gapX_list, vinY_fr_list, vinX_fr_list, vgapY_fr_list, vgapX_fr_list, freq
    else:
        return cp_list, vin_list, gapY_list, gapX_list, vin_fr_list, vgapY_fr_list, vgapX_fr_list, freq
    

def get_simulation_runtime(log_file):
    with open(log_file, "r") as file:
        text = file.read()

    match = re.search(r"Total elapsed time:\s*([\d.]+)\s*seconds", text)

    if match:
        return float(match.group(1))
    else:
        return None
    
def save_result(filename, timestep, missing_electrons,  runtime):
    file_exists = os.path.isfile(filename)

    with open(filename, "a", newline="") as file:
        writer = csv.writer(file)
        if not file_exists:
            writer.writerow(["Timestep (fs)", "Average_Missing_Electrons",  "Runtime (s)"])


        writer.writerow([timestep, missing_electrons,  runtime])
    







# path_d1 = "two_coupled_antennas_with_currents_part1.txt"

# t_samp_d1 = np.arange(0, 574.99e-15, 0.01e-15)
# delay_list_d1, vinY_list_d1, vinX_list_d1, gapY_list_d1, gapX_list_d1, vinY_fr_list_d1, vinX_fr_list_d1, vgapY_fr_list_d1, vgapX_fr_list_d1, freq_d1, b1_currents_d1, b1C_list_d1, b2_currents_d1, b2C_list_d1, read_currents_d1, readC_list_d1 = TD_and_freq_info_from_two_antenna_spice_sweep(path_d1, t_samp_d1, two_source=True)
# delay_list_d1 = np.array(delay_list_d1, dtype=float)
# np.savetxt("extracted_charge_del1.txt", [delay_list_d1, readC_list_d1, b1C_list_d1, b2C_list_d1], fmt='%.18e')

# missing_electrons_d1 = []
# for i in range(len(delay_list_d1)):
#     print("Amplitude:", delay_list_d1[i])
#     print("Ey emission:",  b1C_list_d1[i], "C")
#     print("Ex emission:",  b2C_list_d1[i], "C")
#     print("Readout charge:",  readC_list_d1[i], "C")
#     print("Missing electrons:",  round((b1C_list_d1[i] - b2C_list_d1[i] - readC_list_d1[i])/consts["e"]))
#     missing_electrons_d1.append((b1C_list_d1[i] - b2C_list_d1[i] - readC_list_d1[i])/consts["e"])
    
# data_d1 = np.loadtxt("extracted_charge_del1.txt")
# delay_list_d1 = data_d1[0, :]
# readC_list_d1 = data_d1[1, :]
# b1C_list_d1 = data_d1[2, :]
# b2C_list_d1 = data_d1[3, :]


# path_d2 = "two_coupled_antennas_with_currents_part2.txt"

# t_samp_d2 = np.arange(0, 574.99e-15, 0.01e-15)
# delay_list_d2, vinY_list_d2, vinX_list_d2, gapY_list_d2, gapX_list_d2, vinY_fr_list_d2, vinX_fr_list_d2, vgapY_fr_list_d2, vgapX_fr_list_d2, freq_d2, b1_currents_d2, b1C_list_d2, b2_currents_d2, b2C_list_d2, read_currents_d2, readC_list_d2 = TD_and_freq_info_from_two_antenna_spice_sweep(path_d2, t_samp_d2, two_source=True)
# delay_list_d2 = np.array(delay_list_d2, dtype=float)
# np.savetxt("extracted_charge_del2.txt", [delay_list_d2, readC_list_d2, b1C_list_d2, b2C_list_d2], fmt='%.18e')

# missing_electrons_d2 = []
# for i in range(len(delay_list_d2)):
#     print("Amplitude:", delay_list_d2[i])
#     print("Ey emission:",  b1C_list_d2[i], "C")
#     print("Ex emission:",  b2C_list_d2[i], "C")
#     print("Readout charge:",  readC_list_d2[i], "C")
#     print("Missing electrons:",  round((b1C_list_d2[i] - b2C_list_d2[i] - readC_list_d2[i])/consts["e"]))
#     missing_electrons_d2.append((b1C_list_d2[i] - b2C_list_d2[i] - readC_list_d2[i])/consts["e"])

# data_d2 = np.loadtxt("extracted_charge_del2.txt")
# delay_list_d2 = data_d2[0, :]
# readC_list_d2 = data_d2[1, :]
# b1C_list_d2 = data_d2[2, :]
# b2C_list_d2 = data_d2[3, :]



# path_d3 = "two_coupled_antennas_with_currents_part3.txt"

# t_samp_d3 = np.arange(0, 574.99e-15, 0.01e-15)
# delay_list_d3, vinY_list_d3, vinX_list_d3, gapY_list_d3, gapX_list_d3, vinY_fr_list_d3, vinX_fr_list_d3, vgapY_fr_list_d3, vgapX_fr_list_d3, freq_d3, b1_currents_d3, b1C_list_d3, b2_currents_d3, b2C_list_d3, read_currents_d3, readC_list_d3 = TD_and_freq_info_from_two_antenna_spice_sweep(path_d3, t_samp_d3, two_source=True)
# delay_list_d3 = np.array(delay_list_d3, dtype=float)
# np.savetxt("extracted_charge_del3.txt", [delay_list_d3, readC_list_d3, b1C_list_d3, b2C_list_d3], fmt='%.18e')

# missing_electrons_d3 = []
# for i in range(len(delay_list_d3)):
#     print("Amplitude:", delay_list_d3[i])
#     print("Ey emission:",  b1C_list_d3[i], "C")
#     print("Ex emission:",  b2C_list_d3[i], "C")
#     print("Readout charge:",  readC_list_d3[i], "C")
#     print("Missing electrons:",  round((b1C_list_d3[i] - b2C_list_d3[i] - readC_list_d3[i])/consts["e"]))
#     missing_electrons_d3.append((b1C_list_d3[i] - b2C_list_d3[i] - readC_list_d3[i])/consts["e"])

# data_d3 = np.loadtxt("extracted_charge_del3.txt")
# delay_list_d3 = data_d3[0, :]
# readC_list_d3 = data_d3[1, :]
# b1C_list_d3 = data_d3[2, :]
# b2C_list_d3 = data_d3[3, :]



# path_d4 = "two_coupled_antennas_with_currentspart4.txt"

# t_samp_d4 = np.arange(0, 574.99e-15, 0.01e-15)
# delay_list_d4, vinY_list_d4, vinX_list_d4, gapY_list_d4, gapX_list_d4, vinY_fr_list_d4, vinX_fr_list_d4, vgapY_fr_list_d4, vgapX_fr_list_d4, freq_d4, b1_currents_d4, b1C_list_d4, b2_currents_d4, b2C_list_d4, read_currents_d4, readC_list_d4 = TD_and_freq_info_from_two_antenna_spice_sweep(path_d4, t_samp_d4, two_source=True)
# delay_list_d4 = np.array(delay_list_d4, dtype=float)
# np.savetxt("extracted_charge_del4.txt", [delay_list_d4, readC_list_d4, b1C_list_d4, b2C_list_d4], fmt='%.18e')

# missing_electrons_d4 = []
# for i in range(len(delay_list_d4)):
#     print("Amplitude:", delay_list_d4[i])
#     print("Ey emission:",  b1C_list_d4[i], "C")
#     print("Ex emission:",  b2C_list_d4[i], "C")
#     print("Readout charge:",  readC_list_d4[i], "C")
#     print("Missing electrons:",  round((b1C_list_d4[i] - b2C_list_d4[i] - readC_list_d4[i])/consts["e"]))
#     missing_electrons_d4.append((b1C_list_d4[i] - b2C_list_d4[i] - readC_list_d4[i])/consts["e"])
    
# data_d4 = np.loadtxt("extracted_charge_del4.txt")
# delay_list_d4 = data_d4[0, :]
# readC_list_d4 = data_d4[1, :]
# b1C_list_d4 = data_d4[2, :]
# b2C_list_d4 = data_d4[3, :]

# path_d5 = "two_coupled_antennas_with_currentspart5.txt"

# t_samp_d5 = np.arange(0, 574.99e-15, 0.01e-15)
# delay_list_d5, vinY_list_d5, vinX_list_d5, gapY_list_d5, gapX_list_d5, vinY_fr_list_d5, vinX_fr_list_d5, vgapY_fr_list_d5, vgapX_fr_list_d5, freq_d5, b1_currents_d5, b1C_list_d5, b2_currents_d5, b2C_list_d5, read_currents_d5, readC_list_d5 = TD_and_freq_info_from_two_antenna_spice_sweep(path_d5, t_samp_d5, two_source=True)
# delay_list_d5 = np.array(delay_list_d5, dtype=float)
# np.savetxt("extracted_charge_del5.txt", [delay_list_d5, readC_list_d5, b1C_list_d5, b2C_list_d5], fmt='%.18e')

# missing_electrons_d5 = []
# for i in range(len(delay_list_d5)):
#     print("Amplitude:", delay_list_d5[i])
#     print("Ey emission:",  b1C_list_d5[i], "C")
#     print("Ex emission:",  b2C_list_d5[i], "C")
#     print("Readout charge:",  readC_list_d5[i], "C")
#     print("Missing electrons:",  round((b1C_list_d5[i] - b2C_list_d5[i] - readC_list_d5[i])/consts["e"]))
#     missing_electrons_d5.append((b1C_list_d5[i] - b2C_list_d5[i] - readC_list_d5[i])/consts["e"])
    
# data_d5 = np.loadtxt("extracted_charge_del5.txt")
# delay_list_d5 = data_d5[0, :]
# readC_list_d5 = data_d5[1, :]
# b1C_list_d5 = data_d5[2, :]
# b2C_list_d5 = data_d5[3, :]



# path_d6 = "two_coupled_antennas_with_currentspart6.txt"

# t_samp_d6 = np.arange(0, 574.99e-15, 0.01e-15)
# delay_list_d6, vinY_list_d6, vinX_list_d6, gapY_list_d6, gapX_list_d6, vinY_fr_list_d6, vinX_fr_list_d6, vgapY_fr_list_d6, vgapX_fr_list_d6, freq_d6, b1_currents_d6, b1C_list_d6, b2_currents_d6, b2C_list_d6, read_currents_d6, readC_list_d6 = TD_and_freq_info_from_two_antenna_spice_sweep(path_d6, t_samp_d6, two_source=True)
# delay_list_d6 = np.array(delay_list_d6, dtype=float)
# np.savetxt("extracted_charge_del6.txt", [delay_list_d6, readC_list_d6, b1C_list_d6, b2C_list_d6], fmt='%.18e')

# missing_electrons_d6 = []
# for i in range(len(delay_list_d6)):
#     print("Amplitude:", delay_list_d6[i])
#     print("Ey emission:",  b1C_list_d6[i], "C")
#     print("Ex emission:",  b2C_list_d6[i], "C")
#     print("Readout charge:",  readC_list_d6[i], "C")
#     print("Missing electrons:",  round((b1C_list_d6[i] - b2C_list_d6[i] - readC_list_d6[i])/consts["e"]))
#     missing_electrons_d6.append((b1C_list_d6[i] - b2C_list_d6[i] - readC_list_d6[i])/consts["e"])
    
# data_d6 = np.loadtxt("extracted_charge_del6.txt")
# delay_list_d6 = data_d6[0, :]
# readC_list_d6 = data_d6[1, :]
# b1C_list_d6 = data_d6[2, :]
# b2C_list_d6 = data_d6[3, :]



# path_d7 = "two_coupled_antennas_with_currentspart7.txt"

# t_samp_d7 = np.arange(0, 574.99e-15, 0.01e-15)
# delay_list_d7, vinY_list_d7, vinX_list_d7, gapY_list_d7, gapX_list_d7, vinY_fr_list_d7, vinX_fr_list_d7, vgapY_fr_list_d7, vgapX_fr_list_d7, freq_d7, b1_currents_d7, b1C_list_d7, b2_currents_d7, b2C_list_d7, read_currents_d7, readC_list_d7 = TD_and_freq_info_from_two_antenna_spice_sweep(path_d7, t_samp_d7, two_source=True)
# delay_list_d7 = np.array(delay_list_d7, dtype=float)
# np.savetxt("extracted_charge_del7.txt", [delay_list_d7, readC_list_d7, b1C_list_d7, b2C_list_d7], fmt='%.18e')

# missing_electrons_d7 = []
# for i in range(len(delay_list_d7)):
#     print("Amplitude:", delay_list_d7[i])
#     print("Ey emission:",  b1C_list_d7[i], "C")
#     print("Ex emission:",  b2C_list_d7[i], "C")
#     print("Readout charge:",  readC_list_d7[i], "C")
#     print("Missing electrons:",  round((b1C_list_d7[i] - b2C_list_d7[i] - readC_list_d7[i])/consts["e"]))
#     missing_electrons_d7.append((b1C_list_d7[i] - b2C_list_d7[i] - readC_list_d7[i])/consts["e"])
    
# data_d7 = np.loadtxt("extracted_charge_del7.txt")
# delay_list_d7 = data_d7[0, :]
# readC_list_d7 = data_d7[1, :]
# b1C_list_d7 = data_d7[2, :]
# b2C_list_d7 = data_d7[3, :]




# path_d8 = "two_coupled_antennas_with_currentspart8.txt"

# t_samp_d8 = np.arange(0, 574.99e-15, 0.01e-15)
# delay_list_d8, vinY_list_d8, vinX_list_d8, gapY_list_d8, gapX_list_d8, vinY_fr_list_d8, vinX_fr_list_d8, vgapY_fr_list_d8, vgapX_fr_list_d8, freq_d8, b1_currents_d8, b1C_list_d8, b2_currents_d8, b2C_list_d8, read_currents_d8, readC_list_d8 = TD_and_freq_info_from_two_antenna_spice_sweep(path_d8, t_samp_d8, two_source=True)
# delay_list_d8 = np.array(delay_list_d8, dtype=float)
# np.savetxt("extracted_charge_del8.txt", [delay_list_d8, readC_list_d8, b1C_list_d8, b2C_list_d8], fmt='%.18e')

# missing_electrons_d8 = []
# for i in range(len(delay_list_d8)):
#     print("Amplitude:", delay_list_d8[i])
#     print("Ey emission:",  b1C_list_d8 [i], "C")
#     print("Ex emission:",  b2C_list_d8[i], "C")
#     print("Readout charge:",  readC_list_d8[i], "C")
#     print("Missing electrons:",  round((b1C_list_d8[i] - b2C_list_d8[i] - readC_list_d8[i])/consts["e"]))
#     missing_electrons_d8.append((b1C_list_d8[i] - b2C_list_d8[i] - readC_list_d8[i])/consts["e"])
    
# data_d8 = np.loadtxt("extracted_charge_del8.txt")
# delay_list_d8 = data_d8[0, :]
# readC_list_d8 = data_d8[1, :]
# b1C_list_d8 = data_d8[2, :]
# b2C_list_d8 = data_d8[3, :]










# readout_all = np.concatenate((readC_list_d1, readC_list_d2, readC_list_d3, readC_list_d4, readC_list_d5, readC_list_d6, readC_list_d7, readC_list_d8))

# emission_y_all = np.concatenate((b1C_list_d1, b1C_list_d2, b1C_list_d3, b1C_list_d4, b1C_list_d5, b1C_list_d6, b1C_list_d7, b1C_list_d8))

# emission_x_all = np.concatenate((b2C_list_d1, b2C_list_d2, b2C_list_d3, b2C_list_d4, b2C_list_d5, b2C_list_d6, b2C_list_d7, b2C_list_d8))

# delays_all = np.concatenate((delay_list_d1, delay_list_d2, delay_list_d3, delay_list_d4, delay_list_d5, delay_list_d6, delay_list_d7, delay_list_d8 ))

# sort_ind = np.argsort(delays_all)

# delays_all = delays_all[sort_ind]
# readout_all = readout_all[sort_ind]
# emission_y_all = emission_y_all[sort_ind]
# emission_x_all = emission_x_all[sort_ind]


# fig, ax = plt.subplots(1, 1, figsize=(12, 5))
# ax.plot(delays_all, readout_all, label="Readout")
# ax.set_xlim([0, 40])
# ax.set_ylim([-9.2e-17, -4e-17])
# ax.set_xlabel(r"$\tau$ (fs)", fontsize=14)
# ax.set_ylabel("Readout charge (C)", fontsize=14)
# ax.tick_params(axis='both', labelsize=14)
# plt.savefig("last_fig_cross_correlation_4.svg", format="SVG")
# plt.show()


# avg_missing_electrons_d8 = np.mean(missing_electrons_d8)
# runtime = get_simulation_runtime("two_coupled_antennas_with_currents.log")
# timestep = 0.13  # fs
# print("Simulation runtime:", runtime, "seconds")
# print("Average missing electrons:", avg_missing_electrons_d8)
# print("Timestep:", timestep, "fs")
# save_result("missing_electrons_vs_timestep_results.csv", timestep, avg_missing_electrons_d8 , runtime)



df = pd.read_csv("missing_electrons_vs_timestep_results.csv")

# Create the plot
plt.figure(figsize=(8, 5))
plt.plot(
    df["Timestep (fs)"],
    df["Runtime (s)"],
    marker="o"
)

plt.xlabel("Timestep (fs)")
plt.ylabel("Runtime (s)")
plt.title("Effect of Timestep on Runtime")
plt.grid(True)

plt.show()


















# dt = t_samp_d1[1] - t_samp_d1[0]    # 0.01 fs
# fs = 1.0 / dt                         # ~100 PHz
# fc = 100e12   # [Hz] – adjust between 30–80 THz if needed
# fig, axs = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
# #plt.plot(amp_list_d10, b1C_lis__d1, label="Y emission")
# for i in range(len(delay_list_d1)):
#     order = 6
#     sos = butter(order, fc, btype='low', fs=fs, output='sos')
#     gapX_butter = sosfiltfilt(sos, gapX_list_d1[i]) 
#     gapY_butter = sosfiltfilt(sos, gapY_list_d1[i]) 
#     # ind = np.argmax(t_samp_v1 >1.50e-13)
#     # avg_gapX_voltage = np.mean(gapX_list_d10[i][ind:])
#     axs[1].plot(t_samp_d1, gapX_butter, label="Ex {} V".format(1e12*delay_list_d1[i]))
#     axs[0].plot(t_samp_d1, gapY_butter, label="Ex {} V".format(1e12*delay_list_d1[i]))
#     #plt.plot(t_samp_d10, gapX_lis__d1[i], '--', alpha=0.5)
# #plt.xlim([0.34e-13, 0.48e-13])
# #axs[0].legend(loc=3)
# #axs[1].legend(loc=2)
# axs[0].set_xlabel("Time (s)")
# axs[1].set_xlabel("Time (s)")
# axs[0].set_ylabel("Y gap")
# axs[1].set_ylabel("X gap")
# # HOLD plt.show()

# plt.plot(delay_list_d1, readC_list_d1, label="Readout")
# plt.plot(delay_list_d1, b1C_list_d1, label="Ey emission")
# plt.plot(delay_list_d1, b2C_list_d1, label="Ex emission")
# plt.legend()
# plt.xlabel(r"$\tau$ (fs)")
# # HOLD plt.show()

# for i in range(len(delay_list_d1)):
#     plt.plot(t_samp_d1, b2_currents_d1[i], label="Ex : {} V".format(round(delay_list_d1[i], 3)))
# #plt.xlim([1.8e-13, 2e-13])
# #plt.legend()
# plt.xlabel("Time (s)")
# plt.ylabel("X antenna tunnerling current (A)")
# # HOLD plt.show()


# for i in range(len(delay_list_d1)):
#     plt.plot(t_samp_d1, read_currents_d1[i], label="Ex amp {} V".format(round(delay_list_d1[i]*1e12, 3)))
# plt.xlabel("Time (s)")
# plt.ylabel("Inst readout current (A)")
# #plt.legend()
# # HOLD plt.show()

# plt.hist(missing_electrons_d1)
# # HOLD plt.show()












