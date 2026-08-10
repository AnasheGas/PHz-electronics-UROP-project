import os
from os import path
import subprocess
from zipfile import Path
import PyLTSpice
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from scipy.integrate import simpson


netlist = """
O1 N001 N008 N002 N005 Stub1
C1 N002 N005 0.468µ
O2 N003 N006 N004 N007 roundTline
C2 0 N011 0.468µ
readout1 0 N011 {Rread}
O3 N010 N009 0 N011 Stub2
X§X3 N004 N007 N001 N008 single_antenna_normalized params: amp = {ampY} del = 105 phase = {phaseTop} sigdur = {sigdurTop} d = {dTop}
X§X4 N006 N003 N010 N009 single_antenna_normalized params: amp = {ampX} del = {105 + delt} phase = {phaseTop} sigdur = {sigdurTop} d ={dTop}

* block symbol definitions
.subckt single_antenna_normalized E2 C2 E1 C1
C1 N002 N004 {ceff*Cbody_norm/(cp+1)}
Cgap1 clctr1 btmEmtr1 {Cgap_norm}
C2 clctr1 N002 {Cwirebody_norm}
Rbod1 N002 N003 {0.5*(cp+1)*Rbody_norm}
R1 N004 N005 {0.5*(1-cp)*Rbody_norm}
L1 N003 N004 {0.5*(cp+1)*Lbody_norm} Rser=1m
L2 N005 btmEmtr1 {0.5*(1-cp)*Lbody_norm} Rser=1m
R2 N001 tpEmtr1 {Rin_norm}
C3 N002 N001 {1*Cin_norm}
C4 N004 btmEmtr1 {ceff*Cbody_norm/(1-cp)}
C5 N002 N004 {ceff*Cbody_norm/(cp+1)}
C6 N004 btmEmtr1 {ceff*Cbody_norm/(1-cp)}
O1 N004 clctr1 E1 C1 FudgeLine
O2 E2 C2 N004 clctr1 FudgeLine
B1 tpEmtr1 clctr1 V= gaussian_pulse({amp}, {del}, {sigdur} , {freq}, {phase})
B§EmissionSrc1 btmEmtr1 clctr1 I = 0.99*ifn(2.4E-16, 1.5414u, 6.8309G, 5.3, -1*V(btmEmtr1, clctr1), d, 7) + 0.01*is(1.6E-16, 1, 290, 5.3, -1*V(btmEmtr1 , clctr1), d, 7)
.include config_realunits.yaml
.param cp =-0.71
.func gaussian_pulse(a, b, c, freq, phase)  a*exp(-square((time- b)/c)/2)*cos(2*pi*freq*(time - b)-phase)
.param del = 30
.param freq =  0.28201
.model FudgeLine LTRA(len={0.075u} R={200G} L=5.791T C=13.3)
.param ceff = 0.57
.func fn_v(f) 1 - f  + f*ln(f)/6
.func get_f(E) pwr(q, 3)/(4*pi*epsilon1*pwr(phi, 2))*E
.func ifn(A, a, b, phi, vc, d, g) -1E-15*A*a*square(g*Efield(vc*1E+12, d))/phi*exp(limit(-b*pwr(phi, 3/2)/(g*Efield(vc*1E+12, d))*fn_v(get_f(g*Efield(vc*1E+12, d))), -160, 160))
.func Efield(v,d) abs(max(v/d, 1E-18))
.func is(A, alph, T, phi, vc, d, g) -1E-15*A*alph*square(T)*exp(limit((-phi +0.5*sqrt(g*Efield(vc*1E+12, d)/(pi*epsilon0)))/(Kb*T), -160, 160))
.param Cbody_norm = {Cbody}*1e12
.param Cin_norm = {Cin}*1e12
.param Cgap_norm = {Cgap}*1e12
.param Cwirebody_norm = {Cwirebody}*1e12
.param Lbody_norm = {Lbody}*1e18
.param Rbody_norm = {Rbody}*1e3
.param Rin_norm = {Rin}*1e3
.ends single_antenna_normalized

.tran 0 260 120 0.007
.step param delt 17 18.5 0.25
.param Rread = 448k
.param leff=1.2
.include config_coupled_subckt.yaml
.model roundTline LTRA(len={1.12*1.513u} R={200G} L=5.791T C=13.3)
.model Stub1 LTRA(len={1.14*0.2u} R={200G} L=5.791T C=13.3)
.model Stub2 LTRA(len={leff*0.2u} R={200G} L=5.791T C=13.3)
.opt method=gear abstol=1e-16 reltol=1e-15 gmin=1e-12 Trtol=0.02
.opt chgtol=5e-18 itl1=300 itl2=450 itl4=500 itl5= 800
.param ampY = Ey*dTop
.param ampX = Ex*dTop
.param phaseTop = 0.72*pi
.param sigdurTop = 5.1
.step param dTop LIST 40n 60n 150n 300n
.param Ex = 1.73E-4
.param Ey = 3.54E-4
.backanno
.end



"""

path = r"C:\Users\anash\OneDrive\Desktop\QNN REPO\PHz-electronics-UROP-project\2026_nanoantenna_subckt_files\2_norm_antennasTest.net"
with open(path, "w") as f:
    f.write(netlist)

LTSPICE = r"C:\Users\anash\AppData\Local\Programs\ADI\LTspice\LTspice.exe"



# result = subprocess.run([LTSPICE, "-b", "-Run", path], check=True, capture_output=True, text=True) 





raw = PyLTSpice.RawRead(r"C:\Users\anash\OneDrive\Desktop\QNN REPO\PHz-electronics-UROP-project\2026_nanoantenna_subckt_files\2_norm_antennasTest.raw")


# time_trace = raw.get_trace("time")
# current_trace = raw.get_trace("I(readout1)")

# delt_values = np.arange(17, 18.5 + 0.25, 0.25)
# dTop_values = [40, 60, 150, 300]

# plt.figure(figsize=(8,6))

# for i, dTop in enumerate(dTop_values):

#     charge = []

#     for j, delt in enumerate(delt_values):

#         step = i * len(delt_values) + j

#         t = time_trace.get_wave(step=step)
#         I = current_trace.get_wave(step=step)

#         # Integrate current over time
#         Q = simpson(I, x=t)

#         charge.append(Q)

#     plt.plot(delt_values,
#              charge,
#              marker='o',
#              linewidth=2,
#              label=f"dTop = {dTop} nm")

# plt.xlabel("delt (fs)")
# plt.ylabel("Integrated Current (C)")
# plt.title("Readout Charge vs. CEP Delay")
# plt.grid(True)
# plt.legend()
# plt.show()


current_trace = raw.get_trace("I(readout1)")

delt_values = np.arange(17, 18.5 + 0.25, 0.25)
dTop_values = [40, 60, 150, 300]

plt.figure(figsize=(8,6))

for i, dTop in enumerate(dTop_values):

    peak_current = []

    for j, delt in enumerate(delt_values):

        step = i * len(delt_values) + j

        I = current_trace.get_wave(step=step)

        # Peak current
        peak = np.max(I)

        peak_current.append(peak)

    plt.plot(delt_values,
             peak_current,
             marker='o',
             linewidth=2,
             label=f"dTop = {dTop} nm")

plt.xlabel("delt")
plt.ylabel("Peak Current through Readout Resistor")
plt.title("Peak Readout Current vs. delt")
plt.grid(True)
plt.legend()
plt.show()




# t = raw.get_trace("time").get_wave()
# readR = raw.get_trace("I(readout1)").get_wave()


# plt.plot(t, readR)
# plt.xlabel("t"); plt.ylabel("I(readout1)")
# plt.show()

