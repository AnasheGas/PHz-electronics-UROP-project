# Petahertz Electronic Circuits with Optical Nanoantennas

## Unit model 
![unitcell-image](https://github.com/AdinaBechhofer/PHz-electronics-UROP-project/blob/main/New_design.png?raw=true)
The unit cell circuit model comprises of 2 components: RLC elements which simulate the electromagnetic response of the antenna to the optical excitation. The values of circuit elements can be fit by comparison to MEEP simulations, which solve for electromagnetics, but do not have a model for current emission. 

In the LTspice circuit model we include the current as a voltage controlled current source that is controlled by the voltage across the gap capacitor. The functional form of the current emission is derived from theory and verified by experiment. More on that in the Emission models section of this README. 

The connecting wires for nanoantennas are gold nanostrips which can be modeled as transmission lines. We were able to charecterize the input impedance and speed of the wave in the nanostrip to effectively describe the transmission line in LTspice.  When cascading devices to make larger intergrated circuits with multiple nanoantennas, we need to consider how an antenna couples to the transmission line. We define the coupling parameter $c_p$ as the quantity that splits the resistance and inductance of the antenna to before and after the wire connection. We can compare LT spice simulations of coupled antennas without the current emission to MEEP simulations and use that to nail down the $c_p$ parameter. 

## Memory cell 
The LTspice simulations of a memory cell are of coupled nano antennas with the current emission component in place. The capacitor on the second antenna collects charge and establishes a bias across the second antenna which shifts the operating point of the current emission function. When the memory capacitor is loaded, a small read pulse on the read antenna produces a large current spike that intergrates to total charge transfer. The state read is determained by intergrating the current in the read line during a read pulse and comparing it to the threshold. 

### Run LT spice
Open the `.asc` file in LTspice and make sure that the required `config.yaml` is in the same directory. Before running the simulation, you need to adjust the relative sensitivity of the nodal variables. This is because in the new unit system, one Volt is simulated as $ 1 \times 10^{-12}$ Volts, which is under the default noise tolerance. To update the tolerance click on the tool's menue and select the control panel with the hammer icon. There, change the absolute voltage tolerance to 1e-16 and the absolute current tolerance to 1e-18.

You can export traces to analyze in an external program by navigating to the trace panel, clicking the file menu, and choosing export data as text. 

Additional information on the scripts that do the simulation data analysis can be found in the readme files in the relevant locations. 

### Theory 
#### Units
SPICE by default take SI units. However, due to our as - fs timescales and aF capacitances, SPICE is unable to simulate the dynamics of the nanoantenna in the SI units probably because of a hard-coded tolerance. We ointroduce the following normalized units for spice simulation. 
* Capacitance: $F_n = 10^{12} F$ 
* Inductance: $H_n = 10^{18} H$
* Resistance: $\Omega_n = 10^{3} \Omega$
* Time: $s_n = fs$
* Voltage: $V_n = 10^{-12} V$
* Current: $A_n = 10^{-15} A$

LTspice discards changes that occur over time spans that are too short because those are seen as numerical noise. The normalized units allow us to stretch the time domain from femtoseconds to seconds, which has the required numerical stability that LTspice desires. 

#### Emission models 
Fowler-Nordheim, Shottkey. 

## TODOs:
Here is a [link](https://docs.google.com/document/d/1VlNZkuP5iME8qnPhX7BGFzTnaLO9OkyRKATvZGK03Pk/edit?usp=sharing) to a google doc with an outline for the UROP. 
