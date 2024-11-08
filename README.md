# TMS_RC_Model
This project taken brain structural connectomes as a large resistor capacitor network. The stimulation generated by TMS coils will provide the external charges into certain nodes at beginning. Then we use the electrical equations(KVL) at each node changing with time to build the ODE matrix. In the end we solve this ODEs and get nodes' voltage 

# TMS_RC_Model

## Step 1: TMS E-Field Simulation in SimNIBS

*This section explains how to position the TMS coil in SimNIBS at a specified MNI-space location and run the electric field (E-field) simulation.*

### 1.1 Define MNI Coordinates(Optional if you know the subject coordinates)
In this project we use the empircal TMS-fMRI data provided by Tik's group(Tik, Woletz et al. 2023) 
In their TMS experiments their targeting position is in mni-space with a coordinate: (-42, 28, 21) with an angle of 45 degree against the brain surface
The coil they used is a Figure-8 coil with a diameter of 70mm.
If your empircal data provide another mni coordinate, you can follow the same process by changing your coordinate in the above, and if you already knows the subject coordinates of the coil you can skip for step 1.3.
**Target MNI Coordinates**: `(-42, 28, 21)`
**Software Required**: [SimNIBS](https://simnibs.github.io/simnibs/) (v3.2 or later)

### 1.2 Build and choose the correct head model inside the SimNIBs(FEM software)


### 1.3 Convert MNI to Subject Coordinates
Run:
```shell
mni2subject_coords -m simnibs_examples/ernie/m2m_ernie/ -c -42 28 21  -o subject_TMS_coords.csv

