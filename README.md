# TMS_RC_Model
This project taken brain structural connectomes as a large resistor capacitor network. The stimulation generated by TMS coils will provide the external charges into certain nodes at beginning. Then we use the electrical equations(KVL) at each node changing with time to build the ODE matrix. In the end we solve this ODEs and get nodes' voltage 