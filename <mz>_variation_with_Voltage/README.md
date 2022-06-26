### In this code, the Input voltage is varied from -4 to 4 volts. Due to the numerical accuracy issues, voltage near 0 are avoided.
### For each voltage, the simulation is run for N_runs times. For each run among 30 runs, average mz is calculated and stored in an array of length N_runs, and the final average mz for that voltage is calculated by taking the average of the average mz s stored in the array of length N_runs.
### For accuracy, N_runs can be increased further.
### Issue of this code remains for choosing delta_t. To calculate the average mz s at the voltages near 0, we need to choose delta_t properly. 
### This modification is under process...
