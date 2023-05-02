# Parareal & GetDP (Proof of Concept)

These Python files provide a proof of concept of a Parareal implementation using GetDP to integrate initial value problems.

To try it out, run the Python script `run_parareal_thermalLinear.py`. This requires to set the path to the GetDP executable correctly at the beginning of the script (`getdp_path`). No automatic searching for GetDP is implemented, unlike in the Matlab version. The output is verbose by default, for performance this should be switched off using the available options.


## File Structure 

The files are 

- `parareal.py`: implementation of Parareal in Python taking functions as coarse and solve integrators. 
- `parareal_test.py`: simple test and model usage of `parareal.py` using van der Pol ODE as example
- `getdp_runner.py`: implements `solve_ivp_getdp`, a function to solve initial value problems using GetDP with syntax inspired by `scipy`
- `getdp_read_presolution.py`: helper function for `getdp_runner`; reads a specific type of file (`.pre`) created by GetDP 
- `getdp_create_preresolution.py`: helper function for `getdp_runner`; creates a specific type of file (`.pre`) created by GetDP
- `run_parareal_thermalLinear.py`: main runner file using a simple thermal problem in GetDP to be solved with GetDP, which is located in `getdp_models`
- `getdp_models/`: folder containing the simple GetDP problem solved here, for more information on GetDP, see [the GetDP website](https://www.getdp.info)
- `plots/`: folder where `.png` of plots are saved if you use the default verbose implementation (nice for debugging)
