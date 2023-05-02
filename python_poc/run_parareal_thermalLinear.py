from time import sleep
from matplotlib import pyplot as plt
import numpy as np
from getdp_create_preresolution import getdp_create_preresolution
from getdp_runner import solve_ivp_getdp
from parareal import solve_ivp_parareal


getdp_path ="C:/Program Files/CERNGetDP/cerngetdp_v1.05/getdp.exe"

t_init = 0 
t_end = 1


# Fine model
solver_f = "solve_ivp_getdp"
fun_f = "getdp_models/thermal_3d_linear/thermal.pro"
getdp_opts_f = {
              "exe": getdp_path,
              "mesh": "getdp_models/thermal_3d_linear/thermal.msh",
              "time_step": 0.0001,
              "verbose": "0"
              }

# Coarse model
solver_c = "solve_ivp_getdp"
fun_c = "getdp_models/thermal_3d_linear/thermal.pro"
getdp_opts_c = {
              "exe": getdp_path,
              "mesh": "getdp_models/thermal_3d_linear/thermal.msh",
              "time_step": 0.05,
              "verbose": "5"
              }


# Parareal options 
para_opts = {"n_time_windows": 5, "NIter": 5, "plotEach": True, "isParallel": True}

t_range = np.array([t_init, t_end])

numdof = getdp_create_preresolution(fun_c, getdp_opts_c)
init = np.zeros(numdof)

# Parareal solution
parallel_sol = solve_ivp_parareal(t_range, init, fun_c, solver_c, getdp_opts_c, fun_f, solver_f, getdp_opts_f, para_opts)

# sequential solution
seq_sol = solve_ivp_getdp(fun_f, t_range, init, **getdp_opts_f)

plt.figure()
plt.plot(parallel_sol.t, np.max(parallel_sol.y, axis=0), '*', label='parallel')
plt.plot(seq_sol.t, np.max(seq_sol.y, axis=0), label='seq')
plt.legend()
plt.xlabel("Time")
plt.ylabel("Solution")
plt.show()