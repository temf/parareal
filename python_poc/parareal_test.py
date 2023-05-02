from parareal import solve_ivp_parareal
import matplotlib.pyplot as plt
import numpy as np


# Van der Pol ODE
def van_der_pol_fun(t, y):
    return [y[1], (1 - y[0] ** 2) * y[1] - y[0]]


trange = np.array([0, 40])
y0 = np.array([2.0, 0])

solC = "solve_ivp"
optionsC = {"rtol": 1e-1, "atol": 1e-1, "method": 'RK23'}

solF = 'solve_ivp'
optionsF = {"rtol": 1e-10, "atol": 1e-10, "method": 'RK45'}

paraOpt = {"n_time_windows": 5, "NIter": 5, "plotEach": True, "isParallel": True}

parallel_sol = solve_ivp_parareal(trange, y0, van_der_pol_fun, solC, optionsC, van_der_pol_fun, solF, optionsF, paraOpt)

plt.figure()
plt.plot(parallel_sol.t, parallel_sol.y[0, :], '-*', parallel_sol.t, parallel_sol.y[1, :])
plt.xlabel('t')
plt.ylabel('y')
plt.legend(['x(t)', 'dx/dx(t)'])
plt.show()
