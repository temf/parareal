from scipy.integrate import solve_ivp
import scipy
from getdp_runner import solve_ivp_getdp
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
import multiprocessing
from joblib import Parallel, delayed
from datetime import datetime
import os

plt.style.use('seaborn')


def regenerate_plot(ax, sol, n_time_windows, trange_time_windows, colormap, output_folder):
    """
    Redraws the matplotlib plot in ax using solution sol.
    :param ax: axis object to update
    :param sol: ODE solution object to get solution from
    :param n_time_windows: number of time intervals of parareal
    :param trange_time_windows: the time range of each time interval
    :param colormap: the colormap to be used
    """
    for i in range(n_time_windows):
        ax.plot(sol[i].t, np.max(sol[i].y, axis=0), 
                color=colormap[i])
        # draw vertical lines to mark time intervals, redrawn every time since remove_lines removes all lines
        if i != 0:
            ax.axvline(trange_time_windows[i, 0], ls='--', color='grey')

    plt.pause(0.1)  # needed for the dynamic update
    now = datetime.now().strftime("%H_%M_%S")
    plt.savefig(os.path.join(output_folder, "plots", now))


def remove_lines(ax):
    """
    Remove all lines in the axis object ax.
    :param ax: axis object to remove lines from
    """
    while ax.lines:
        ax.lines[0].remove()


def solve_ivp_parareal(trange, y0, fun_c, solver_c, optionsC, fun_f, solver_f, optionsF, paraOpt):
    """
    Solve initial value problem using Parareal.
    :param trange: numpy.array([starting time, end time])
    :param y0: numpy array with initial solution
    :param fun_c: function to integrate for coarse integrator
    :param solver_c: method for coarse integrator
    :param optionsC: options for coarse integrator
    :param fun_f: function to integrate for fine integrator (i.e., the .pro file)
    :param solver_f: method for fine integrator
    :param optionsF: options for fine integrator
    :param paraOpt: Parareal options
    """
    
    # get options
    if "n_time_windows" in paraOpt:
        n_time_windows = paraOpt["n_time_windows"]
    else:
        # default
        n_time_windows = 8

    if "NIter" in paraOpt:
        NIter = paraOpt["NIter"]
    else:
        # default
        NIter = 20

    if "NProc" in paraOpt:
        NProc = paraOpt["NProc"]
    else:
        # default
        NProc = multiprocessing.cpu_count()

    if "plotEach" in paraOpt:
        plotEach = paraOpt["plotEach"]
    else:
        # default
        plotEach = False

    if "isParallel" in paraOpt:
        isParallel = paraOpt["isParallel"]
    else:
        # default
        isParallel = False

    if "abs_tol" in paraOpt:
        abs_tol = paraOpt["abs_tol"]
    else:
        # default
        abs_tol = 1E-10
          
    if "rel_tol" in paraOpt:
        rel_tol = paraOpt["rel_tol"]
    else:
        # default
        rel_tol = 1E-10

    if "output_folder" in paraOpt:
        output_folder = paraOpt["output_folder"]
    else:
        # default
        output_folder = os.getcwd()

    if (y0.shape[0] == 1 and y0.shape[1] != 1):
        y0 = np.transpose(y0)

    if plotEach:
        colors = plt.cm.viridis(np.linspace(0, 0.9, n_time_windows))  # do not use complete yellow spectrum of viridis

        fig, ax = plt.subplots(2)
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle("Evolution of coarse and fine solution")
        ax[0].set_title("Coarse Solution")
        ax[1].set_title("Fine Solution")

        ax[0].set_xlabel("Time in (s)")
        ax[0].set_ylabel("Solution in a.u.")

        ax[1].set_xlabel("Time in (s)")
        ax[1].set_ylabel("Solution in a.u.")

        ax[0].plot([],[])
        ax[1].plot([],[])
        plt.pause(1)


    y_tot = np.zeros((n_time_windows + 1, len(y0)))
    err = np.zeros((NIter, 1))

    y_cnew = np.zeros((n_time_windows + 1, len(y0)))  # matrix containing the new coarse solution
    y_f0 = np.zeros((n_time_windows + 1, len(y0)))  # matrix containing the fine solution for proc i+1 at t_0
    y_fend = np.zeros((n_time_windows + 1, len(y0)))  # matrix containing the fine solution for proc i at t_end
    y_tot = np.zeros((n_time_windows + 1, len(y0)))  # matrix containing the initial values
    sol_c = []
    sol_f = []

    y_tot[0, :] = y0

    if len(trange) == 2: # same time step size for each 
        equiTRange = np.linspace(trange[0], trange[1], n_time_windows + 1)
        idx = np.array((range(0, n_time_windows), range(1, n_time_windows + 1))).transpose()
        trange_time_windows = equiTRange[idx]
    else: 
        trange_time_windows = np.array(list(zip(trange[:-1], trange[1:])))
    
    k = 0

    while ((k < 1 or err[k - 1] > 1)):
        sol_c = []
        sol_f = []
        k = k + 1

        print("\nPARAREAL: START %d/%d ITERATION \n" % (k, NIter))
        print('\nPARAREAL: start coarse grid\n')

        for i in range(0, n_time_windows):

            print('Parareal: coarse integration in t = [%f, %f]' % (trange_time_windows[i, 0], trange_time_windows[i, 1]))

            # solve_ivp unfortunately does not allow args..
            sol_c.append(globals()[solver_c](fun_c, trange_time_windows[i, :], y_tot[i], **optionsC))
            # sol_c.append(eval(sol_c))
            # ATTENTION deepcopy needed otherwise points to the same list object
            y_cold = deepcopy(y_cnew[i + 1, :])

            y_cnew[i + 1, :] = sol_c[i].y[:, -1]

            y_tot[i + 1, :] = y_cnew[i + 1, :] + y_fend[i + 1, :] - y_cold

        if plotEach:
            remove_lines(ax[0])
            regenerate_plot(ax[0], sol_c, n_time_windows, trange_time_windows, colors, output_folder)

        print('\nPARREAL: start fine grid\n')
        sol_f = []

        if isParallel:
            sol_f = Parallel(n_jobs=NProc)(
                delayed(globals()[solver_f])(fun_f, trange_time_windows[i, :], y_tot[i], **optionsF) for i in range(n_time_windows))
        else:
            # this could be a parallel loop!
            for i in range(0, n_time_windows):

                f_sol_helper = globals()[solver_f](fun_f, trange_time_windows[i, :], y_tot[i], **optionsF)
                sol_f.append(f_sol_helper)

        if plotEach:
            remove_lines(ax[1])
            regenerate_plot(ax[1], sol_f, n_time_windows, trange_time_windows, colors, output_folder)

        # reconstruct global t and y vectors
        t = np.array([trange[0]])

        y = np.array([y0])
        y = np.transpose(y)

        for i in range(0, len(sol_f)):
            t = np.concatenate((t, sol_f[i].t[1:]))
            y = np.concatenate((y, sol_f[i].y[:, 1:]), axis=1)

        # compare solutions at start and end of each time interval
        for i in range(1, n_time_windows):
            y_fend[i, :] = sol_f[i - 1].y[:, -1]
            y_f0[i, :] = sol_f[i].y[:, 0]

        y_fend[-1, :] = sol_f[-1].y[:, -1]

        err[k - 1] = para_AbsRel_norm(y_fend[0:-1, :], y_f0[0:-1, :], abs_tol, rel_tol)

        print("PARAREAL: ERROR after %d ITERATIONS is %g \n" % (k, err[k - 1]))

        if (k > NIter):
            print("Parareal DID NOT CONVERGE after %d iterations with an error norm of %f\n" % (k, err[k - 1]))

    print("PARAREAL CONVERGED SUCCESSFULLY after %d iterations with an error norm of %f\n" % (k, err[k - 1]))

    ode_sol_obj = scipy.integrate._ivp.ivp.OdeResult(t=t, y=y, err=err)
    return ode_sol_obj


def para_AbsRel_norm(x, y, AbsTol, RelTol):
    # inspired by Hairer, Norsett, Wanner
    # Solving Ordinary Equations 1, page 167
    sc = AbsTol + abs(x) * RelTol
    M = np.divide(abs(x - y), sc)
    twoNorm = np.sqrt(np.sum(np.square(M), axis=1)) / np.sqrt(x.shape[1])
    infNorm = max(twoNorm)  # infinity norm of the 2norm
    return infNorm
