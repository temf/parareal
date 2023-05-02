import os
import subprocess
import numpy as np
import tempfile

import scipy
import time

from getdp_read_preresolution import getdp_read_preresolution

def solve_ivp_getdp(fun, t_span, y0, **getdp_options):
    """
    Solve initial value problem using GetDP.
    :param fun: function to integrate (i.e., GetDP pro_file)
    :param t_span: 2-member sequence with interval of integration (t_start, t_end)
    :param y0: array_like, shape (n,), initial state
    :param getdp_options: options for GetDP
    """

    wall_time_start = time.time()

    pro_file = fun

    if not os.path.isfile(pro_file):
        raise Exception("GetDP pro file not found.")
    
    if "complex_getdp" in getdp_options:
        complex_getdp = getdp_options["complex_getdp"]
    else:
        complex_getdp = False

    if "exe" in getdp_options:
        getdp_exe = getdp_options["exe"]

        if subprocess.run([getdp_exe, '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode:
            raise Exception("GetDP executable cannot be called, please check path.")
    else:
        raise Exception("GetDP executable not provided.")

    if "mesh" in getdp_options:
        mesh_file = getdp_options["mesh"]

        if not os.path.isfile(mesh_file):
            raise Exception("GetDP mesh not found.")
    else:
        raise Exception("Please provide a Gmsh mesh.")

    if "pre_processing" in getdp_options:
        pre_processing = getdp_options["pre_processing"]
    else:
        pre_processing = "#1"
        # print("PRE-RESOLUTION: Using default GetDP pre-processing.\n")

    if "postop" in getdp_options:
        postop = getdp_options["post-postop"]
    else:
        postop = "#1"
        # print("PRE-RESOLUTION: Using default GetDP post-operatin.\n")

    if "verbose" in getdp_options:
        verbose = getdp_options["verbose"]
    else:
        verbose = "0"

    if "nl_iteration" in getdp_options:
        nl_iteration = getdp_options["nl_iteration"]
    else:
        nl_iteration = 1

    init = y0
    if not isinstance(init, np.ndarray):
        raise Exception("Initial value not provided as a numpy array.")

    time_start = t_span[0]
    time_end = t_span[1]

    if "time_step" in getdp_options:
        time_step = getdp_options["time_step"]
    else:
        if "num_time_step" in getdp_options: 
            time_step = (time_end - time_start)/getdp_options["num_time_step"]
        else:
            raise Exception("Neither time step size nor number of time steps provided.")

    with tempfile.TemporaryDirectory() as tmp_dir:

        path, file_and_ext = os.path.split(pro_file)
        filename = os.path.splitext(file_and_ext)[0]

        tmpname = os.path.join(tmp_dir, filename)

        resdir = os.path.join(tmp_dir, "res")
        prefile = os.path.join(tmp_dir, filename + ".pre")
        resfile = os.path.join(tmp_dir, filename + ".res")

        exe_list = [getdp_exe, pro_file,
                      "-pre", pre_processing,
                      "-msh", mesh_file,
                      "-name", tmpname,
                      "-v", verbose]

        if subprocess.run(exe_list).returncode:
            raise Exception("Pre-processing failed")

        # try to read number of DoF
        num_dof = getdp_read_preresolution(prefile)
        if num_dof != init.size:
            raise Exception("Initial value has wrong size")
        
        # create initial data 
        create_resolution(resfile, time_start, init, complex_getdp)

        exe_list = [getdp_exe, pro_file, 
                "-restart",
                "-msh", mesh_file,
                "-name", tmpname,
                "-res", resfile,
                "-v", verbose,
                "-setnumber", "t_end", str(time_end),
                "-setnumber", "t_step", str(time_step),
                "-setnumber", "Flag_nl_iteration", str(nl_iteration)
                ]
        
        if subprocess.run(exe_list).returncode:
            raise Exception("Processing failed")
        
        # print(f"solve_ivp_getdp: reading solution from t = {time_start} s to t = {time_end} s")

        time_res, y_vals_res = read_resolution(resfile, num_dof)

        # return same as scipy would
        ode_sol_obj = scipy.integrate._ivp.ivp.OdeResult(t=time_res, y=y_vals_res)

        wall_time_end = time.time()
        time_needed = wall_time_end - wall_time_start
        print(f"solve_ivp_getdp: done from t = {time_start} s to t = {time_end} s in {time_needed} s")

        return ode_sol_obj


def create_resolution(res_file, time, val_array, complex_getdp):
    """
    Create resolution file. 
    :param res_file: .res file path
    :param time: time value (only one step)
    :param val_array: array_like, shape (n,), current state
    :param complex_getdp: boolean, GetDP compiled with complex Petsc
    """

    if np.isnan(np.min(val_array)):
        raise Exception("At least on value is NaN!")

    # atm: only one dofdata object treated
    header_string = f"$ResFormat /* solve_ivp_getdp copyright E. Schnaubelt, CERN/TU Darmstadt \n" + \
        "1.1 0\n" + \
        "$EndResFormat\n" + \
        "$Solution  /* DofData 0 */ \n" + \
        "0 %.16f 0 0" % (time) # dofdata-number time-value time-imag-value time-step-number
    
    footer_string="$EndSolution"

    with open(res_file, 'w') as f:
        if complex_getdp: 
            np.savetxt(f, np.c_[ val_array, np.zeros(len(val_array)) ], header=header_string, footer=footer_string, comments="")
        else:
            np.savetxt(f, val_array, header=header_string, footer=footer_string, comments="")

def read_resolution(res_file, num_dof): 
    """
    Read resolution file. 
    :param res_file: .res file path
    :param num_dof: number of degrees of freedom
    """

    with open(res_file, 'r') as f:
        fct_vals = []

        t = []

        t_idx = 0 
        old_t_idx = 0

        for line in f:
            if line.find("$Solution") >= 0:
                # dofdata-number time-value time-imag-value time-step-number
                meta_data = f.readline().split()

                dofdata = int(meta_data[0])
                time = float(meta_data[1])
                t_step = float(meta_data[3])

                # print(meta_data)

                if old_t_idx < t_step: # check if time step advances
                    old_t_idx = t_step
                elif old_t_idx > t_step: 
                    raise Exception(f"GetDP runner: error reading file {res_file}. Time step {t_step} is stored after {old_t_idx}.")
                else: 
                    pass
                    # warnings.warn(f"GetDP runner reading file {res_file}. Time step {t_step} == {old_t_idx} already read before.")
                
                t.append(time) 
                fct_vals.append([])

                for _ in range(num_dof):
                    fct_vals[t_idx].append(float(f.readline().split()[0]))
                    
                t_idx = t_idx + 1

            elif line.find("$ResFormat") >= 0:
                res_format = f.readline().split()[0]
                if res_format != "1.1": 
                     raise Exception(f"GetDP runner. unknown res file format version: {res_format}")

        fct_vals = np.array(fct_vals).transpose()
        t = np.array(t)

        if np.isnan(np.min(fct_vals)):
            raise Exception("At least on value is NaN!")
        
        return t, fct_vals