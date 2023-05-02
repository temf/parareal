import subprocess
import tempfile
import os
from getdp_read_preresolution import getdp_read_preresolution


def getdp_create_preresolution(fun, getdp_options):
    """
    Call GetDP to create a .pre file.
    :param fun: function to call GetDP with (i.e., the .pro file)
    :param getdp_options: dictionary of options for GetDP
    """
    pro_file = fun

    if not os.path.isfile(pro_file):
        raise Exception("GetDP pro file not found.")

    if "exe" in getdp_options:
        getdp_exe = getdp_options["exe"]

        if subprocess.run([getdp_exe, '--version']).returncode:
            raise Exception(f"GetDP executable at {getdp_exe} cannot be called, please check path.")
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

    if "verbose" in getdp_options:
        verbose = getdp_options["verbose"]
    else:
        verbose = "0"

    # TODO: optional arguments for setnumber or setstring

    with tempfile.TemporaryDirectory() as tmp_dir:

        path, file_and_ext = os.path.split(pro_file)
        filename = os.path.splitext(file_and_ext)[0]

        tmpname = os.path.join(tmp_dir, filename)

        resdir = os.path.join(tmp_dir, "res")
        prefile = os.path.join(tmp_dir, filename + ".pre")
        resfile = os.path.join(tmp_dir, filename + ".res")

        exe_string = [getdp_exe, pro_file, 
                     "-pre", pre_processing,
                     "-msh", mesh_file,
                     "-name", tmpname,
                     "-v", verbose]

        if subprocess.run(exe_string).returncode:
            raise Exception("Pre-processing failed")

        # try to read number of DoF
        num_dof = getdp_read_preresolution(prefile)
        
        return num_dof
