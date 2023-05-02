def getdp_read_preresolution(preres_file):
    """
    Read pre-resolution of GetDP. At the moment, it only reads the number of Dof.
    :param preres_file: .pre file to parse
    :return: integer number of degrees of freedom
    """
    with open(preres_file, 'r') as pre_file:

        for line in pre_file:
            if line.find("$DofData") >= 0:
                # Ignore the first four lines after DofData as data not needed atm
                pre_file.readline()
                pre_file.readline()
                pre_file.readline()
                pre_file.readline()

                # parse number of dof by splitting at whitespace
                dof_line = pre_file.readline().split()

                num_dof = int(dof_line[1])

                return num_dof


