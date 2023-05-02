
/* Thermal Thin Shell Approximation - Copyright (C) 2022 E. Schnaubelt,
   CERN/TU Darmstadt */

/* This file lists the parameters used by GetDP and/or Gmsh. */

// ------------------ Geometrical parameters -----------------------------------
w_Bare = 0.00101;
h_Bare = 0.00161;
w_Ins  = 0.00002;

// ------------------ Solver parameters --------------------------------------

// time stepping
t_init  = 0; // starting time
t_theta = 1; // theta for theta time stepping scheme
T       = 0.1; // period

DefineConstant[
  // Use thin shell approximation?
  Flag_thinShell = {0, Choices{0,1},
    Name "TSA/Use thin shell approx"},

  // # of elements inside thin shell
  N_ele = {2, Min 0, Max 20, Name "TSA/# thin shell elements",
    Visible Flag_thinShell == 1},

  // number of nodes of structured mesh
  numNodesHeight = {5, Min 2, Max 20, Name "Mesh/# nodes along thickness of block"},

  numNodesWidth = {3, Min 2, Max 20, Name "Mesh/# nodes along width of block"},

  numNodesWidthIns = {3, Min 2, Max 20, Name "Mesh/# nodes along insulation",
    Visible Flag_thinShell == 0},

  numLayersHeight = {3, Min 2, Max 20, Name "Mesh/# nodes along height of block"},

  t_end = {0.02, Name "Mesh/# nodes along height of block"}, 
  
  t_step = {0.005, Name "Mesh/# nodes along height of block"}
];

// ------------------ Ids for physical regions ---------------------------------
SHELL       = 10000;
SHELL_MINUS = 20000;
SHELL_PLUS  = 30000;
SHELL_AUX   = 40000;
SHELL_BND   = 50000;

BARE_1     = 10;
BARE_2     = 20;
BND_RIGHT  = 30;
BND_LEFT   = 40;
INSULATION = 50;
