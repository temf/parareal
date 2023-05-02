/* Thermal Thin Shell Approximation - Copyright (C) 2022 E. Schnaubelt,
   CERN/TU Darmstadt */

// geometry definitions for case 1

// the mesh is created by the Mesh command inside the geo file
Solver.AutoMesh = -1;

Include "thermal_param.pro"; // parameters needed by Gmsh and GetDP

/* Thermal Thin Shell Approximation - Copyright (C) 2022 E. Schnaubelt,
   CERN/TU Darmstadt */

// ----------------------- POINTS ----------------------------------------------

Point(1) = {0, 0, 0};
Point(2) = {w_Bare + Flag_thinShell * w_Ins/2, 0, 0};
Point(3) = {w_Ins + w_Bare, 0, 0};
Point(4) = {w_Bare + w_Ins + w_Bare, 0, 0};
Point(5) = {0, h_Bare, 0};
Point(6) = {w_Bare + Flag_thinShell * w_Ins/2, h_Bare, 0};
Point(7) = {w_Ins + w_Bare, h_Bare, 0};
Point(8) = {w_Bare + w_Ins + w_Bare, h_Bare, 0};

// --------------------- HORIZONTAL LINES --------------------------------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {5, 6};
Line(5) = {6, 7};
Line(6) = {7, 8};

// ---------------------- VERTICAL LINES ---------------------------------------
Line(10) = {1, 5};
Line(11) = {2, 6};
Line(12) = {3, 7};
Line(13) = {4, 8};

// --------------- SURFACES ----------------------------------------------------
Curve Loop(1) = {1, 11, -4, -10};
Curve Loop(2) = {2, 12, -5, -11};
Curve Loop(3) = {3, 13, -6, -12};

Surface(1) = {1};
Surface(2) = {2};
Surface(3) = {3};

// -------------- EXTRUSION ----------------------------------------------------

num[] = Extrude {0,0, 1E-3} { Surface{1, 2, 3}; Layers{(numLayersHeight - 1)}; 
  Recombine; };


// -------------- MESHING PARAMS -----------------------------------------------
Transfinite Curve{10, 11, 12, 13} = numNodesHeight;
Transfinite Curve{1, 3, 4, 6}     = numNodesWidth;
Transfinite Curve{2, 5}           = numNodesWidthIns;

Transfinite Surface{1, 2, 3};
Recombine Surface "*";

// -------------- PHYSICAL REGIONS ---------------------------------------------
Physical Surface("Bnd right", BND_RIGHT) = {70};
Physical Surface("Bnd left", BND_LEFT)   = {34};

Physical Volume("Bare 1", BARE_1)         = {1};
Physical Volume("Bare 2", BARE_2)         = {3};
Physical Volume("Insulation", INSULATION) = {2};

// --------------- MESH --------------------------------------------------------
Mesh 3;
Save "thermal.msh";
