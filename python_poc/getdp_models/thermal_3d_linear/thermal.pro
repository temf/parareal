/* Thermal Thin Shell Approximation - Copyright (C) 2022 E. Schnaubelt,
   CERN/TU Darmstadt */

Include "thermal_param.pro";

// ----------------- Strings for naming results --------------------------------
If (Flag_thinShell)
  fileEnding = "_thinShell";
Else
  fileEnding = "_meshedIns";
EndIf

// ----------------- Group definitions -----------------------------------------
Group {
  Bare = Region[{BARE_1, BARE_2}]; // Bare region

  If(!Flag_thinShell)
    Insulation = Region[ INSULATION ]; // Insulation region, NA for TSA
  Else
    Insulation = Region[ {} ];
  EndIf

  Bnd_dirichlet = Region[ BND_RIGHT ]; // Dirichlet boundary on the left

  Bnd_neumann = Region[ BND_LEFT ]; // Neumann boundary on the right

  // shell regions
  Shell= Region[{SHELL}];

  // Auxiliary domains to group regions
  Vol_The = Region[ {Bare, Insulation} ]; // Full thermal domain incl insulation
  Vol_TheAndShell = Region[ {Bare, Insulation} ];

  DomainOneSideOfShell_up = ElementsOf[ {BARE_1, Shell}, OnOneSideOf Shell ];

  DomainOneSideOfShell_down = ElementsOf[ {BARE_2, Shell}, OnOneSideOf Shell ];
}

// ----------------- Function definitions -----------------------------------------
Function {
  // --------- Correction Factors ----------------------------------------------

  // Since the thin shell approach changes the geometry of the problem and by
  // this the thermal mass as well, correction factors can be used. Their
  // influence will vanish as the thin layer becomes smaller.

  If(Flag_thinShell)
    // correct heat capacity
    f_heatCap[] = w_Bare  / (w_Bare + 0.5 * w_Ins);

    // correct thermal conductivity
    f_kappa[] = (w_Bare + 0.5 * w_Ins)/(w_Bare);
  Else
    // no correction needed for meshed insulation
    f_heatCap[] = 1;
    f_kappa[]   = 1;
  EndIf

  // ------------ Material Functions -------------------------------------------

  // thermal conductivity of the bare part w/o correction in W/mK
  kappa_SCSt[] = 0.1;

  // thermal conductivity of the bare part after correction in W/mK
  kappa_Bare[] = f_kappa[] * kappa_SCSt[];

  // thermal conductivity of the insulation in W/mK
  kappa_Insulation[] = 0.0025;

  // heat capacity of bare part w/o correction in J/m^3 K
  heatCap_SCSt[] = 1E6;

  // heat capacity of the bare part after correction in J/m^3 K
  heatCap_Bare[] = f_heatCap[] * heatCap_SCSt[];

  // heat capacity of insulation in J/m^3 K
  heatCap_Insulation[] = 1E10;

  // final piece-wise defined function
  kappa[Bare]       = kappa_Bare[];
  kappa[Insulation] = kappa_Insulation[];

  heatCap[ Bare ]       = heatCap_Bare[];
  heatCap[ Insulation ] = heatCap_Insulation[];


  // ------------ 1D FE Matrices -----------------------------------------------

  // element width of 1D FE element
  delta = w_Ins / (N_ele == 0 ? 1 : N_ele);

  kappa_stiffness[] = (($1 == $2)? 1 : -1 ) * kappa_Insulation[]/delta;

  kappa_mass[] = (($1 == $2)? 2: 1) * delta * kappa_Insulation[]/6;

  heatCap_mass[] = (($1 == $2)? 2: 1) * delta * heatCap_Insulation[]/6;

  // ------------ Excitation ---------------------------------------------------
  impHeatFlux[] = 5E6;
}

// ----------------- Initial and boundary conditions ---------------------------
Constraint {
  // Dirichlet and initial conditions for non-thin shell temperature
  { Name Dirichlet_the ;
    Case {
      { Region Bnd_dirichlet; Value 4.2; Type Assign; }
      { Region Vol_The ;      Value 4.2 ; Type Init; }
    }
  }

  { Name init_the ;
    Case {
      { Region Vol_The ; Value 4.2 ; Type Init; }
    }
  }

  // initialize thin shell unknowns
  { Name Init_shell ;
    Case {
      { Region Shell ; Value 4.2;  Type Init;}
    }
  }
}

// ----------------- Gaussian integration --------------------------------------
Integration {
  { Name Int ;
    Case {
      { Type Gauss ;
        Case {
          { GeoElement Quadrangle ; NumberOfPoints  7 ; }
          { GeoElement Hexahedron ; NumberOfPoints  6 ; }
        }
      }
    }
  }
}

// ----------------- Jacobian --------------------------------------------------
Jacobian {
  { Name Vol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }
  { Name Sur ;
    Case {
      { Region All ; Jacobian Sur ; }
    }
  }
}


// ----------------- Function space definition ---------------------------------
FunctionSpace {

  // temperature for meshed insulation
  { Name Hgrad_T_meshedIns; Type Form0 ;
    BasisFunction {
      { Name sn ; NameOfCoef Tn ; Function BF_Node ;
        Support Region[{Vol_The, Bnd_neumann}] ; Entity NodesOf[All] ; }
    }

    Constraint {
      { NameOfCoef Tn ; EntityType NodesOf ; NameOfConstraint Dirichlet_the ; }
    }
  }

  // temperature for TSA
  { Name Hgrad_T_thinShell; Type Form0 ;
    BasisFunction {
        // contributions outside the thin shell
        { Name sn ; NameOfCoef Tn ; Function BF_Node ;
            Support Region[{Vol_The, Bnd_neumann}] ;
            Entity NodesOf[All, Not Shell]; }

        { Name sdn_up ; NameOfCoef Tdn_up ; Function BF_Node ;
          Support DomainOneSideOfShell_up ; Entity NodesOf[ Shell ] ; }

        { Name sdn_down ; NameOfCoef Tdn_down ; Function BF_Node ;
          Support DomainOneSideOfShell_down ; Entity NodesOf[ Shell ] ; }
    }

    // Subspace to easily refer to plus/minus side of shell
    SubSpace {
      { Name Shell_down ; NameOfBasisFunction { sdn_down } ; }
      { Name Shell_up ;   NameOfBasisFunction { sdn_up } ; }
    }

    // Dirichlet and link contraints
    Constraint {
        { NameOfCoef Tn ; EntityType NodesOf ;
          NameOfConstraint Dirichlet_the ; }

        { NameOfCoef Tdn_up ; EntityType NodesOf ;
          NameOfConstraint init_the ; }

        { NameOfCoef Tdn_down ; EntityType NodesOf ;
          NameOfConstraint init_the ; }
    }
  }

  For i In {1:N_ele-1}
    { Name Hgrad_T_thinShell~{i} ; Type Form0 ;
      BasisFunction {
        { Name sn ; NameOfCoef Tn~{i} ; Function BF_Node ;
        Support Shell ; Entity NodesOf[ All ] ; }
      }

      Constraint {
        { NameOfCoef Tn~{i} ; EntityType NodesOf ;
          NameOfConstraint Init_shell ; }
      }
    }
  EndFor
}

// ----------------- Formulation -----------------------------------------------
Formulation {
  If(!Flag_thinShell) // explicitly meshed
    { Name Thermal_T ; Type FemEquation ;
      Quantity {
          { Name T ; Type Local ; NameOfSpace Hgrad_T_meshedIns ; }
      }
      Equation {
        Integral { [ kappa[] * Dof{d T} , {d T} ] ;
          In Vol_The; Integration Int ; Jacobian Vol ; }

        Integral { DtDof[ heatCap[] * Dof{T}, {T} ];
          In Vol_The; Integration Int; Jacobian Vol;  }

         Integral { [- impHeatFlux[] * ($Time < 0.7? Jn[3, 100 * $Time]: 1/100 * Exp[1/10000 * ($Time - 0.7)])  , {T} ] ;
          In Bnd_neumann; Integration Int ; Jacobian Sur ; }
      }
    } // end explicitly meshed

  Else // Link constraint ansatz

    { Name Thermal_T ; Type FemEquation ;
      Quantity {
        { Name T ; Type Local ; NameOfSpace Hgrad_T_thinShell ; }
        // TSA contributions
        { Name Ti~{0}; Type Local; NameOfSpace Hgrad_T_thinShell[Shell_down]; }
        For i In {1:N_ele-1}
          { Name Ti~{i} ; Type Local ; NameOfSpace Hgrad_T_thinShell~{i}; }
        EndFor
        { Name Ti~{N_ele}; Type Local;
          NameOfSpace Hgrad_T_thinShell[Shell_up]; }
      }
      Equation {
        Integral { [ kappa[] * Dof{d T} , {d T} ] ;
          In Vol_The; Integration Int ; Jacobian Vol ; }

        Integral { DtDof[ heatCap[] * Dof{T}, {T} ]; In Vol_The;
          Integration Int; Jacobian Vol;  }

        Integral { [- impHeatFlux[] , {T} ];
         In Bnd_neumann; Integration Int; Jacobian Sur;}

        For i In {0:N_ele-1} // loop over 1D FE elements
          For k In {1:2} // row of 1D FE matrix
            For l In {1:2} // column of 1D FE matrix
              Integral { [ kappa_mass[k, l] *
                Dof{d Ti~{i + k - 1}} , {d Ti~{i + l - 1}} ];
                In Shell; Integration Int; Jacobian Sur; }

              Integral { DtDof[ heatCap_mass[k, l] *
                Dof{Ti~{i + k - 1}} , {Ti~{i + l - 1}} ];
                In Shell; Integration Int; Jacobian Sur; }

              Integral { [ kappa_stiffness[k, l] *
                Dof{Ti~{i + k - 1}} , {Ti~{i + l - 1}} ];
                In Shell; Integration Int; Jacobian Sur;}
            EndFor
          EndFor
      EndFor
      }
    }
    EndIf // end Link constraint ansatz
}

// ----------------- Resolution ------------------------------------------------
Resolution {
  { Name Thermal_T ;
    System {
      { Name Sys_The ; NameOfFormulation Thermal_T ; }
    }
    Operation {

      // PostOperation[resetMaxTemp]; // clear maximum temperature output file

      InitSolution Sys_The; // init. the solution using init. constraints

      // PostOperation[PrintMaxTemp]; // get maximum temperature (= init. temp)

      TimeLoopTheta [t_init, t_end, t_step, t_theta] {
          // Generate and solve system, get maximum temperature
          Generate Sys_The ; Solve Sys_The;
			SaveSolution Sys_The ;

          // PostOperation[PrintMaxTemp];
      }

      // SaveSolution Sys_The ;

      // create temperature field maps
      // PostOperation[Map_T];
    }
  }
}

// ----------------- PostProcessing --------------------------------------------

PostProcessing {

  { Name DoNothing; }

  { Name Thermal_T ; NameOfFormulation Thermal_T ;
    PostQuantity {

      // Temperature map
      { Name T ; Value {
          Term { [ {T} ] ; In Vol_The ; Jacobian Vol ; }
        }
      }

      // Heat flux
      { Name Flux ; Value {
          Term { [kappa[] *  {d T} ] ; In Vol_The ; Jacobian Vol ; }
        }
      }

     // Maximum temperature in bare part from register 1 (saved in
     // post-operation by StoreMaxInRegister)
     { Name Tmax; Value{  Term{ Type Global; [#1]; In Bare; Jacobian Vol;} } }
    }
  }
}

// ----------------- PostOperation ---------------------------------------------
PostOperation PrintMaxTemp UsingPost Thermal_T {
  // Get maximum in bare region and store in register 1
  Print[ T, OnElementsOf Bare, StoreMaxInRegister 1, Format Table ,
    LastTimeStepOnly 1, File "resQOI/dummy.txt"] ;

  // Save maximum temperature in file by appending to existing file
  Print[Tmax, OnRegion Bare, Format TimeTable,
    File StrCat["resQOI/Tmax", fileEnding, ".txt"],
    LastTimeStepOnly 1, AppendToExistingFile 1];
}

// Clean the file saving maximum temperatures
PostOperation resetMaxTemp UsingPost Thermal_T {
  Echo["", Format Table, File StrCat["resQOI/Tmax", fileEnding, ".txt"] ];
}

// Plot field maps
PostOperation Map_T UsingPost Thermal_T {
  Print[ T, OnElementsOf Vol_The,
    File StrCat["resField/T", fileEnding, ".pos" ] ] ;

  Print[ Flux, OnElementsOf Vol_The,
    File StrCat["resField/Flux", fileEnding, ".pos" ] ] ;
}
