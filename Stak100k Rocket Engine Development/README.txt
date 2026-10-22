================================================================================
  SOLID ROCKET ENGINE DESIGN TOOL  -  EngineDesignV2
  Space Propulsion - Big Assignment  |  Politecnico di Milano
================================================================================

OVERVIEW
--------

The workflow is split into four main scripts (see main.m):
  1. engine_graindesignFinal.m     -- grain + structural design
  2. engine_montecarlo.m           -- Monte Carlo analysis (do not run, it takes a lot of time!)
  2b. engine_resultsMontecarlo     -- displays the result of the report Monte Carlo analysis
  3. engine_thermalCoatingDesign.m -- thermal/ablative design
  4. hydraulics.m                  -- water and pressurant system sizing


HOW TO RUN
----------
1. Open MATLAB and set the working directory to this folder.
2. execute main.m, altogether or by sections
   EXRA NOTE: the scripts called in main can also be executed one by one, 
              this can only be done in the specific order reporterd in main, which ensures the correct code sequentiality 


DEPENDENCIES
------------
  - NASA CEA executable accessible from the MATLAB path (required by src/CEA/)
  - No additional toolboxes required beyond the Optimization Toolbox (used by
    solverCoolantChannelFminCon.m, which is an optional variant)


FOLDER STRUCTURE
----------------
  EngineDesignV2/
  ├── main.m                          top-level run script
  ├── PathsInit.m                     adds all src/ and data/ to MATLAB path
  ├── engine_graindesignFinal.m       Stage 1: grain + nozzle + structure
  ├── engine_montecarlo.m             Stage 2: Monte Carlo run
  ├── engine_resultsMontecarlo.m      Stage 2: MC results showing / post-processing / plots
  ├── engine_thermalCoatingDesign.m   Stage 3: thermal / TBC / ablative
  ├── hydraulics.m                    Stage 4: pressurant circuit sizing
  |                                   In case of single scritps execution 3 and 4 need to be executed in sequence! 
  |                                   This is already taken care of in the main.m corresponding section
  |                                   
  ├── data/                           input data and saved workspaces
  ├── src/
  │   ├── CEA/                        NASA CEA wrappers
  │   ├── general/                    shared utilities
  │   ├── geometry/                   grain mesh and burning-area evolution
  │   ├── RK4/                        time-integration solvers
  │   ├── RAOnozzle/                  Rao nozzle contour optimisation
  │   ├── Montecarlo/                 MC / quasi-MC 
  │   ├── design/                     model comparison utilities
  │   ├── sizing/                     structural and hydraulic sizing
  │   │   ├── pressureVessel/
  │   │   └── idraulics/
  │   ├── complete_cooling/           cooling package 
  │   └── plots/                      visualisation scripts

================================================================================
