# 1. modify these files.
#Grid.h                        
#Grid_RemoveMassFromGrid.C     
#Grid_ComputeCoolingTime.C     
#Grid_GrackleWrapper.C         
#Grid_CorrectForRefinedFluxes.C
#Grid_SolveHydroEquations.C    
#Grid_IdentifySpeciesFields.C  
#typedefs.h                    
#fortran.def                   
#GrackleReadParameters.C       
#GrackleWriteParameters.C      
#Grid_AddFeedbackSphere.C      
#global_data.h                
#CoolData.h                   

# 2. modify initial condition files
### cosmology simulation
#CosmologySimulationInitialize.C         
#Grid_CosmologySimulationInitializeGrid.C
### when PartitionNestedGrids == 1, modify them
#NestedCosmologySimulationInitialize.C          
#Grid_NestedCosmologySimulationInitializeGrid.C 
### collapse test
#CollapseTestInitialize.C
#Grid_CollapseTestInitializeGrid.C

# If using SN dust models, modify
#Star_HitEndpoint.C                
#Star_SetFeedbackFlag.C            
#Star_CalculateFeedbackParameters.C
### add particle attribute "FaintSN" (Under construction)
#Star.h
#Star_AssignFinalMassFromIMF.C
#StarRoutines.C

# (OPTION) use old version of hydro files with fixed gamma
#Grid.h                    
#Grid_SolveHydroEquations.C
### add these files
#Grid_SolvePPM_DE_vg.C        Grid_SolvePPM_DE.C         
#Grid_xEulerSweep_vg.C        Grid_xEulerSweep.C         
#Grid_yEulerSweep_vg.C        Grid_yEulerSweep.C         
#Grid_zEulerSweep_vg.C        Grid_zEulerSweep.C         
#inteuler_vg.F                inteuler.F                 
#flux_hll_vg.F                flux_hll.F                 
#flux_hllc_vg.F               flux_hllc.F                
#euler_vg.F                   euler.F                    
#pgas2d_dual_vg.F             pgas2d_dual.F              
#calcdiss_vg.F                calcdiss.F                 
#twoshock_vg.F                twoshock.F                 
#flux_twoshock_vg.F           flux_twoshock.F            
#intprim_vg.F                 intprim.F
#intpos_vg.F                  intpos.F
#pgas2d_vg.F                  pgas2d.F
### too add them, modify
#Make.config.objects
#euler_sweep.h


###Before you use modification, type
#make machine-{MACHINE FILE}
#make grackle-yes
#make max-baryons-100
#make max-tasks-per-node-128
#make opt-high
