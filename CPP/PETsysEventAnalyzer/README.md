
 general PETsys group analyser / event collector
 
 last edit Feb-9, 2022 (EOC)
 

 # Execute:
 --------------------------
 
 ./build/apps/PETsysEventAnalyzer <setup> <subdir> <DataLabel> <dT> <DataType> <verbose> <Nrows>
 
 
 # Supported options:
 --------------------------
 
    setup                   BoxSi_proto2.1
                            BoxSi_proto2.2
                            PlasticScintillator_SensL8x8
                            VariableThresholdMeasurements
                            
    DataType                single
                            group
                            events
                            single_PSD
                            group_PSD
                            
    Nrows                   defaul is all events ( Nrows=-1 )

 
 
 
  # Execution examples:
 --------------------------

    *MC:*
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.2 Geant4Sims/YuvalSimulations G4sim_Cf252_50.0mCi 1_sec events

    *data - from PETsys:*
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.2 vth_1PE/ToF_measurements background 10_sec single
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.2 vth_1PE/ToF_measurements nToF_Cf252 10_sec single
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.2 vth_1PE/ToF_measurements nToF_Cf252 5000_sec single 1 100000
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.2 vth_1PE/Cf252_data Cf252 100_sec single
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.2 vth_1PE/KETEK_QDC_calibration_using_gamma_sources only_KETEK_vth12e_5_Bkg 100sec single


    *data - from PETsys with dual readout splitter for PSD:*
    ./build/apps/PETsysEventAnalyzer PSD_1.3 Cf252 Cf252_AMCRYS_UPS113NG 10_sec single_PSD 5


    *data - only events:*
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.2 vth_1PE/Cf252_data Cf252 100_sec events
    
    
    *other data examples - from PETsys:*
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.2 VariableThresholdMeasurements noSource_Vth4 10sec single
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.2 vth_1PE/Cf252_data Cf252  10_sec single
        
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.1 vth_2PE vth12e_10_Na22 10sec single
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.1 Na22_10uC 10_sec single 3
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.1 Na22_10uC 10_sec single 1 100
    ./build/apps/PETsysEventAnalyzer BoxSi_proto2.1 Cf252 10sec group 3
    ./build/apps/PETsysEventAnalyzer PlasticScintillator_SensL8x8 PlasticScintillator_Cf252_1mm_vth1_20_vth2_20 30_sec 3
    ./build/apps/PETsysEventAnalyzer PlasticScintillator_SensL8x8 PlasticScintillator_cosmic_vth1_20_vth2_20 30_sec 3


