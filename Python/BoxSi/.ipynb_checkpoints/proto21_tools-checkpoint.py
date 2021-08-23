ps  = 1; 
ns  = 1.e3*ps; 
us  = 1.e3*ns; 
ms  = 1.e3*us; 
sec = 1.e3*ms;


    
# ------------------------------------------------------------------------------------- #
def process_SiPMgroups_to_events(groups, NSiPM_min=3, dT=100, fdebug=0):#{
    '''
    input:
    --------
    groups      pandas.DataFrame
    NSiPM_min   minimal number of allowed SiPMs that fired per event
    
    return:
    --------
    events  pandas.DataFrame
    '''    
    import sys;  sys.path.insert(0, '/Users/erezcohen/Desktop/PETsys/Software/PETsysAnalysis/Python'); 
    from PETsys_analysis_tools import *;

    dataset = groups

    event_idx = 0
    N = dataset['N(SiPMs)'][0]
    Qtot, t, channels, detector = 0, [], [], 0
    events = pd.DataFrame() #events = dict()
    for i,SiPMdata in dataset.iterrows():#{
        if (i%(len(dataset)/10)==0): print '%.1f'%(100*float(i)/len(dataset)),'%'

        if fdebug>1: print SiPMdata['N(SiPMs)'],SiPMdata['n(SiPM)'],N,SiPMdata['channel'],SiPMdata['charge']

        # stop iterating if we reached a new event
        if (SiPMdata['N(SiPMs)'] != N) or (i==len(dataset)-1):#{        
            t_min = np.min(t); #t_mean = np.mean(t)
            # determine which detector was hit
            # print (np.min(channels),np.max(channels))
            if (64 <= np.min(channels)) and (np.max(channels)<=127):
                detector = 12;
                detector_type = 'KETEK 3x3 - 1';
            elif (128 <= np.min(channels)) and (np.max(channels)<=192):
                detector = 3;
                detector_type = 'KETEK 3x3 - 2';
            elif (895 <= np.min(channels)) and (np.max(channels)<=1023):
                detector = 6;
                detector_type = 'SensL 6x6 - 1';
            elif (767 <= np.min(channels)) and (np.max(channels)<=895):
                detector = 9;
                detector_type = 'SensL 6x6 - 2';
            else:
                detector = 0; # unknown
                detector_type = 'unknown';

            # if event is valid, plug into a Pandas.DataFrame
            valid_event = False
            if (N>=NSiPM_min) and (Qtot>0): #{
                valid_event = True
            #}
            if valid_event:#{
                events = events.append([pd.DataFrame(data={'event':[event_idx],'N(SiPMs)':[N],
                              'detector':[detector],'det.type':[detector_type],
                              'time[ns]':[t_min/ns],'Qtot':[Qtot]})], ignore_index=True)
                event_idx = event_idx+1
            #}
            # initialize next event-group
            N = SiPMdata['N(SiPMs)']
            Qtot, t, channels, detector = 0, [], [], 0        
        #}

        # iterate over all SiPMs that fired in this event
        Qtot = Qtot + lineariseChargeDeposited(SiPMdata['charge'])
        t.append(SiPMdata['time'])
        channels.append(SiPMdata['channel'])
    #}
    
    KETEK_events = events[(events['detector']==3) | (events['detector']==12)]
    SensL_events = events[(events['detector']==6) | (events['detector']==9)]
    print(len(events),'events during',dT,'seconds, with average rate %.2f+/-%.2f Hz'%(len(events)/dT,sqrt(len(events))/dT))
    print('%.1f'%(100.*len(KETEK_events[KETEK_events.detector==12])/len(events))+'% in detector 12 (KETEK)')
    print('%.1f'%(100.*len(KETEK_events[KETEK_events.detector==3])/len(events))+'% in detector 3 (KETEK)')
    print('%.1f'%(100.*len(SensL_events[SensL_events.detector==6])/len(events))+'% in detector 6 (SensL)')
    print('%.1f'%(100.*len(SensL_events[SensL_events.detector==9])/len(events))+'% in detector 9 (SensL)')
    print('%.1f'%(100.*len(SensL_events[SensL_events.detector==0])/len(events))+'% in unknown detector')

    
    return events,KETEK_events,SensL_events
#}
# ------------------------------------------------------------------------------------- #



# ------------------------------------------------------------------------------------- #
def ADC_2_Edep_approximately(ADC, ScintillatorAreaCoverage=1, ScintillationYield_phPerMeV=8000):
    '''
    calibration for plastic scintillators that emit around 8000 ph/MeV
    ScintillatorAreaCoverage is how much of the scintillator is read by the SiPM array
    e.g. for KETEK 3x3 on a 51x51 scintillator size it is 1/4
    
    from "DrawCoincidenceSpectrum.ipynb" for LYSO which yields 39,900 ph/MeV
    calibration polynomial for channel 482 : Edep [keV] = 27.63 x ADC - 26.68
    calibration polynomial for channel 533 : Edep [keV] = 28.43 x ADC - 9.065
    '''
    import numpy as np
    ScintillationFactor = 39900./ScintillationYield_phPerMeV; # conversion from LYSO to plastic
    Edep_MeV = ((28.*ADC - 15)/ScintillatorAreaCoverage*ScintillationFactor/1000);
    Edep_MeV_err = np.sqrt(np.square(0.5*ADC) + np.square(15))/ScintillatorAreaCoverage*ScintillationFactor/1000;

    return Edep_MeV,Edep_MeV_err
# ------------------------------------------------------------------------------------- #

