import numpy as np, pandas as pd, matplotlib.pyplot as plt
import sys; sys.path.insert(0, '/Users/erezcohen/larlite/UserDev/mySoftware/MySoftwarePackage/mac/'); 
from my_tools import *; from plot_tools import *

main_data_path = '/Users/erezcohen/Desktop/PETsys/data/'



'''
coincidence FoM labels
------------------------
[page 37 of software guide]

    mh n1, the number of events that belong to the same group of this particular event (for the event in the trigger region of higher ID number). If mh n is 1, then no other events were detected within a radius of 100mm and a 100 ns time window (default grouping spacial and time windows).

    mh j1, the number ID (from 0 to mh n-1) of this event in the group it belongs to (ordered from higher to lower energy, for the event in the trigger region of higher ID number). For example, an event with mh j set to 1 means this 
event has the second highest energy in the group it belongs to.

    Time of detection in picoseconds (for the event in the trigger region of higher ID number).

    Energy of the pulse, for the event in the trigger region of higher ID number (in QDC mode, in arbitrary units of charge; in TOT mode, in nanoseconds). If the energy calibration file (see section 3.1.4) exists and is referenced in the reference configuration file, the energy value will also have arbitrary units of charge.

    Absolute channel ID, for the event in the trigger region of higher ID number (see appendix A).

    mh n2, same as mh n1, but for the coincidence event in the trigger region of lower ID number.

    mh j2, same as mh n2, but for the coincidence event in the trigger region of lower ID number.

    Time of detection in picoseconds (for the event in the trigger region of lower ID number).

    Energy of the pulse, for the event in the trigger region of lower ID number (in QDC mode, in ADC units; in TOT mode, in nanoseconds).

    Absolute channel ID, for the event in the trigger region of lower ID number
'''


def lineariseChargeDeposited(Q):
    '''
    # Charge conversion to energy
    [PETsys TOFPET2 ASIC Evaluation Kit - Software User Guide v2019.09, p. 25]
    ...PETsys Electronics has certified that using 
    $$ E=P_0 \cdot P_1^{Q^{P2}} + P_3 \cdot Q - P_0$$
    with 
    $$ (P_0, P_1, P_2, P_3) = (8.00000, 1.04676, 1.02734, 0.31909) $$
    allows to approximate the non-linear response of the ASIC's QDC for all channels,
    up to a level of 1-2% on the energy resolution of the 511 keV photopeak.
    '''
    P0, P1, P2, P3 = 8.00000, 1.04676, 1.02734, 0.31909
    linearisedQ = P0*(np.power(P1, np.power(Q,P2))) + P3*Q - P0
    return linearisedQ