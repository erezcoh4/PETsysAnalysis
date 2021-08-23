import numpy as np, pandas as pd, matplotlib.pyplot as plt
import sys; sys.path.insert(0, '/'); 
from my_tools import *; from plot_tools import *
from datetime import date

main_data_path = '/Users/erezcohen/Desktop/PETsys/data/'

import mpl_scatter_density # adds projection='scatter_density'
from matplotlib.colors import LinearSegmentedColormap
# colormap for kde scatter plots
white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#ffffff'),
    (1e-20, 'lightblue'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#fde624'),
    (1, 'red'),
], N=256)




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



def calibrate_charge_2_keV_using_Na22_data(bins=np.linspace(0,58,150), do_plot=True,                                          
                                           data_path='/Users/erezcohen/Desktop/data/PETsys/Prototype1/Na22_source_0.8uC_3Dec2020/',                                          
                                           Na22fname='Na22_1000sec_single.dat', dT_Na22 = 1000.0,
                                           Bkgfname='Background_100sec_single.dat', dT_Bkg = 100.0 
                                          ):
    '''
    calibrate arbitrary ADC to energy deposition keV
    using Na22 data
    Jan-19, 2021
    
    return:
    Epoly - calibration polynomial per channel
    x_keV - calibrated bin centers (in keV)
    axes - array of two axes
    '''    
    Na22_singles = pd.read_csv(data_path+Na22fname,delimiter='\t',names=['time','charge','channel'])
    Bkg_singles = pd.read_csv(data_path+Bkgfname,delimiter='\t',names=['time','charge','channel'])
    
    from scipy.signal import find_peaks

    channels = np.unique(Na22_singles.channel)
    print ('detected events in the following channels: ',channels)
    Na22singles,Bkgsingles = dict(),dict()
    singlesPerChnnl,axes = dict(),dict()
    hNa22SigBkg,hBkg,hNa22Sig,hNa22Sig_err = dict(),dict(),dict(),dict()    
    x = bins[1:]
    # calibrate charge to energy
    peaks_Na22sing = dict()
    Epoly,x_keV = dict(),dict()

    for ch,chIdx in zip(channels,[1,2]): #{   
        # Na22 signal+background
        Na22singles[ch] = Na22_singles[Na22_singles.channel==ch]    
        hNa22SigBkg[ch],edges = np.histogram(lineariseChargeDeposited(Na22singles[ch].charge),bins=bins);
        # only background
        Bkgsingles[ch] = Bkg_singles[Bkg_singles.channel==ch]    
        hBkg[ch],edges = np.histogram(lineariseChargeDeposited(Bkgsingles[ch].charge),bins=bins);
        # Na22 signal
        hNa22Sig[ch] = np.array(hNa22SigBkg[ch] - hBkg[ch]*(dT_Na22/dT_Bkg))
        hNa22Sig_err[ch] = np.sqrt( np.array(hNa22SigBkg[ch] + hBkg[ch]*(dT_Na22/dT_Bkg)))     

        # calibrate using the Na22 source
        peaks_Na22sing[ch],tmp = find_peaks(hNa22Sig[ch][1:-1],prominence=(0.05*np.max(hNa22Sig[ch]), None),distance=30)
        print ('in channel ',ch,',peaks at bins number',peaks_Na22sing[ch]    )
        # do not use backscattering at 255.5 keV, 
        # use photo-peak at 511 keV, de-excitation gamma at 1274.5
        peaks2use = peaks_Na22sing[ch][1:];
        Q_au = np.array(x[peaks2use])
        E_keV = np.array([511,1274.5])
        print ('peaks2use:',peaks2use,'Q_au:',Q_au)
        popt = np.polyfit( Q_au, E_keV , 1 )
        Epoly[ch] = np.poly1d(popt)
        x_keV[ch] = Epoly[ch](x)
    #}
    
    # print resulting calibration polynomial parameters
    print ('calibration polynomial for channel',channels[0],':',Epoly[channels[0]])
    print ('calibration polynomial for channel',channels[1],':',Epoly[channels[1]])
    
    if do_plot:#{
        # apply calibration to Na22
        for ch,chIdx in zip([channels[0]],[1,2]):    
            # signal+background
            Na22singles[ch]['Energy [keV]'] = Epoly[channels[0]](Na22singles[ch].charge)
        fig=plt.figure(figsize=(14,6))
        for ch,chIdx in zip(channels,[1,2]):    
            ax=fig.add_subplot(1,2,chIdx)
            plt.errorbar( x_keV[ch], hNa22Sig[ch],  yerr=hNa22Sig_err[ch] ,
                         linewidth=0.4,marker='s',markersize=5,
                         linestyle='-.',color='forestgreen',markeredgecolor='black',
                         capthick=1,capsize=3 );    
            set_axes(ax=ax,x_label='$E_{dep}$ [keV]',y_label='number of events',title='$^{22}$Na',                     
                     fontsize=17,do_add_grid=True,do_add_legend=False,legend_loc='best');
            axes[ch] = ax;

        plt.tight_layout()   
    #}    
    return Epoly,x_keV,axes




def fit_gaussian_to_data(x=[],do_add_plot=False,ax=None,xlims=None,bins=None):
    '''
    fit gaussian to data
    '''
    if xlims is not None:
        x2fit = x[(x>np.min(xlims))&(x<np.max(xlims))];
    else:
        x2fit = x

    print (len(x2fit),'events fitting to Gaussian')
    (mu, sigma) = norm.fit(x2fit)        
    print ('mu:',mu,',sigma:', sigma)
    # uncertainty, conservative 5/20% - need to revise this!
    if len(x2fit)>1000:  sigma_err = 0.05*sigma
    else: sigma_err = 0.2*sigma
        
    if do_add_plot:
        if bins is None:
            bins=np.linspace(np.percentile(x2fit,1)-4*sigma,np.percentile(x2fit,99)+4*sigma,50)            
        x_fit = np.linspace(np.min(bins),np.max(bins),200)
        y_fit = mlab.normpdf( x_fit, mu, sigma)
        histo,histbins,patches = plt.hist( x, bins=bins,  normed=1, edgecolor='black',label='data, %d events'%len(x2fit));
        plt.plot(x_fit, y_fit, '--r', linewidth=2,
                 label='gaussian fit, $\sigma=%.1f\\pm%.1f$'%(sigma,sigma_err))




def time_difference_in_2SiPM_coincidence_measurement(coinc_data=None,Epoly=None,
                                                     channels = [482,533], # channel numbers in data
                                                     calibration_channels=[482,533], # channel numbers in calibration polynomial
                                                     channelIdcs = [2,1], # channel indices order in coincidence data
                                                     bins=np.linspace(0,750,50),
                                                     Edeplims_keV=None,dt_lims=(-2e3,2e3),
                                                     dt_bins=None, detector_labels=None):
    '''
    Epoly - calibration polynomial from ADC [a.u.] to energy [keV], 
            has two keys which are the two channel numbers, ordered by channleIdcs in coinc_data...
            
    last update Dec-9, 2020
    '''
    if channels is None:#{
        channels = np.unique([coinc_data.ch1,coinc_data.ch2])
        print ('detected events in the following channels: ',channels)
    #}
    if calibration_channels is None:#{
        calibration_channels = channels;
    #}
    if detector_labels is None:#{
        detector_labels = ['Left','Right']
    #}
    # omit events with unreasonable charge deposition
    coinc_data = coinc_data[(coinc_data['Q1']>0)&(coinc_data['Q2']>0)]
    
    if Edeplims_keV is None:
        coinc_2_use = coinc_data
    else:
        if len(Edeplims_keV)==2:
            min1=np.min(Edeplims_keV); min2=min1; 
            max1=np.max(Edeplims_keV); max2=max1;
        elif len(Edeplims_keV)==4:
            min1 = Edeplims_keV[0];
            max1 = Edeplims_keV[1];
            min2 = Edeplims_keV[2];
            max2 = Edeplims_keV[3];

        coinc_2_use = coinc_data[ (min1
                                   < Epoly[calibration_channels[0]](lineariseChargeDeposited(coinc_data['Q%d'%channelIdcs[0]])))
                                & (Epoly[calibration_channels[0]](lineariseChargeDeposited(coinc_data['Q%d'%channelIdcs[0]])) 
                                   < max1)
                                & (min2
                                   < Epoly[calibration_channels[1]](lineariseChargeDeposited(coinc_data['Q%d'%channelIdcs[1]])))
                                & (Epoly[calibration_channels[1]](lineariseChargeDeposited(coinc_data['Q%d'%channelIdcs[1]])) 
                                   < max2)]        

    print (len(coinc_2_use),'events in coincidence window and Edep limits')
    delta_t = coinc_2_use.t1 - coinc_2_use.t2
    hEdep,hEdep_err,Edep = dict(),dict(),dict()

    
    fig=plt.figure(figsize=(16,12));
    for i,ch,calib_ch,chIdx,subplotTitle in zip([1,2],channels,calibration_channels,channelIdcs,detector_labels):

        Qlin = lineariseChargeDeposited(coinc_data['Q%d'%chIdx])
        Edep[ch] = Epoly[calib_ch](Qlin)
        hEdep[ch],edges = np.histogram(Edep[ch],bins=bins);        
        hEdep_err[ch] = np.sqrt(hEdep[ch])

        ax=fig.add_subplot(2,2,chIdx)
        plt.errorbar( edges[:-1], hEdep[ch],  yerr=hEdep_err[ch] ,
                     linewidth=0.4,marker='s',markersize=5,
                     linestyle='-.',color='royalblue',markeredgecolor='black',
                     capthick=1,capsize=3 );

        if Edeplims_keV is not None:
            if i==1: plt.plot([min1,min1],[0,np.max(hEdep[ch])],'--r', [max1,max1],[0,np.max(hEdep[ch])],'--r')
            else: plt.plot([min2,min2],[0,np.max(hEdep[ch])],'--r', [max2,max2],[0,np.max(hEdep[ch])],'--r')
        set_axes(ax=ax,x_label='$E_{dep}$ [keV]',y_label='number of coincidence events',
                 xlim=(np.min(bins),np.max(bins)),title=subplotTitle,
                 fontsize=17,do_add_grid=True,do_add_legend=False,legend_loc='best');

    # kde scatter plot
    ax=fig.add_subplot(2,2,3, projection='scatter_density')
    density = ax.scatter_density(np.array(Edep[channels[0]]),np.array(Edep[channels[1]]),  
                       cmap=white_viridis, dpi=8, downres_factor=20)
    fig.colorbar(density, label='')
    if Edeplims_keV is not None:
        plt.plot([min1,min1],[min2,max2],'--r',
                 [max1,max1],[min2,max2],'--r',
                 [min1,max1],[min2,min2],'--r',
                 [min1,max1],[max2,max2],'--r')

    # cosmetics
    set_axes(ax=ax,x_label='$E_{dep}^{R}$ [keV]',xlim=(np.min(bins),np.max(bins)),
             y_label='$E_{dep}^{L}$ [keV]',ylim=(np.min(bins),np.max(bins)),
             fontsize=17,do_add_grid=True,do_add_legend=False,legend_loc='best');     
        
        
    ax=fig.add_subplot(2,2,4)
    fit_gaussian_to_data(x=delta_t,xlims=dt_lims,do_add_plot=True,ax=ax,bins=dt_bins)
    # cosmetics
    set_axes(ax=ax,x_label='$\Delta t$ [ps]',
             y_label='frequency',
             title='$\Delta t$ in coincidence events',
             fontsize=17,do_add_grid=True,do_add_legend=True,legend_loc='best');     
    plt.tight_layout()       
    
    
    
    
    
    
def get_all_differences_R_L(t_L,t_R):#{
    # compute time differences between all combinations of t_L and t_R
    mesh = np.array(np.meshgrid(t_R,t_L))
    combinations = mesh.T.reshape(-1, 2)
    dt = np.diff(combinations)
    return dt
#}



# source activity at time of measurement
def Na22_activity_per_day( d_t = date(2021,1,9) ):
    from datetime import date

    d_t0 = date(2021,1,6)
    Na22_activity_t0 = 10 # uC
    Na22_t_half_days = 2.6 * 365 
    tau_Na22_days    = Na22_t_half_days/np.log(2)
    dt_days = abs((d_t - d_t0).days)
    
    Na22_activity_t = Na22_activity_t0 * exp(- dt_days/tau_Na22_days ) # uC
    # also convert from uC to Bqrl
    uC2Bqrl = (1.0e-6)/(2.703e-11)
    Na22_eventRateHz_t = Na22_activity_t * uC2Bqrl # events/sec
    return Na22_activity_t, Na22_eventRateHz_t


# Cf252 source activity at time of measurement
def Cf252_activity_per_day( d_t = date(2021,1,9) ):
    # return
    # -------
    # Cf252_fissionRateHz_t, Cf252_activity_t, Cf252_eventRateHz_t
    from datetime import date

    d_t0 = date(2021,7,1)
    Na22_activity_t0 = 5 # uC
    Na22_t_half_days = 2.645 * 365 
    tau_Na22_days    = Na22_t_half_days/np.log(2)
    dt_days = abs((d_t - d_t0).days)
    
    Cf252_activity_t = Na22_activity_t0 * exp(- dt_days/tau_Na22_days ) # uC
    # also convert from uC to Bqrl
    uC2Bqrl = (1.0e-6)/(2.703e-11)
    Cf252_eventRateHz_t = Cf252_activity_t * uC2Bqrl # events/sec
    # only 3.1% of Cf252 disintegrations are spontaneous fission
    # [https://www.sciencedirect.com/topics/medicine-and-dentistry/californium-252]
    Cf252_fissionRateHz_t = Cf252_eventRateHz_t * 0.031  
    return Cf252_fissionRateHz_t, Cf252_activity_t, Cf252_eventRateHz_t



# rate of gamma rays emitted from the source per second
def Na22_activity_Bqrl_to_gamma_rate_Hz( Na22_activity_Bqrl=1 ):
    '''
    [https://ehs.umich.edu/wp-content/uploads/2016/04/Sodium-22.pdf]
    GAMMA ENERGIES
    - 511.0 keV (179.8% abundance/annihilation)
    - 1274.5 keV (99.9%)
    '''
    return Na22_activity_Bqrl*(179.8/100. + 99.9/100.)
    