#ifndef SiPMChannelPair_HPP
#define SiPMChannelPair_HPP

#include <iostream>

#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <string>
#include <string.h>
#include <fstream>
#include <sstream>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <numeric>
#include <functional>

#include <SiPMHit.hpp>
#include <detectorEvent.hpp>


class SiPMChannelPair
{
public:
    
    
    
    // time in ms in this class
    SiPMChannelPair (){fN = 0;};
    SiPMChannelPair (int eventID,
                     std::vector<int>    channels,
                     std::vector<double> t_ms,
                     std::vector<double> Q,
                     int N,
                     std::string setup="unknown",
                     int verbose=0);
    

    void          SetVerbose (int v)                { verbose = v;};
    
    void                     Print ();
    
    // Getters
    int                      GetNhits ()                       {return fN;};
    int                    GetEventID ()                       {return feventID;};
    int                   GetDetector ()                       {return fdetector;};
//    double             GetEventCharge ()                       {return fQtot;};
    bool                   IsGoodPair ()                       {return fIsGoodPair;};
    double               GetEventTime ()                       {return ftime_ms;}; // return time in ms
    
    int                 GetHitChannel (int hitIdx)             {return fchannels.at(hitIdx);};
    std::vector<int>      GetChannels ()                       {return fchannels;};
    
    double                 GetHitTime (int hitIdx)             {return ft_ms.at(hitIdx);};
    std::vector<double>   GetHitsTime ()                       {return ft_ms;};
    
    double               GetHitCharge (int hitIdx)             {return fQ.at(hitIdx);};
    std::vector<double> GetHitsCharge ()                       {return fQ;};
    

    std::vector<SiPMHit>      GetHits ()                       {return fHits;};
    double                    GetSiPM ()                       {return fSiPM;};
    double                 GetQprompt ()                       {return fQprompt;};
    double                  GetQtotal ()                       {return fQtotal;};
    double             GetTailToTotal ()                       {return fTailToTotal;};
    

    // Setters
    void                  SetSetup (std::string setup)      {fsetup = setup;};
    void                  SetNhits ( int N )                {fN = N;};
    void                SetEventID ( int eventID )          {feventID = eventID;};
    void               SetDetector ( int detector )         {fdetector = detector;};
    // define the event time as the time of the first SiPM hit in the event
    void              SetEventTime ()                       { ftime_ms = *(std::min_element(ft_ms.begin(), ft_ms.end())); }
    // fix the event time
    void         SetEventGivenTime ( double time_ms)        { ftime_ms = time_ms; }
    // define the event time as the average time of all the SiPM hits
//    void            SetEventCharge ()                       { fQtot = std::accumulate(fQ.begin(), fQ.end(), 0.0);  }
////    // fix the event charge
//    void       SetEventGivenCharge ( double Qtot)        { fQtot = Qtot; }
    
    // set prompt and total charge
    void                   SetSiPM ();
    void                SetQprompt ();
    void                 SetQtotal ();
    void            SetTailToTotal ();

    // add or remove hit
    void                    AddHit ( double hit_time_ms, double hit_charge, int hit_channel );
    void                 RemoveHit ( int hitIdx );
    
    // auxiliary
    //    int     DetermineWhichDetector (std::vector<int> channels, std::string setup="unknown");
    //    int    ChannelNumberToDetector ( int channel, std::string setup="unknown") ;
    int         DetermineWhichSiPM ( std::vector<int> channels, std::string setup="unknown" );
    int        ChannelNumberToSiPM ( int channel, std::string setup="unknown" );
    
    
    // Time window for collection of a SiPM channel-pair
    // According to Luis, there should be 2-3 ns time variation between the channels
    // so we open a time window of +/- 5 ns ( = 5.e-6 ms )
    double       GetPairTimeWindow ()                        { return 5.e-6; } ;

    // check if event passes minimal requirements
    bool        CheckIfPairIsGood ( int NSiPMs = 2 );
    int         verbose;

    
    
    // SiPM signal is split between two channels in the ASIC.
    //
    // Using table-10 of the FEB/s 2ro datasheet
    // we convert the channel number to SiPM
    //
    std::vector<std::vector<int>> SiPM_ch1_ch2_map = {
        {1, 16, 18},
        {9, 17, 19},
        {2, 21, 20},
        {10,23, 22},
        {3, 25, 24},
        {11,26, 28},
        {4 ,29, 30},
        {12,31, 32},
        {5 ,34, 33},
        {13,36, 35},
        {6 ,37, 39},
        {14,41, 40},
        {7 ,43, 42},
        {15,44, 46},
        {8 ,45, 47},
        {16,48, 49},
        {17, 1,  0},
        {25,10,  3},
        {18,12,  5},
        {26,14,  7},
        {19, 8, 15},
        {27, 6,  4},
        {20,13,  2},
        {28,11, 27},
        {21,38,  9},
        {29,62, 58},
        {22,61, 59},
        {30,50, 57},
        {23,52, 51},
        {31,54, 53},
        {24,56, 55},
        {32,60, 63}
    };

    
    
private:
    
    bool                fIsGoodPair;
    int                 feventID;
    int                 fSiPM;
    std::vector<int>    fchannels;  // channel numbers
    std::vector<double> ft_ms;      // time of each channel fire in ms
    std::vector<double> fQ;         // charge of each channel fire in ADC units
    double              ftime_ms;   // time of event in ms (convert here to ms since ns gives too-large numbers)
    double              fQprompt;
    double              fQtail;
    double              fQtotal;
    double              fTailToTotal;
    int                 fN;
    int                 fdetector;
    std::string         fsetup;
    
    detectorEvent         detEv;
    std::vector<SiPMHit>  fHits;    // charge of each channel fire in ADC units
    
};

#endif

