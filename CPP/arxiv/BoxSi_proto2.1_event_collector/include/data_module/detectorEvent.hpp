#ifndef DETECTOREVENT_HPP
#define DETECTOREVENT_HPP

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

class detectorEvent
{
public:
    
    // time in ms in this class
    detectorEvent(){fN = 0;};
    detectorEvent(int eventID, std::vector<int> channels, std::vector<double> t_ms, std::vector<float> Q, int N);
    
    
    int     DetermineWhichDetector (std::vector<int> channels);
    void                     Print ();
    
    // Getters
    int                      GetNhits ()                       {return fN;};
    int                    GetEventID ()                       {return feventID;};
    int                   GetDetector ()                       {return fdetector;};
    int                GetEventCharge ()                       {return fQtot;};
    bool                  IsGoodEvent ()                       {return fIsGoodEvent;};
    double               GetEventTime ()                       {return ftime_ms;}; // return time in ms
    
    int                 GetHitChannel (int hitIdx)             {return fchannels.at(hitIdx);};
    std::vector<int>      GetChannels ()                       {return fchannels;};
    
    double                 GetHitTime (int hitIdx)             {return ft_ms.at(hitIdx);};
    std::vector<double>   GetHitsTime ()                       {return ft_ms;};
    
    float                GetHitCharge (int hitIdx)             {return fQ.at(hitIdx);};
    std::vector<float>  GetHitsCharge ()                       {return fQ;};
    

    
    

    // Setters
    void                SetEventID (int eventID)            {feventID = eventID;};
    void               SetDetector (int detector)           {fdetector = detector;};
    // define the event time as the time of the first SiPM hit in the event
    void              SetEventTime ()                       { ftime_ms = *(std::min_element(ft_ms.begin(), ft_ms.end())); }
    // define the event time as the average time of all the SiPM hits
//    void              SetEventTime ()                       { ftime_ms = std::accumulate(ft_ms.begin(), ft_ms.end(), 0.0) / ft_ms.size(); }
    void            SetEventCharge ()                       { fQtot = std::accumulate(fQ.begin(), fQ.end(), 0);  }

    // add or remove hit
    void                    AddHit ( double hit_time_ms, float hit_charge, int hit_channel );
    void                 RemoveHit ( int hitIdx );
    // auxiliary
    int    ChannelNumberToDetector (int channel);
    // check if event passes minimal requirements
    bool        CheckIfEventIsGood (int NSiPMs_min=2);

    
private:
    
    bool            fIsGoodEvent;
    int                 feventID;
    std::vector<int>    fchannels;  // channel numbers
    std::vector<double> ft_ms;       // time of each channel fire in ms
    std::vector<float>  fQ;         // charge of each channel fire in ADC units
    double  ftime_ms;               // time of event in ms (convert here to ms since ns gives too-large numbers)
    float   fQtot;
    int     fN;
    int     fdetector;
    
    
};

#endif
