#ifndef SIPMHIT_HPP
#define SIPMHIT_HPP

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

class SiPMHit
{
public:
    
    // time in ms in this class
    SiPMHit ()                                      {fhitID = -1; fchannel = -1; ft_ms = -1.; fQ = -1.;};
    SiPMHit (int hitID, int channel,
             double t_ms, double Q);
    
    
    void                     Print ();
    
    // Getters
    int                   GetHitID ()               {return fhitID;};
    int                 GetChannel ()               {return fchannel;};
    double                 GetTime ()               {return ft_ms;};
    double               GetCharge ()               {return fQ;};
    

    
    // Setters
    void                  SetHitID (int id)         {fhitID=id;};
    void                SetChannel (int ch)         {fchannel=ch;};
    void                   SetTime (double t)       {ft_ms=t;};
    void                 SetCharge (double Q)       {fQ=Q;};
    
    
private:
    
    int     fhitID;
    int     fchannel;   // channel numbers
    double  ft_ms;      // time of in ms
    double  fQ;         // charge in ADC units
    
    
};

#endif

