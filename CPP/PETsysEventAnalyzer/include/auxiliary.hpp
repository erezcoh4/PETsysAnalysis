#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

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
#include <algorithm>
#include <math.h>       /* fabs */

#include <detectorEvent.hpp> // std::runtime_error
#include <SiPMChannelPair.hpp> // std::runtime_error


#define ps 1.
#define ns 1000.
#define us 1000000.
#define ms 1000000000.

#define PrintEmptyLine(){ std::cout << std::endl;}
#define PrintLine(){ std::cout << "------------------------------------------------" << std::endl;}
#define PrintXLine(){ std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;}
#define PrintTVector3(v){ std::cout <<std::setprecision(1)<<std::fixed << #v<< ": ("<<v.X()<<","<<v.Y()<<","<<v.Z()<<")"<< std::endl;}
#define Print4TVector3(v){ std::cout <<std::setprecision(4)<<std::fixed << #v<< ": ("<<v.X()<<","<<v.Y()<<","<<v.Z()<<")"<< std::endl;}

using namespace std;

class auxiliary
{
private:
    
    
public:
    
    int verbose;
    std::string     fsetup, fsubdir;
    std::string     csvpath;
    std::string     csvfilename;
    std::ofstream   csvfile;

    // contructor
    auxiliary(){
        SetSetup("unknown");
        SetCSVPath();
    };
    
    // prints
    void Debug (int verobosity_level, std::string text) { if ( verbose > verobosity_level ) std::cout << text << std::endl; }
    
    
    
    // csv files
    // read group csv-file and return a vector of vector of doubles
    std::vector< std::vector<double> >  ReadCSV ( std::string filename, int Nrows=-1 );
    // results are vector of "rows", each includes the following:
    // if DataType == single,   PETsysData vector = {'time [ps]','charge','channel'}
    // if DataType == group,    PETsysData vector = {'N(SiPMs)','n(SiPM)','time [ps]','charge','channel'}
    // if DataType == group/PSD,PETsysData vector = {'N(SiPMs)','n(SiPM)','time [ps]','charge','channel'}
    
    std::vector< detectorEvent >   ReadEventsCSV ( std::string filename, int Nrows=-1 );
    // if DataType == events,   Events vector     = {'eventID','N(SiPMs)','time[ms]','Qtot[a.u.]','detector','channels'}

    
    
    
    // open and write event csv-file
    void   StreamEventsToCSV ( std::vector<detectorEvent> events,
                              std::string filename = "events",
                              std::string header = "eventID,N(SiPMs),time[ms],Qtot[a.u.],detector,channels");
    void   StreamEventsToCSV ( std::vector<SiPMChannelPair> eventsWithPSD,
                              std::string filename = "eventsWithPSD",
                              std::string header = "eventID,SiPM,time[ms],Q(prompt)[a.u.],Q(total)[a.u.],tail/total,detector,channels");
    void       OpenEventsCSV ( std::string filename, std::string header);
    
    void     WriteEventToCSV ( detectorEvent    event );
    void     WriteEventToCSV ( SiPMChannelPair  event );
    void    write_events_csv ( std::vector<double> values );
    void    close_events_csv ();
    
    void        PrintSummary (std::vector<detectorEvent>);
    
    
    
    
    
    // Setters
    void          SetVerbose (int v)                { verbose = v;};
    void            SetSetup (std::string setup)    { fsetup = setup;   SetCSVPath();};
    void           SetSubdir (std::string subdir)   { fsubdir = subdir; SetCSVPath();};
    void          SetCSVPath ()                     { csvpath = "/Users/erezcohen/Desktop/data/PETsys/"+fsetup+"/"+fsubdir+"/";};
    
    
    // collect events from SiPM group data
    std::vector<detectorEvent>                    CollectEvents (std::vector< std::vector<double> > PETsysData,
                                                                 std::string DataType="single");
    std::vector<detectorEvent>         CollectEventsFromSingles (std::vector< std::vector<double> > PETsysData);
    std::vector<detectorEvent>          CollectEventsFromGroups (std::vector< std::vector<double> > PETsysData);
    
    std::vector<SiPMChannelPair>                 CollectHitsPSD (std::vector< std::vector<double> > PETsysData,
                                                                 std::string DataType="single_PSD");
    std::vector<SiPMChannelPair>  CollectHitPairsPSDFromSingles (std::vector< std::vector<double> > PETsysData);
    std::vector<SiPMChannelPair>   CollectHitPairsPSDFromGroups (std::vector< std::vector<double> > PETsysData);

    // separate events that include SiPMs from different detectors and filter "good" events
    std::vector<detectorEvent>          SeparateAndFilterEvents (std::vector<detectorEvent> collected_events);
    std::vector<SiPMChannelPair>        SeparateAndFilterEvents (std::vector<SiPMChannelPair> collected_events);

    
    std::vector<double>    CollectDetectionTimeDifferencesArray (std::vector<detectorEvent> events,
                                                                 double dt_max_ms = 10);             // maximal time we consider for dt is 10 ms
    
    
    void                           ExtractSinglesDoublesTriples (std::vector<detectorEvent> events,
                                                                 double R_gate = 100.e-6,           // gate for doubles and triples - 100 ns
                                                                 double RA_gate_delay = 100.e-3);   // gate delay for accidentals - 100 us
    double                                         LinearizeQDC ( double Q=0 ) ;
    void                             StreamTimeDifferencesToCSV (std::vector<double> time_differences, std::string filename);
    
    
    
    
    
    
    // trigonometry
    double           rad2deg (double angle_rad) {return angle_rad*180./3.1415;};
    double           deg2rad (double angle_deg) {return angle_deg*3.1415/180.;};
};




#endif

