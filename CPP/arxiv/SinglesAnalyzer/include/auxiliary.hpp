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
#include </Users/erezcohen/Desktop/PETsys/Software/PETsysAnalysis/CPP/GroupAnalyzer/include/detectorEvent.hpp> // std::runtime_error

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
    std::string fsetup;
    std::string csvpath;
    std::string csvfilename;
    std::ofstream csvfile;

    // contructor
    auxiliary(){
        SetSetup("unknown");
        SetCSVPath();
    };
    
    // prints
    void Debug (int verobosity_level, std::string text) { if ( verbose > verobosity_level ) std::cout << text << std::endl; }
    
    // csv files
    // read group csv-file
    std::vector< std::vector<double> > read_csv (std::string filename);

    // open and write event csv-file
    void   StreamEventsToCSV ( std::vector<detectorEvent> events,
                              std::string filename = "events",
                              std::string header = "eventID,N(SiPMs),time[ms],Qtot[a.u.],detector,channels");
    void     open_events_csv ( std::string filename, std::string header);
    
    void  write_event_to_csv ( detectorEvent event );
    void    write_events_csv ( std::vector<double> values );
    void    close_events_csv ();
    
    void        PrintSummary ();
    
    
    // Setters
    void          SetVerbose (int v)            { verbose = v;};
    void            SetSetup (std::string setup){ fsetup = setup; SetCSVPath();};
    void          SetCSVPath ()                 { csvpath = "/Users/erezcohen/Desktop/data/PETsys/"+fsetup+"/";};
    
    
    
    // trigonometry
    double           rad2deg (double angle_rad) {return angle_rad*180./3.1415;};
    double           deg2rad (double angle_deg) {return angle_deg*3.1415/180.;};

    
    // collect events from SiPM group data
    std::vector<detectorEvent>                CollectEvents (std::vector< std::vector<double> > SiPMgroups);
    
    // separate events that include SiPMs from different detectors and filter "good" events
    std::vector<detectorEvent>      SeparateAndFilterEvents (std::vector<detectorEvent> collected_events);
    
    
    std::vector<double> CollectDetectionTimeDifferencesArray (std::vector<detectorEvent> events,
                                                                 double dt_max_ms = 10 ); // maximal time we consider for dt is 1 ms
    void                           StreamTimeDifferencesToCSV (std::vector<double> time_differences, std::string filename);
};




#endif

