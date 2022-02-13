

#include <detectorEvent.hpp>
#include <auxiliary.hpp>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int main(int argc, char **argv){

    bool doTimeDifferences  = false;
    bool doMonitorEvents    = false;
    bool doExtractSDT       = false; // extract singles, doubles, triples
    
    std::string     setup = "BoxSi_proto2.2";
    std::string    subdir = "vth_1PE";
    std::string DataLabel = "Cf252";
    std::string     dTStr = "10_sec";
    std::string  DataType = "single";
    int           verbose = 1;
    int             Nrows = -1; // number of rows to process
    
    if (argc>1) {
        setup = argv[1];
        if (argc>2) {
            subdir = argv[2];
            if (argc>3) {
                DataLabel = argv[3];
                if (argc>4) {
                    dTStr = argv[4];
                    if (argc>5) {
                        DataType = argv[5];
                        if (argc>6) {
                            verbose = std::stoi(argv[6]);
                            if (argc>7) {
                                Nrows = std::stoi(argv[7]);
                            }
                        }
                    }
                }
            }
        }
    }
    auxiliary * aux = new auxiliary();
    aux -> SetVerbose( verbose );
    aux ->  SetSetup ( setup );
    aux -> SetSubdir ( subdir );
    
    
    std::string filelabel = DataLabel + "_" + dTStr;
    std::string filename  = aux->csvpath + "/" + filelabel + "_" + DataType;
    std::vector<detectorEvent> events;
    
    if (verbose>3){
        std::cout
        << "filelabel: "    << filelabel
        << std::endl
        << "filename: "     << filename
        << std::endl
        << "DataType: "     << DataType
        << std::endl
        << "Nrows: "        << Nrows
        << std::endl;
    }
    
    if (DataType.compare("single")==0 || DataType.compare("group")==0){
        // if its a PETsys file (single / group), process the SiPM hits and generate a CSV of "events"
        filename  = aux->csvpath + "/" + filelabel + "_" + DataType + ".dat";
        
        // (1) read SiPM (PETsys) data
        std::cout   << std::endl << "(1) read SiPM " << DataType << " data from file" << std::endl << filename << std::endl;
        std::vector< std::vector<double> > PETsysData = aux->ReadCSV( filename, Nrows );
        
        // (2) collect detector events
        std::cout << std::endl << "(2) collect detector events" << std::endl << std::endl;
        std::vector<detectorEvent> collected_events = aux->CollectEvents( PETsysData, DataType );
            
        // (3) filter good and events separate events from multiple detectors
        std::cout << std::endl << "(3) filter good and events separate events from multiple detectors" << std::endl<< std::endl;
        events = aux->SeparateAndFilterEvents(collected_events);
        
        // (4) stream separated events into output csv file
        std::cout << std::endl << "(4) stream separated events into output csv file" << std::endl<< std::endl;
        aux -> StreamEventsToCSV( events , filelabel + "_events" );

    } else if (DataType.compare("events")==0){
        // if the input file is already sorted to "events"
        // we can skip the process of forming events,
        // and directly go to the next stages
        filename  = aux->csvpath + "/" + filelabel + "_" + DataType + ".csv";
        events = aux->ReadEventsCSV( filename, Nrows );
    }
    
    
    
    
    // (5) compile array of time difference from each detection to form Rossi-alpha distribution
    // this takes some time and is not required usually
    if (doTimeDifferences){
        std::cout
        << std::endl
        << "(5) compile array of time difference from each detection to form Rossi-alpha distribution"
        << std::endl
        << std::endl;
        std::vector<double> time_differences = aux -> CollectDetectionTimeDifferencesArray( events );
        aux -> StreamTimeDifferencesToCSV( time_differences, filelabel + "_time_differences" );
    }
    
    // (6) extract s,d,t
    if (doExtractSDT){
        aux -> ExtractSinglesDoublesTriples(events);
    }
    
    // (7) print interesting event characteristics for monitoring
    if (doMonitorEvents){
        aux->PrintSummary( events );
    }
    
    std::cout << std::endl << "Done." << std::endl<< std::endl;
    
    return 0;
}



