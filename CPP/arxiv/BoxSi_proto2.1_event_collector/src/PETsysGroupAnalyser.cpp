/**

 BoxSi proto2.1 group analyser / event collector
 
 last edit Feb-18, 2021 (EOC)

 operate:
 -------------
 
 ./build/apps/PETsysGroupAnalyser <SourceStr="Cf252"> <dT="10sec"> <verbose=1>
 ./build/apps/PETsysGroupAnalyser Cf252 10sec 3
 
 */

#include <data_module/detectorEvent.hpp>
#include <auxiliary.hpp>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int main(int argc, char **argv){

    
    std::string SourceStr = "Cf252", dTStr = "10sec";
    int verbose = 1;
    if (argc>1) {
        SourceStr = argv[1];
        if (argc>2) {
            dTStr = argv[2];
            if (argc>2) {
                verbose = std::stoi(argv[3]);
            }
        }
    }
    std::string filelabel = SourceStr+"_"+dTStr;
    auxiliary * aux = new auxiliary();
    aux -> SetVerbose(verbose);
    
    // (1) read SiPM group data
    std::cout << std::endl << "(1) read SiPM group data from file" << std::endl << aux->csvpath + "/" + filelabel + "_group.dat" << std::endl;
    std::vector< std::vector<double> > SiPMgroups = aux->read_csv( aux->csvpath+"/" + filelabel + "_group.dat" );
    
    // (2) collect detector events
    std::cout << std::endl << "(2) collect detector events" << std::endl << std::endl;
    std::vector<detectorEvent> collected_events = aux->CollectEvents(SiPMgroups);
    
    
    // (3) loop over the events
    //      (A) if they include channels from multiple detectors, separate them
    //      (B) filter only "good" events
    std::cout << std::endl << "(3) loop over the events, if they include channels from multiple detectors, separate them" << std::endl<< std::endl;
    std::vector<detectorEvent> events = aux->SeparateAndFilterEvents(collected_events);
    
    // (4) stream separated events into output csv file
    std::cout << std::endl << "(4) stream separated events into output csv file" << std::endl<< std::endl;
    aux -> StreamEventsToCSV( events , filelabel + "_events" );
    
    // (5) compile array of time difference from each detection to form Rossi-alpha distribution
    std::cout << std::endl << "(5) compile array of time difference from each detection to form Rossi-alpha distribution" << std::endl<< std::endl;
    std::vector<double> time_differences = aux -> CollectDetectionTimeDifferencesArray( events );
    aux -> StreamTimeDifferencesToCSV( time_differences, filelabel + "_time_differences" );
        
    std::cout << std::endl << "Done." << std::endl<< std::endl;
    
    return 0;
}



