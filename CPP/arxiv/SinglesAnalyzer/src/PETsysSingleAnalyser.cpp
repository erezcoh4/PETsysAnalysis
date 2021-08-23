/**

 general PETsys group analyser / event collector
 
 last edit Mar-8, 2021 (EOC)

 operate:
 -------------
 
 ./build/apps/PETsysGroupAnalyser <setup> <DataLabel> <dT> <verbose>
 
 examples:
 
 ./build/apps/PETsysGroupAnalyser BoxSi_proto2.1 Cf252 10sec 3
 ./build/apps/PETsysGroupAnalyser PlasticScintillator_SensL8x8 PlasticScintillator_Cf252_1mm_vth1_20_vth2_20 30_sec 3
 ./build/apps/PETsysGroupAnalyser PlasticScintillator_SensL8x8 PlasticScintillator_cosmic_vth1_20_vth2_20 30_sec 3

 */

#include </Users/erezcohen/Desktop/PETsys/Software/PETsysAnalysis/CPP/GroupAnalyzer/include/detectorEvent.hpp> // std::runtime_error
#include <auxiliary.hpp>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int main(int argc, char **argv){

    
    std::string     setup = "BoxSi_proto2.1";
    std::string DataLabel = "Na22_10uC";
    std::string     dTStr = "10_sec";
    int verbose = 1;
    if (argc>1) {
        setup = argv[1];
        if (argc>2) {
            DataLabel = argv[2];
            if (argc>3) {
                dTStr = argv[3];
                if (argc>4) {
                    verbose = std::stoi(argv[4]);
                }
            }
        }
    }
    std::string filelabel = DataLabel+"_"+dTStr;
    auxiliary * aux = new auxiliary();
    aux -> SetVerbose( verbose );
    aux -> SetSetup( setup );
    
    // (1) read SiPM group data
    std::cout   << std::endl << "(1) read SiPM singles data from file" << std::endl
                << aux->csvpath + "/" + filelabel + "_single.dat" << std::endl;
    std::vector< std::vector<double> > SiPMsingles = aux->read_csv( aux->csvpath+"/" + filelabel + "_single.dat" );
    
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
    // this takes some time and is not required usually
    bool doTimeDifferences = false;
    if (doTimeDifferences){
        std::cout << std::endl << "(5) compile array of time difference from each detection to form Rossi-alpha distribution" << std::endl<< std::endl;
        std::vector<double> time_differences = aux -> CollectDetectionTimeDifferencesArray( events );
        aux -> StreamTimeDifferencesToCSV( time_differences, filelabel + "_time_differences" );
    }
    
    std::cout << std::endl << "Done." << std::endl<< std::endl;
    
    return 0;
}



