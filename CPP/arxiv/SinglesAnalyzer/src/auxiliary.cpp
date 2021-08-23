#include <iostream>
#include <auxiliary.hpp>

// csv files
// read group csv-file
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector< std::vector<double> > auxiliary::read_csv(std::string filename){
    // Reads a CSV file into a vector of <string, vector<int>> pairs where
    // each pair represents <column name, column values>
    
    // Create a vector of <string, int vector> pairs to store the result
    std::vector< std::vector<double> > result;
    std::vector<double> line_vals;
    
    
    // Create an input filestream
    std::ifstream myFile(filename);
    
    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");
    
    // Helper vars
    std::string line, colname;
    double val;
    // Read data, line by line
    while(std::getline(myFile, line))
    {
        // std::cout << "line:" << line << std::endl;
        // Create a stringstream of the current line
        std::stringstream ss(line);
        line_vals.clear();
        // Extract each integer
        while(ss >> val){
            
            line_vals.push_back(val);
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
        }
        result.push_back(line_vals);
        if ( verbose > 2 ){
            for (size_t i=0; i<line_vals.size(); i++){
                std::cout << line_vals.at(i) << ",";
            }
            std::cout << std::endl;
        }
    }
    
    // Close file
    myFile.close();
    
    return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// open and write event csv-file
void auxiliary::open_events_csv ( std::string filename,std::string header ){
    Debug(0,"auxiliary::open_events_csv()");
    
    csvfilename = csvpath + "/" + filename + ".csv";
    csvfile.open( csvfilename );
    csvfile << header << std::endl;
    
    Debug(0,"opened output csv: \n" + csvfilename);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::write_events_csv ( std::vector<double> values ){
    for (size_t i=0; i<values.size()-1; i++){
        csvfile << values.at(i) << ",";
    }
    csvfile << values.at(values.size()-1);
    csvfile << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::close_events_csv (){
    csvfile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::PrintSummary (){
    std::cout
    << std::setprecision(1) << std::fixed
    << std::endl << "--------------------"
    << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<detectorEvent> auxiliary::CollectEvents(std::vector< std::vector<double> > SiPMgroups){
    std::vector<detectorEvent> events;
    
    //  csv columns: 'N(SiPMs)','n(SiPM)','time','charge','channel'
    std::vector<double> t_ms;
    std::vector<int>    channels;
    std::vector<float>  Q;
    
    int eventNumber = 0;
    for (size_t rowIdx=0; rowIdx<SiPMgroups.size(); rowIdx++){
        
        int     N       = int(SiPMgroups.at(rowIdx).at(0));
        int     n       = int(SiPMgroups.at(rowIdx).at(1));
        // convert to time in ms here since the long numbers (PETsys group data time are given in ps) are too long for the computer to digest
        double  time_ms = double(SiPMgroups.at(rowIdx).at(2))/1.e9;
        float   charge  = SiPMgroups.at(rowIdx).at(3); // ToDo: add lineariseChargeDeposited
        int     ch      = SiPMgroups.at(rowIdx).at(4);
        
        
        if (verbose>2){
            std::cout<< "N: " << N << ", n: " << n << ", ch: " << ch  << std::setprecision(14) << ", time: " << time_ms << " ms" << std::endl ;
        }
        
        t_ms.push_back(time_ms);
        channels.push_back(ch);
        Q.push_back(charge);
        
        if (n==N-1) {
            // create detector event
            detectorEvent event( eventNumber, channels, t_ms , Q, N, fsetup) ;
            if (verbose>2) event.Print();
            events.push_back(event);
            eventNumber++;
            Q.clear();
            t_ms.clear();
            channels.clear();
        }
    }
    return events;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<detectorEvent> auxiliary::SeparateAndFilterEvents(std::vector<detectorEvent> collected_events){
    
    std::vector<detectorEvent> separated_events;
    
    for (auto event: collected_events){
        int eventID = event.GetEventID();
        
        if (event.GetDetector()==0){ // separate the event
            std::vector<detectorEvent> splitted_events;
            if (verbose>1){
                std::cout << "separating event " << event.GetEventID() << " of detector " << event.GetDetector() << std::endl;
                event.Print();
            }
            
            detectorEvent original_event = detectorEvent(event.GetEventID(),
                                                         event.GetChannels(),
                                                         event.GetHitsTime(),
                                                         event.GetHitsCharge() ,
                                                         event.GetNhits(),
                                                         fsetup );
            
            std::vector<int> event_channels = event.GetChannels();
            std::vector<int> HitsRetained;
            
            // define the correct event detector as the one for the first channel
            int event_detector = event.ChannelNumberToDetector(event_channels.at(0));
            int separated_event_id = 0;
            if (verbose>1){
                std::cout << "event detector is " << event_detector << std::endl;
            }
            
            for (auto hitIdx=0; hitIdx<original_event.GetNhits(); hitIdx++ ) {
                int hit_channel     = event_channels.at(hitIdx);
                int hit_detector    = original_event.ChannelNumberToDetector( hit_channel, fsetup );
                double hit_time_ms  = original_event.GetHitTime(hitIdx); // time of event in ms
                float hit_charge    = original_event.GetHitCharge(hitIdx);

                // if this hit does not belong to the same detector, separate it
                if ( hit_detector != event_detector){
                    if (verbose>1){
                        std::cout << "separating hit " << hitIdx << ", channel " << hit_channel << " detector " << hit_detector << std::endl;
                    }
                
                    // now add this hit into another event, either a new or an existing one
                    Debug(1,"check if hits in this detector was already splitted into another event");
                    // check if this detector was already separated into another event
                    // and if not, create a new separated event
                    bool added_hit_to_another_detector = false;
                    for (auto &splitted_event:splitted_events){
                        if ( hit_detector == splitted_event.GetDetector() ){
                            splitted_event.AddHit( hit_time_ms, hit_charge, hit_channel );
                            added_hit_to_another_detector = true;
                            
                            if(verbose>1) std::cout << "added hit to a separated event " << splitted_event.GetEventID() << std::endl;
                        }
                    }
                    if (added_hit_to_another_detector == false){
                        separated_event_id ++;
                        detectorEvent separated_event = detectorEvent();
                        separated_event.SetSetup( fsetup );
                        separated_event.SetDetector( hit_detector );
                        separated_event.SetEventID( event.GetEventID()*100 + separated_event_id );
                        separated_event.AddHit( hit_time_ms, hit_charge, hit_channel );
                        splitted_events.push_back( separated_event );
                        if(verbose>1) std::cout << "created a new separated event " << separated_event.GetEventID() << std::endl;
                    }
                    
                }
                else{
                    // finally, remove this hit from the hit collection in this event
                    Debug(1,"retain this hit from the original hit collection of this event");
                    HitsRetained.push_back(hitIdx);

                    if (verbose>1) {
                        std::cout
                        << "channel " << event_channels.at(hitIdx)
                        << " is from detector " << event.ChannelNumberToDetector( event_channels.at(hitIdx) )
                        << std::endl;
                    }
                }
            }
            
            detectorEvent reducedEvent = detectorEvent();
            reducedEvent.SetEventID( original_event.GetEventID()*100 );
            reducedEvent.SetDetector( event_detector );
            for (int hitIdx:HitsRetained) {
                reducedEvent.AddHit( original_event.GetHitTime(hitIdx), original_event.GetHitCharge(hitIdx), original_event.GetHitChannel(hitIdx));
            }
            splitted_events.push_back( reducedEvent );
            
            if (verbose>1){
                std::cout << "splitted" << std::endl;
                original_event.Print();
                std::cout << "into" << std::endl;
                for (auto splitted_event:splitted_events) splitted_event.Print();
            }
            
            // after separating all the events from the one which detector was "unknown",
            // plug the separated evennts into the results
            for (auto splitted_event:splitted_events){
                if ( splitted_event.CheckIfEventIsGood() ) {
                    separated_events.push_back(splitted_event);
                }
            }
            splitted_events.clear();
        } // end if (event.GetDetector()==0)
        // in case there was no need to split this event into multiple ones,
        // plug the separated evennts into the results
        else {
            separated_events.push_back(event);
        }
        if (verbose>2) {
            std::cout << "done stepping through event " << eventID << std::endl;
            std::cout << " ------------------------------------------------------- " << std::endl;
        }
        
    }
    return separated_events;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::StreamEventsToCSV (std::vector<detectorEvent> events, std::string filename, std::string header){

    open_events_csv(filename, header);
    for (auto event: events){
        write_event_to_csv( event );
    }
    close_events_csv();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::write_event_to_csv ( detectorEvent event ){
    csvfile
    << event.GetEventID()      << ","
    << event.GetNhits()        << ","
    << std::setprecision(12) << event.GetEventTime()    << ","
    << event.GetEventCharge()  << ","
    << event.GetDetector()     << ",";
    
    csvfile << "[";
    for (auto ch:event.GetChannels()){
        csvfile << ch << ";";
    }
    csvfile << "]";
    csvfile << std::endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<double> auxiliary::CollectDetectionTimeDifferencesArray(std::vector<detectorEvent> events, double dt_max_ms){
    // we convert to ms and work in ms here since the numbers in ns are too large for the machine to digest
    std::vector<double> time_differences_ms;
    std::vector<double> event_times_ms;
    
    for (auto event: events) {
        double t_ms = double(event.GetEventTime());
        event_times_ms.push_back(t_ms);
    }
    std::sort( event_times_ms.begin() , event_times_ms.end() );
        
    if (verbose>2){
        std::cout << "event times: [";
        for (auto t: event_times_ms) std::cout << t << ",";
        std::cout << "]" << std::endl;
    }
    
    for (size_t evtIdx=0; evtIdx<events.size(); evtIdx++) {
        double t = event_times_ms.at(evtIdx);
        // calculate time difference to next events within a few ms away, or 400 events
        // the reason for the limit on 400 events is that a typical fission source we are using emits < 37k fissions per sec,
        // which means a fission every 27 us
        // so even if we detect 10 neutrons and gammas per fission (practically impossible due to our inefficeincy) - we read
        // on average an event every 2.7 us.
        // 100 events are thus 270 micro-seconds, and 400 events is around 1 ms
        for (size_t next_evtIdx = evtIdx+1; next_evtIdx < std::min(evtIdx+400, events.size()); next_evtIdx++){
            double dt = event_times_ms.at(next_evtIdx) - t;
            if ((dt==0) && (verbose>0)){
                std::cout
                << "CollectDetectionTimeDifferencesArray() dt=0 between events "
                << evtIdx << " and " << next_evtIdx << endl;
            }
            if (dt < dt_max_ms){
                time_differences_ms.push_back(dt);
            } else { // if we got farther than a few ms away (dt_max_ms) no need to continue the loop
                break;
            }
        }
    }
    if (verbose>2){
        std::cout << "event time differences: [";
        for (auto dt: time_differences_ms) std::cout << dt << ",";
        std::cout << "]" << std::endl;
    }

    if (verbose>0){
        std::cout << "collected time difference vector of size " << time_differences_ms.size() << " from " << events.size() << " events ";
    }
    return time_differences_ms;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::StreamTimeDifferencesToCSV (std::vector<double> time_differences_ms, std::string filename){

    csvfilename = csvpath + "/" + filename + ".csv";
    csvfile.open( csvfilename );
    csvfile << "dt[ms]" << std::endl;
    Debug(0,"opened output csv: \n" + csvfilename);
    for (auto dt: time_differences_ms){
        csvfile << dt << std::endl;
    }
    // csvfile << std::endl;
    csvfile.close();
}
