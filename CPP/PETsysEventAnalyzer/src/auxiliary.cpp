#include <iostream>
#include <auxiliary.hpp>
#define ms_to_us 1e3
#define ms_to_ns 1e6

// csv files
// read PetSys (single/group) csv-file
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector< std::vector<double> > auxiliary::ReadCSV( std::string filename, int Nrows ){
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
    int lineNumber=0;
    // Read data, line by line
    while(std::getline(myFile, line)) {
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
        lineNumber++;
        // stop after N lines
        if (Nrows>0 && lineNumber>Nrows) {
            break;
        }
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


// csv files
// read Events-csv file that was either created by this app or using a MC
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector< detectorEvent > auxiliary::ReadEventsCSV( std::string filename, int Nrows ){
    // read a CSV file of
    //    {'eventID','N(SiPMs)','time[ms]','Qtot[a.u.]','detector','channels'}
    std::vector< detectorEvent > events;
    detectorEvent event;
    
    std::vector<double> line_vals;
    
    // Create an input filestream and make sure the file is open
    std::ifstream myFile(filename);
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");
    
    // Helper vars
    std::string line, colname;
    double val;
    int lineNumber=0;
    // Read first line - header
    std::getline(myFile, line);
    
    // Read data, line by line
    while(std::getline(myFile, line)) {
        std::stringstream ss(line);
        if ( verbose > 5 ) std::cout << line << std::endl;
        event = detectorEvent();
        
        // read line values
        line_vals.clear();
        while(ss >> val){
            line_vals.push_back(val);
            if(ss.peek() == ',') ss.ignore();
        }
        // plug into "event"
        event.SetEventID            ( line_vals.at(0) );
        event.SetNhits              ( line_vals.at(1) );
        event.SetEventGivenTime     ( line_vals.at(2) );
        event.SetEventGivenCharge   ( line_vals.at(3) );
        event.SetDetector           ( line_vals.at(4) );
        // and we do not set the individual SiPM hits at this stage,
        // as this implementation is seemingly a little challanging and unnecessary...
        
        
        events.push_back( event );
        lineNumber++;
        // stop after N lines
        if ( verbose > 2 )                  event.Print();
        if (Nrows>0 && lineNumber>(Nrows-1) )    break;
    }
    // Close file
    myFile.close();
    return events;
    
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// open and write event csv-file
void auxiliary::OpenEventsCSV ( std::string filename,std::string header ){
    Debug(0,"auxiliary::OpenEventsCSV()" + csvfilename);
    
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
void auxiliary::PrintSummary(std::vector<detectorEvent> events){
    std::cout
    << std::setprecision(1) << std::fixed
    << std::endl << "--------------------"
    << std::endl;
    
    for (auto event:events) {
        bool didPrintSomething = false;
        // check if there are more than 64 hits in this event
        if (event.GetHits().size() > 64) {
            std::cout << event.GetHits().size() << " hits in event " << event.GetEventID() << std::endl;
            didPrintSomething = true;
        }

        // check if there is more than a single hit from each channel
        for (auto hit1:event.GetHits()) {
            for (auto hit2:event.GetHits()){
                if ((hit1.GetHitID() != hit2.GetHitID()) && (hit1.GetChannel() == hit2.GetChannel())) {
                    std::cout
                    << "in event "  << event.GetEventID()
                    << "hit "       << hit1.GetHitID() << " (time= " << hit1.GetTime() << " ms) "
                    << " and hit "  << hit2.GetHitID() << " (time= " << hit2.GetTime() << " ms) "
                    << " are from the same channel (" << hit1.GetChannel() << ")"
                    << std::endl;
                    didPrintSomething = true;
                }
            }
        }
        if (didPrintSomething) std::cout << "--------------------------------------------------------" << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<detectorEvent> auxiliary::CollectEvents(std::vector< std::vector<double> > PETsysData,
                                                    std::string DataType){
    std::vector<detectorEvent> events;
    if (DataType.compare("single")==0) {
        events = CollectEventsFromSingles( PETsysData );
    }
    else if (DataType.compare("group")==0) {
        events = CollectEventsFromGroups( PETsysData );
    }
    if (verbose>3) {
        std::cout << "done collecting events from singles:" << std::endl;
        for (auto event: events) event.Print();
        std::cout << std::endl;
    }
    return events;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<SiPMChannelPair> auxiliary::CollectHitsPSD(std::vector< std::vector<double> > PETsysData,
                                                    std::string DataType){
    std::vector<SiPMChannelPair> events;
    if (DataType.compare("single_PSD")==0) {
        events = CollectHitPairsPSDFromSingles( PETsysData );
    }
    else if (DataType.compare("group_PSD")==0) {
        events = CollectHitPairsPSDFromGroups( PETsysData );
    }
    if (verbose>3) {
        std::cout << "done collecting events from singles:" << std::endl;
        for (auto event: events) event.Print();
        std::cout << std::endl;
    }
    return events;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<detectorEvent> auxiliary::CollectEventsFromSingles(std::vector< std::vector<double> > PETsysData){
    std::vector<detectorEvent> events;
    // each row is a hit
    // for each hit, the columns are: 'time','charge','channel'
    std::vector<double> t_ms;
    std::vector<int>    channels;
    std::vector<double> Q;
    detectorEvent event;
    
    int eventID = 0;
    
    for (size_t rowIdx=0; rowIdx<PETsysData.size(); rowIdx++){
        // convert to time in ms here since the long numbers (PETsys group data time are given in ps) are too long for the computer to digest
        double  time_ms = double(PETsysData.at(rowIdx).at(0))/1.e9;
        double  charge  = PETsysData.at(rowIdx).at(1); // ToDo: add lineariseChargeDeposited
        int     channel = PETsysData.at(rowIdx).at(2);
        
        
        int detector = event.ChannelNumberToDetector( channel , fsetup );
        bool AddedHitToEvent = false;
        int Nevents = (int)events.size();
        
        // We first step over all Existing events in our list
        // and check if we can add this hit to one of them

        for (size_t evtIdx = std::max( 0, Nevents - 20 ); evtIdx < events.size(); evtIdx++) {
            detectorEvent & event = events.at(evtIdx);
            if (    (event.GetDetector() == detector)
                &&  (fabs(time_ms - event.GetEventTime()) < event.GetEventTimeWindow()) ) {
                event.AddHit( time_ms, charge, channel );
                AddedHitToEvent = true;
                break;
            }
        }
        
        // If there were no events to which this hit was added,
        // we create a new event

        if ( AddedHitToEvent == false ){
            // this means that we do not add hits to this event
            // but rather start a new event
            detectorEvent new_event;
            new_event.  SetEventID ( eventID );
            new_event.    SetSetup ( fsetup );
            new_event. SetDetector ( detector );
            new_event.      AddHit ( time_ms, charge, channel );
            events.push_back( new_event );
            eventID++;
        }
        

        if (verbose>2){
            std::cout
            << "channel: " << channel  << ", detector " << detector
            << std::setprecision(14) << ", time: " << time_ms << " ms"
            << ", charge: " << charge << " [a.u.]"<< std::endl ;
        }
    }

    return events;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<detectorEvent> auxiliary::CollectEventsFromGroups(std::vector< std::vector<double> > PETsysData){
    std::vector<detectorEvent> events;
    // each row is a hit
    // for each hit, the columns are: 'N(SiPMs)','n(SiPM)','time','charge','channel'
    std::vector<double> t_ms;
    std::vector<int>    channels;
    std::vector<double> Q;
    
    int eventNumber = 0;
    
    for (size_t rowIdx=0; rowIdx<PETsysData.size(); rowIdx++){
        
        int     N       = int(PETsysData.at(rowIdx).at(0));
        int     n       = int(PETsysData.at(rowIdx).at(1));
        // convert to time in ms here since the long numbers (PETsys group data time are given in ps) are too long for the computer to digest
        double  time_ms = double(PETsysData.at(rowIdx).at(2))/1.e9;
        double  charge  = PETsysData.at(rowIdx).at(3); // ToDo: add lineariseChargeDeposited
        int     ch      = PETsysData.at(rowIdx).at(4);
        
        
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
std::vector<SiPMChannelPair> auxiliary::CollectHitPairsPSDFromSingles(std::vector< std::vector<double> > PETsysData){
    //
    // collect pairs of channel hits from the same SiPM
    // according to table-10 of the FEB/s 2ro data-sheet,
    // and couple them together to get the prompt-charge and the total-charge
    // for the same SiPM hit
    // The pairing is done if the pair of channels saw a hit (crossed-threshold)
    // within 10 ns (1-3 ns variation between channels is expected due to wiring [Luis F., Feb-2022])
    //
    // For simplicity and to keep consistency with existing code above,
    // this routine is implemented in two steps:
    // (1) Collect "detector-events" like auxiliary::CollectEventsFromSingles()
    // (2) Split pairs of channels that read the same SiPM, into SiPMChannelPair objects
    //
    // return:
    // a vector of <SiPMChannelPair> objects to be written into a file
    //
    // Feb-16, 2022
    //
    
    std::vector<SiPMChannelPair>    hitpairs;
    // each row is a hit
    // for each hit, the columns are: 'time','charge','channel'
    std::vector<double> t_ms;
    std::vector<int>    channels;
    std::vector<double> Q;

    detectorEvent   event;
    SiPMChannelPair hitpair;
    
    int eventID = 0;
    for (size_t rowIdx=0; rowIdx<PETsysData.size(); rowIdx++){
        // convert to time in ms here since the long numbers (PETsys group data time are given in ps) are too long for the computer to digest
        double  time_ms = double(PETsysData.at(rowIdx).at(0))/1.e9;
        double  charge  = LinearizeQDC( PETsysData.at(rowIdx).at(1)) ;
        int     channel = PETsysData.at(rowIdx).at(2);
        int SiPM_number = hitpair.ChannelNumberToSiPM( channel );
        
        if (verbose>3){
            std::cout << "row " << rowIdx
            << ", time "        << std::setprecision(12) << time_ms << " ms "
            << ", channel "     << channel
            << ", SiPM_number " << SiPM_number
            << std::endl;
        }
        
        
        
        int detector = event.ChannelNumberToDetector( channel , fsetup );
        bool AddedHitToEvent = false;
        int Nhitpairs = (int)hitpairs.size();
        
        // We first step over all Existing events in our list
        // and check if we can add this hit to one of them
        for (size_t evtIdx = std::max( 0, Nhitpairs - 20 );
             evtIdx < hitpairs.size();
             evtIdx++) {
            
            SiPMChannelPair & hitpair = hitpairs.at(evtIdx);
            double dt_ms = fabs(time_ms - hitpair.GetEventTime());
            if (verbose>3){
                std::cout    << " Comparing hit to event " << hitpair.GetEventID()
                << " (det. "  << hitpair.GetDetector()
                << ", dt "   << std::setprecision(3) << dt_ms*1.e6 << " ns "
                << ", SiPM " << hitpair.GetSiPM() << ")"
                << std::endl;
            }

            if (    (hitpair.GetDetector() == detector)
                &&  (dt_ms < hitpair.GetPairTimeWindow())
                &&  (hitpair.GetSiPM() == SiPM_number ) 
                ) {
                hitpair.AddHit( time_ms, charge, channel );
                hitpair. CheckIfPairIsGood ();
                hitpair.        SetQprompt ();
                hitpair.         SetQtotal ();
                hitpair.    SetTailToTotal ();

                AddedHitToEvent = true;
                
                if (verbose>3){
                    std::cout    << " Added hit to event " << hitpair.GetEventID()
                    << std::endl;
                    std::cout    << "break;" << std::endl;
                }
                break;
            }
        }
        
        // If there were no events to which this hit was added,
        // we create a new event
        if ( AddedHitToEvent == false ){

            SiPMChannelPair new_hitpair;
            new_hitpair.        SetEventID ( eventID );
            new_hitpair.        SetVerbose ( verbose );
            new_hitpair.          SetSetup ( fsetup );
            new_hitpair.       SetDetector ( detector );
            new_hitpair.            AddHit ( time_ms, charge, channel );
            new_hitpair.           SetSiPM ();
            
            if (verbose>3){
                std::cout    << " Creating a new event of a hit-pair "
                << std::endl;
                new_hitpair.Print();
            }
            
            hitpairs.push_back( new_hitpair );
            eventID++;
        }
        

        if (verbose>2){
            std::cout
            << "Done stepping through "
            << "channel " << channel  << " (detector " << detector << ")" << std::endl ;
            PrintEmptyLine();
        }
    }

    return hitpairs;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<SiPMChannelPair> auxiliary::CollectHitPairsPSDFromGroups(std::vector< std::vector<double> > PETsysData){
    //
    // collect pairs of channel hits from the same SiPM
    // according to table-10 of the FEB/s 2ro data-sheet,
    // and couple them together to get the prompt-charge and the total-charge
    // for the same SiPM hit
    // The pairing is done if the pair of channels saw a hit (crossed-threshold)
    // within 10 ns (1-3 ns variation between channels is expected due to wiring [Luis F., Feb-2022])
    //
    // For simplicity and to keep consistency with existing code above,
    // this routine is implemented in two steps:
    // (1) Collect "detector-events" using auxiliary::CollectEventsFromGroups()
    // (2) Split pairs of channels that read the same SiPM, into SiPMChannelPair objects
    //
    // return:
    // a vector of <SiPMChannelPair> objects to be written into a file
    //
    // Feb-13, 2022
    // NOTE:
    // This routine is not finished, as
    // it was deemed deprated in favour of CollectHitPairsFromSinglesWithPSD()
    //
    std::vector<SiPMChannelPair>    sipm_channel_pair_events;
    std::vector<detectorEvent>      events;
    
    // each row is a hit
    // for each hit, the columns are: 'N(SiPMs)','n(SiPM)','time','charge','channel'
    std::vector<double> t_ms;
    std::vector<int>    channels;
    std::vector<double> Q;
    std::vector<double> Q_prompt;
    std::vector<double> Q_total;
    
    int eventNumber = 0;
    
    for (size_t rowIdx=0; rowIdx<PETsysData.size(); rowIdx++){
        
        int     N       = int(PETsysData.at(rowIdx).at(0));
        int     n       = int(PETsysData.at(rowIdx).at(1));
        // convert to time in ms here since the long numbers
        // (PETsys group data time are given in ps) are too long for the computer to digest
        double  time_ms = double(PETsysData.at(rowIdx).at(2))/1.e9;
        double  charge  = LinearizeQDC( PETsysData.at(rowIdx).at(3)) ;
        int     ch      = PETsysData.at(rowIdx).at(4);
        
        if (verbose>2){
            std::cout
            << "N: "        << N
            << ", n: "      << n
            << ", ch: "     << ch
            << std::setprecision(14)
            << ", time: " << time_ms << " ms"
            << std::endl ;
        }
        
        t_ms.push_back(time_ms);
        channels.push_back(ch);
        Q.push_back(charge);
        
        if (n==N-1) {
            // create detector event
            
            SiPMChannelPair sipm_channel_pair_event( eventNumber, channels, t_ms , Q, N, fsetup, verbose );
            sipm_channel_pair_events.push_back( sipm_channel_pair_event );
            
            detectorEvent event( eventNumber, channels, t_ms , Q, N, fsetup) ;
            if (verbose>2) event.Print();
            events.push_back(event);
            eventNumber++;
            Q.clear();
            t_ms.clear();
            channels.clear();
        }
    }
    return sipm_channel_pair_events;
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
            if ( event.CheckIfEventIsGood() ) {
                separated_events.push_back(event);
                
                if (verbose>2) {
                    std::cout << "Kept (retained) event " << event.GetEventID() << std::endl;
                    event.Print();
                }
            } else {
                if (verbose>2) {
                    std::cout << "Omitted (filtered out) event " << event.GetEventID() << std::endl;
                    event.Print();
                }
            }
        }
        if (verbose>2) {
            std::cout << "done splitting and filtering event " << eventID << std::endl;
            std::cout << "------------------------------------------------------- " << std::endl;
        }
        
    }
    return separated_events;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<SiPMChannelPair> auxiliary::SeparateAndFilterEvents(std::vector<SiPMChannelPair> collected_events){
    
    std::vector<SiPMChannelPair> filtered_events;
    
    for (auto event: collected_events){
        int eventID = event.GetEventID();

        if (verbose>3){
            std::cout    << " SeparateAndFilterEvents() "
            << std::endl << " collected_event #" << eventID
            << std::endl;
            event.Print();
        }

        
        if ( event.CheckIfPairIsGood() ) {
            filtered_events.push_back(event);
            if (verbose>2) {
                std::cout << "*Kept* (retained) event " << event.GetEventID() << std::endl;
                event.Print();
            }
        } else {
            if (verbose>2) {
                std::cout << "*Omitted* (filtered out) event " << event.GetEventID() << std::endl;
                event.Print();
            }
        }
        if (verbose>2) {
            std::cout << "done splitting and filtering event " << eventID << std::endl;
            std::cout << "------------------------------------------------------- " << std::endl;
        }
        
    }
    return filtered_events;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::StreamEventsToCSV (std::vector<detectorEvent> events,
                                   std::string filename,
                                   std::string header){

    OpenEventsCSV(filename, header);
    Debug(0,"writing " + std::to_string(events.size()) + " events");
    for (auto event: events){
        WriteEventToCSV( event );
    }
    close_events_csv();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::StreamEventsToCSV (std::vector<SiPMChannelPair> events,
                                   std::string filename,
                                   std::string header){

    OpenEventsCSV(filename, header);
    Debug(0,"writing " + std::to_string(events.size()) + " events");
    for (auto event: events){
        WriteEventToCSV( event );
    }
    close_events_csv();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::WriteEventToCSV ( detectorEvent event ){
    // "eventID,N(SiPMs),time[ms],Qtot[a.u.],detector,channels"
    csvfile
    << event.GetEventID()      << ","
    << event.GetNhits()        << ","
    << std::setprecision(15) << event.GetEventTime()    << ","
    << event.GetEventCharge()  << ","
    << event.GetDetector()     << ",";
    
    csvfile << "[";
    for (auto ch:event.GetChannels()){
        csvfile << ch << ";";
    }
    csvfile << "]";
    csvfile << std::endl;
    
    // debug...
    //    if (event.GetEventID()==340217 ||
    //        event.GetEventID()==340218 ||
    //        event.GetEventID()==357817 ||
    //        event.GetEventID()==357818 ||
    //        event.GetEventID()==369732 ||
    //        event.GetEventID()==369733 ) {
    //
    //        std::cout << "Found event " << event.GetEventID() << std::endl;
    //        event.Print();
    //
    //    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::WriteEventToCSV ( SiPMChannelPair pair ){
    // "eventID,SiPM,time[ms],Q(prompt)[a.u.],Q(total)[a.u.],tail/total,detector,channels"
    csvfile
    << pair.GetEventID()        << ","
    << pair.GetSiPM()           << ","
    << std::setprecision(15)    << pair.GetEventTime()    << ","
    << pair.GetQprompt()        << ","
    << pair.GetQtotal()         << ","
    << pair.GetTailToTotal()    << ","
    << pair.GetDetector()       << ",";
    
    csvfile << "[";
    for (auto ch:pair.GetChannels()){
        csvfile << ch << ";";
    }
    csvfile << "]";
    csvfile << std::endl;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<double> auxiliary::CollectDetectionTimeDifferencesArray(std::vector<detectorEvent> events,
                                                                    double dt_max_ms){
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
        // 100 events are thus 270 micro-seconds, and 200 events is around 0.5 ms
        for (size_t next_evtIdx = evtIdx+1; next_evtIdx < std::min(evtIdx+200, events.size()); next_evtIdx++){
            double dt = event_times_ms.at(next_evtIdx) - t;
            
            if ((dt==0) && (verbose>2)){
                std::cout
                << "CollectDetectionTimeDifferencesArray() dt=0 between events "
                << evtIdx << " (" << t << " ms)  and " << next_evtIdx << " (" << event_times_ms.at(next_evtIdx) << " ms)" << endl;
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
        std::cout
        << "collected time difference vector of size " << time_differences_ms.size() << " from " << events.size() << " events "
        << std::endl;
    }
    return time_differences_ms;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void auxiliary::ExtractSinglesDoublesTriples(std::vector<detectorEvent> events,
                                             double R_gate,
                                             double RA_gate_delay){
    Debug(0,"ExtractSinglesDoublesTriples()");
    // extract number of singles, doubles, and triples from event data set
    int     Nsingles=0, Ndoubles_direct=0, Ntriples_direct=0;
    double              Ndoubles_multip=0, Ntriples_multip=0;
    
    std::vector<double> t; // event time stamp in ms
    std::vector<double> dt_adjacent; // adjacent time-differences
    
    // (1) form the timestamp array
    for (size_t evtIdx=0; evtIdx < events.size(); evtIdx++) {
        t.push_back( events.at(evtIdx).GetEventTime() );
    }
    // (2) sort timestamp array
    std::sort(t.begin(), t.end());
    
    // (3) form the adjacent time differences array
    for (size_t evtIdx=0; evtIdx < events.size(); evtIdx++) {
        if (evtIdx < events.size()-1){
            dt_adjacent.push_back( t.at(evtIdx+1) - t.at(evtIdx) );
        }
    }
    
    
    // singles is the number of events
    Nsingles = (int)t.size();
    
    // direct counting of doubles and triples, from the adjacent time differences array
    for (size_t evtIdx=0; evtIdx < events.size()-1; evtIdx++) {
        if ( dt_adjacent.at(evtIdx) < R_gate ) {
            Ndoubles_direct ++ ;
            if (verbose>0){
                std::cout
                << "dt(event " << evtIdx                << " - event "<< evtIdx+1 << ") = "
                << dt_adjacent.at(evtIdx)*ms_to_ns      << " ns, "
                << "counting as a double! "
                << "( t(" << evtIdx     << ")=" << t.at(evtIdx)*ms_to_ns   << " ns, "
                << "t(" << evtIdx+1   << ")=" << t.at(evtIdx+1)*ms_to_ns << " ns )"
                << std::endl;
            }
            if ( evtIdx < events.size()-2 ){
                if (( dt_adjacent.at(evtIdx) + dt_adjacent.at(evtIdx+1)) < R_gate ) {
                    Ntriples_direct ++ ;
                }
            }
        }
    }
    
    if (verbose>0){
        std::cout
        << "collected "
        << " timestamp array of size "                   << t.size()
        << ", adjacent time difference vector of size "   << dt_adjacent.size()
        << std::endl
        << Nsingles << " singles, "
        << std::endl
        << "Direct counting: "
        << Ndoubles_direct << " doubles, "
        << Ntriples_direct << " triples, "
        << std::endl;
    }

    // multiplicity arithmatics
    int                 mul_R,      mul_RA;     // multiplicity
    std::vector<int>    mul_ctr_R,  mul_ctr_RA; // multioplicity counter
    // initialize
    int vMax = 50;
    for (int v=0; v<vMax; v++) {
        mul_ctr_R.push_back(0);
        mul_ctr_RA.push_back(0);
    }
    for (size_t evtIdx=0; evtIdx < events.size(); evtIdx++) {
        double t1 = t.at(evtIdx);
        
        std::vector<double> t_to_end;
        t_to_end.assign( t.begin()+evtIdx+1 ,t.end());
        
        mul_R = mul_RA = 0;
        for (auto t2: t_to_end) {
            if ((t1 <= t2)                  && (t2 <= t1 + R_gate )) {
                mul_R ++;
            }
            if ((t1 + RA_gate_delay <= t2)  && (t2 <= t1 + RA_gate_delay + R_gate )) {
                mul_RA ++;
            }
        }
        mul_ctr_R.at( mul_R ) ++;
        mul_ctr_RA.at( mul_RA ) ++;
        
    }
    
    // extract s,d,t from multiplicity
    double sum_Pv=0.,sum_v_Pv=0., sum_v_v1_Pv=0., sum_Qv=0., sum_v_Qv=0., sum_v_v1_Qv=0.;
    for (int v=0; v<vMax; v++) {
        sum_Pv      += mul_ctr_R.at(v);
        sum_v_Pv    += v * mul_ctr_R.at(v);
        sum_v_v1_Pv += v * (v-1) * mul_ctr_R.at(v);
        sum_Qv      += mul_ctr_RA.at(v);
        sum_v_Qv    += v * mul_ctr_RA.at(v);
        sum_v_v1_Qv += v * (v-1) * mul_ctr_RA.at(v);
    }
    Ndoubles_multip = (sum_v_Pv - sum_v_Qv);
    Ntriples_multip = (0.5 * sum_v_v1_Pv - 0.5 * sum_v_v1_Qv - (sum_v_Qv / sum_Qv) *(sum_v_Pv - sum_v_Qv) );

    if (verbose>0){
        std::cout
        << "Using multiplicity arithmetics: "
        << Ndoubles_multip << " doubles, "
        << Ntriples_multip << " triples, "
        << std::endl;
    }

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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double auxiliary::LinearizeQDC ( double Q ) {
    // linearise deposited charge (QDC) in the TOFPET2
    // Charge conversion to energy
    // [PETsys TOFPET2 ASIC Evaluation Kit - Software User Guide v2019.09, p. 25]
    // PETsys Electronics has certified that using
    //$$ E=P_0 \cdot P_1^{Q^{P2}} + P_3 \cdot Q - P_0$$
    // with
    //$$ (P_0, P_1, P_2, P_3) = (8.00000, 1.04676, 1.02734, 0.31909) $$
    //allows to approximate the non-linear response of the ASIC's QDC for all channels,
    //up to a level of 1-2% on the energy resolution of the 511 keV photopeak.
    double P[4] = {8.00000, 1.04676, 1.02734, 0.31909};
    double linearisedQ = P[0]*( pow(P[1], pow(Q,P[2]))) + P[3]*Q - P[0];
    
    return linearisedQ;
}
