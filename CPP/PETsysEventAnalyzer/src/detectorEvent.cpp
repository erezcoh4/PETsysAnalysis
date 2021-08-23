#include <iostream>

#include <detectorEvent.hpp>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
detectorEvent::detectorEvent(int eventID,
                             std::vector<int> channels,
                             std::vector<double> t_ms, // time in ms
                             std::vector<double> Q,
                             int N,
                             std::string setup){
    
    feventID = eventID;
    for (size_t chIdx=0; chIdx < channels.size(); chIdx++){
        
        fchannels.push_back(channels.at(chIdx)) ;
        ft_ms.push_back(t_ms.at(chIdx)) ;
        fQ.push_back(Q.at(chIdx)) ;
        
        fHits.emplace_back( chIdx, channels.at(chIdx) , t_ms.at(chIdx), Q.at(chIdx) );
    }
    SetSetup(setup);
    SetEventTime();
    SetEventCharge();
    fN = N;
    DetermineWhichDetector( fchannels, setup );
    CheckIfEventIsGood();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int detectorEvent::DetermineWhichDetector(std::vector<int> channels, std::string setup){
    int ch_min = *std::min_element(channels.begin(),channels.end());
    int ch_max = *std::max_element(channels.begin(),channels.end());
    
    int detector_ch_min = ChannelNumberToDetector( ch_min, setup );
    int detector_ch_max = ChannelNumberToDetector( ch_max, setup );
    if (detector_ch_min != detector_ch_max){
        fdetector = 0;
    }
    else{
        fdetector = detector_ch_min;
    }
    return fdetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int detectorEvent::ChannelNumberToDetector(int channel, std::string setup) {
    // supported setups:
    //
    // BoxSi_proto2.1
    // PlasticScintillator_SensL8x8
    // BoxSi_proto2.2
    //
    int detector = 0;
    
    if (setup.compare("BoxSi_proto2.1")==0){ // Feb-17, 2021
        if ((0 <=channel)  && (channel<=127)){
            // 'KETEK 3x3 in port 1'
            detector = 12;
        }
        else if ((128 <= channel) && (channel<=255)){
            // 'KETEK 3x3 in port 2'
            detector = 3;
        }
        else if ((895 <= channel) && (channel<=1023)){
            // 'SensL 6x6 in port 8'
            detector = 6;
        }
        else if ((767 <= channel) && (channel<=895)){
            // 'SensL 6x6 - 2'
            detector = 9;
        }
        else{
            // unknown
            detector = 0;
        }
    }
    
    if (setup.compare("PlasticScintillator_SensL8x8")==0){ // Feb-25, 2021
        if ((895 <= channel) && (channel<=1023)){
            // 'SensL 6x6 - 1' coupled to scintillator
            detector = 1;
        }
        else if ((0 <= channel) && (channel<=127)){
            // 'SensL 6x6 - 2' coupled to nothing
            detector = -1;
        }
    }
    
    if (setup.compare("BoxSi_proto2.2")==0){ // Aug-9, 2021
        if ((512 <=channel)  && (channel <= 639)){
            // SensL-1 in port 5, side 12
            detector = 12;
        }
        else if ((640 <= channel) && (channel <= 767)){
            // KETEK-1 in port 6, side 9
            detector = 9;
        }
        else if ((768 <= channel) && (channel <= 895)){
            // SensL-2 in port 7, side 6
            detector = 6;
        }
        else if ((896 <= channel) && (channel <= 1023)){
            // KETEK-1 in port 8, side 3
            detector = 3;
        }
        else{
            // unknown
            detector = 0;
        }
    }
    return detector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void detectorEvent::AddHit ( double hit_time_ms, double hit_charge, int hit_channel  ){
    fchannels.push_back(hit_channel);
    ft_ms.push_back(hit_time_ms); // hit time in ms
    fQ.push_back(hit_charge);
    fN = (int)fchannels.size();
    
    fHits.emplace_back( int(fHits.size()), hit_channel , hit_time_ms, hit_charge );
    
    SetEventCharge();
    SetEventTime();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void detectorEvent::RemoveHit ( int hitIdx ){
    fN--;
    fchannels.erase(fchannels.begin() + hitIdx);
    ft_ms.erase(ft_ms.begin() + hitIdx);
    fQ.erase(fQ.begin() + hitIdx);
    SetEventCharge();
    SetEventTime();
    // std::cout << "detectorEvent::RemoveHit fchannels.size():" << fchannels.size() << ", ft.size():" << ft.size() << std::endl;
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool detectorEvent::CheckIfEventIsGood(int NSiPMs_min){
    fIsGoodEvent = false;
    if ((fN >= NSiPMs_min) && (fQtot>0)) {
        fIsGoodEvent = true;
    }
    return fIsGoodEvent;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void detectorEvent::Print(){
    
    std::cout
    << "event: "        << feventID
    << ", N(SiPMs): "   << fN       << ", Qtot: " << fQtot
    << std::setprecision(10)
    << ", event time: " << ftime_ms << " ms"
    << ", detector: "   << fdetector
    << ", IsGood: "     << CheckIfEventIsGood()
    << std::endl << ", channels: [";
    for (size_t chIdx=0; chIdx < fchannels.size(); chIdx++){
        std::cout << fchannels.at(chIdx) << ",";
    }
    std::cout << "]" << std::endl;
    // check if there are more than 64 hits in this event
    if (fHits.size() > 64) {
        std::cout << fHits.size() << " hits in event " << feventID << std::endl;
    }

    // check if there is more than a single hit from each channel
    for (auto hit1:fHits) {
        for (auto hit2:fHits){
            if ((hit1.GetHitID() != hit2.GetHitID()) && (hit1.GetChannel() == hit2.GetChannel())) {
                std::cout
                << "hit "       << hit1.GetHitID() << " (time= " << hit1.GetTime() << " ms) "
                << " and hit "  << hit2.GetHitID() << " (time= " << hit2.GetTime() << " ms) "
                << " are from the same channel (" << hit1.GetChannel() << ")" << std::endl;
            }
            
        }
    }
    std::cout << std::endl;
}
