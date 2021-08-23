#include <iostream>

#include <data_module/detectorEvent.hpp>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
detectorEvent::detectorEvent(int eventID,
                             std::vector<int> channels,
                             std::vector<double> t_ms, // time in ms
                             std::vector<float> Q,
                             int N){
    
    feventID = eventID;
    for (size_t chIdx=0; chIdx < channels.size(); chIdx++){
        fchannels.push_back(channels.at(chIdx)) ;
        ft_ms.push_back(t_ms.at(chIdx)) ;
        fQ.push_back(Q.at(chIdx)) ;
    }
    SetEventTime();
    SetEventCharge();
    fN = N;
    DetermineWhichDetector( fchannels );
    CheckIfEventIsGood();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int detectorEvent::DetermineWhichDetector(std::vector<int> channels){
    int detector = 0;
    int ch_min = *std::min_element(channels.begin(),channels.end());
    int ch_max = *std::max_element(channels.begin(),channels.end());
    if ((64 <=ch_min)  && (ch_max<=127)){
        // 'KETEK 3x3 - 1'
        detector = 12;
    }
    else if ((128 <= ch_min) && (ch_max<=192)){
        // 'KETEK 3x3 - 2'
        detector = 3;
    }
    else if ((895 <= ch_min) && (ch_max<=1023)){
        // 'SensL 6x6 - 1'
        detector = 6;
    }
    else if ((767 <= ch_min) && (ch_max<=895)){
        // 'SensL 6x6 - 2'
        detector = 9;
    }
    else{
        // unknown
        detector = 0;
    }
    
    fdetector = detector;
    return detector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int detectorEvent::ChannelNumberToDetector(int channel){
    int detector = 0;
    if ((64 <=channel)  && (channel<=127)){
        // 'KETEK 3x3 - 1'
        detector = 12;
    }
    else if ((128 <= channel) && (channel<=192)){
        // 'KETEK 3x3 - 2'
        detector = 3;
    }
    else if ((895 <= channel) && (channel<=1023)){
        // 'SensL 6x6 - 1'
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
    return detector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void detectorEvent::AddHit ( double hit_time_ms, float hit_charge, int hit_channel  ){
    fchannels.push_back(hit_channel);
    ft_ms.push_back(hit_time_ms); // hit time in ms
    fQ.push_back(hit_charge);
    fN = (int)fchannels.size();
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
void detectorEvent::Print(){
    
    std::cout
    << "event: "        << feventID
    << ", N(SiPMs): "   << fN       << ", Qtot:" << fQtot
    << ", event time: " << ftime_ms << " ms"
    << ", detector: "   << fdetector
    << ", IsGood: "     << fIsGoodEvent
    << ", channels: [";
    for (size_t chIdx=0; chIdx < fchannels.size(); chIdx++){
        std::cout << fchannels.at(chIdx) << ",";
    }
    std::cout << "]";
    std::cout << std::endl;
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool detectorEvent::CheckIfEventIsGood(int NSiPMs_min){
    fIsGoodEvent = false;
    if ((fN >= NSiPMs_min) && (fQtot>0)) {
        fIsGoodEvent = true;
    }
    return fIsGoodEvent;
}
