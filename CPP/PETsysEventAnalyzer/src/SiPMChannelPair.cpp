#include <iostream>

#include <detectorEvent.hpp>
#include <SiPMChannelPair.hpp>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SiPMChannelPair::SiPMChannelPair(int eventID,
                                 std::vector<int> channels,
                                 std::vector<double> t_ms, // time in ms
                                 std::vector<double> Q,
                                 int N,
                                 std::string setup,
                                 int verbose){
    
    feventID = eventID;
    for (size_t chIdx=0; chIdx < channels.size(); chIdx++){
        
        fchannels.push_back(channels.at(chIdx)) ;
        ft_ms.push_back(t_ms.at(chIdx)) ;
        fQ.push_back(Q.at(chIdx)) ;
        
        fHits.emplace_back( chIdx, channels.at(chIdx) , t_ms.at(chIdx), Q.at(chIdx) );
    }
    SetSetup            ( setup );
    SetVerbose          ( verbose );
    DetermineWhichSiPM  ( fchannels, setup );
    
    // use detectorEvent class to determie which detector is this
    fdetector = detEv.DetermineWhichDetector( fchannels, setup );
    
    SetEventTime        ();
    fN = N;
    CheckIfPairIsGood();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int SiPMChannelPair::DetermineWhichSiPM(std::vector<int> channels, std::string setup){
    // for each channel of the two channels, determine the SiPM
    // if they match - this is a good "pair" of the same SiPM
    // if they don't match - return "-1" as a flag to disregard this pair
    int SiPM_1 = ChannelNumberToSiPM ( channels.at(0), setup );
    int SiPM_2 = -1;
    std::cout << "channels.size(): " << channels.size() << std::endl;
    if ( channels.size()>1 ){
        SiPM_2 = ChannelNumberToSiPM ( channels.at(1), setup );
    }
    
    if (verbose>3){
        std::cout << "DetermineWhichSiPM() "
        << std::endl
        << "SiPM 1 " << SiPM_1
        << ", SiPM 2 " << SiPM_2
        << std::endl;
    }
    if (SiPM_1==SiPM_2) {
        fSiPM = SiPM_1;
    }
    else {
        fSiPM = -1;
    }
    return fSiPM;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiPMChannelPair::SetSiPM() {
    fSiPM = -1;
    if (fchannels.size()){
        fSiPM = ChannelNumberToSiPM( fchannels.at(0) );
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int SiPMChannelPair::ChannelNumberToSiPM(int channel, std::string setup ) {
    // determine the SiPM number in the ASIC,
    int ASICchannel = channel%64;
    int SiPM_number = -1;
    
    for (auto SiPM_ch1_ch2 : SiPM_ch1_ch2_map){
//        std::cout << "channel "
//        << ASICchannel
//        << " vs. ch. "
//        << SiPM_ch1_ch2.at(1) << " & " << SiPM_ch1_ch2.at(2)
//        << std::endl;
        
        if ( (ASICchannel == SiPM_ch1_ch2.at(1)) || (ASICchannel == SiPM_ch1_ch2.at(2))){
            SiPM_number = SiPM_ch1_ch2.at(0);
        }
    }
    
    if (verbose>10){
        std::cout << "ChannelNumberToSiPM() " << std::endl;
        std::cout << "setup: " << setup << std::endl;
        std::cout << "ASIC channel " << ASICchannel << std::endl;
        
        for (auto SiPM_ch1_ch2 : SiPM_ch1_ch2_map){
            std::cout
            << "SiPM "      << SiPM_ch1_ch2.at(0)
            << ", channel " << SiPM_ch1_ch2.at(1)
            << " and "      << SiPM_ch1_ch2.at(2)
            << std::endl;
        }
        std::cout << "channel " << channel << ", SiPM number: " << SiPM_number << std::endl;
    }
    return SiPM_number;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiPMChannelPair::AddHit ( double hit_time_ms, double hit_charge, int hit_channel  ){
    fchannels.push_back(hit_channel);
    ft_ms.push_back(hit_time_ms); // hit time in ms
    fQ.push_back(hit_charge);
    fN = (int)fchannels.size();
    
    fHits.emplace_back( int(fHits.size()), hit_channel , hit_time_ms, hit_charge );
    SetEventTime();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiPMChannelPair::RemoveHit ( int hitIdx ){
    fN--;
    fchannels.erase(fchannels.begin() + hitIdx);
    ft_ms.erase(ft_ms.begin() + hitIdx);
    fQ.erase(fQ.begin() + hitIdx);
    SetEventTime();
    // std::cout << "detectorEvent::RemoveHit fchannels.size():" << fchannels.size() << ", ft.size():" << ft.size() << std::endl;
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool SiPMChannelPair::CheckIfPairIsGood(int NSiPMs){
    fIsGoodPair = false;
       
    if ((fN == NSiPMs) &&
        (fQtotal>0) &&
        (fdetector!=0)
        ) {
        fIsGoodPair = true;
    }
    return fIsGoodPair;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiPMChannelPair::SetQprompt (){
    // Define the "prompt" charge as the smaller one of the two charge values
    // in the SiPM channel pair.
    //
    // We define this, instead of fixing the
    // "prompt" and "total" charge to specific channels,
    // because we are not sure at this point (Feb-2022)
    // which of the channels is integrated in the small time window
    // and which in the large one
    //
    fQprompt = -9999;
    CheckIfPairIsGood();
    if (fIsGoodPair){
        fQprompt = *(std::min_element(fQ.begin(), fQ.end()));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiPMChannelPair::SetQtotal (){
    fQtotal = 9999;
    CheckIfPairIsGood();
    if (fIsGoodPair){
        fQtotal = *(std::max_element(fQ.begin(), fQ.end()));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiPMChannelPair::SetTailToTotal (){
    fQtail = ( fQtotal - fQprompt );
    fTailToTotal = fQtail / fQtotal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiPMChannelPair::Print(){
    
    std::cout
    << "event: "        << feventID
    << ", N(channels): "<< fN
    << std::setprecision(15)
    << ", event time: " << ftime_ms*1e6 << " ns"
    << ", detector: "   << fdetector
    << ", SiPM: "       << fSiPM
    << ", Q(prompt): "  << fQprompt
    << ", Q(total): "   << fQtotal
    << ", tail/total: "  << fTailToTotal
    
    << ", IsGood: "     << CheckIfPairIsGood()
    << std::endl
    << "Channels: [";
    for (size_t chIdx=0; chIdx < fchannels.size(); chIdx++){
        std::cout << fchannels.at(chIdx) << ",";
    }
    std::cout << "]" << std::endl;
    // check if there are more than 64 hits in this event
    if (fHits.size() > 64) {
        std::cout << fHits.size() << " hits in event " << feventID << std::endl;
    }
    
    std::cout << "Channel hits: ";
    for (size_t chIdx=0; chIdx < fchannels.size(); chIdx++){
        std::cout     << fchannels.at(chIdx)
        << " ( SiPM " << ChannelNumberToSiPM(fchannels.at(chIdx))
        << " , time " << std::setprecision(15)
        << ft_ms.at(chIdx)*1e6 << " ns), ";
    }
    std::cout << std::endl;

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
}

