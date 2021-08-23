#include <iostream>

#include <SiPMHit.hpp>

SiPMHit::SiPMHit(int hitID, int channel, double t_ms, double Q){
    fhitID  =hitID;
    fchannel=channel;
    ft_ms   =t_ms;
    fQ      = Q;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SiPMHit::Print(){
    
    std::cout
    << "hit: "          << fhitID
    << ", time: "       << ft_ms    << " ms"
    << ", channel: "    << fchannel << ","
    << ", charge: "     << fQ       << ",";
    std::cout << std::endl;
    
}
