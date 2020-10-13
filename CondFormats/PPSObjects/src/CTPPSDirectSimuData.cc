#include "CondFormats/PPSObjects/interface/CTPPSDirectSimuData.h"
#include <iostream>

// Constructors
CTPPSDirectSimuData::CTPPSDirectSimuData():
    useEmpiricalApertures(0),
    empiricalAperture45(""),
    empiricalAperture56(""),
    timeResolutionDiamonds45(""),
    timeResolutionDiamonds56(""){}

// Destructor
CTPPSDirectSimuData::~CTPPSDirectSimuData() {}

// Getters
bool CTPPSDirectSimuData::getUseEmpiricalApertures() const{return useEmpiricalApertures;}
const std::string& CTPPSDirectSimuData::getEmpiricalAperture45() const{return empiricalAperture45;}
const std::string& CTPPSDirectSimuData::getEmpiricalAperture56() const{return empiricalAperture56;}
const std::string& CTPPSDirectSimuData::getTimeResolutionDiamonds45() const{return timeResolutionDiamonds45;}
const std::string& CTPPSDirectSimuData::getTimeResolutionDiamonds56() const{return timeResolutionDiamonds56;}
bool CTPPSDirectSimuData::getUseTimeEfficiencyCheck() const{return useTimeEfficiencyCheck;}
const std::string& CTPPSDirectSimuData::getEffTimePath() const{return effTimePath;}
const std::string& CTPPSDirectSimuData::getEffTimeObject45() const{return effTimeObject45;}
const std::string& CTPPSDirectSimuData::getEffTimeObject56() const{return effTimeObject56;}


// Setters
void CTPPSDirectSimuData::setUseEmpiricalApertures(bool b){useEmpiricalApertures=b;}
void CTPPSDirectSimuData::setEmpiricalAperture45(std::string s){empiricalAperture45=s;}
void CTPPSDirectSimuData::setEmpiricalAperture56(std::string s){empiricalAperture56=s;}
void CTPPSDirectSimuData::setTimeResolutionDiamonds45(std::string s){timeResolutionDiamonds45=s;}
void CTPPSDirectSimuData::setTimeResolutionDiamonds56(std::string s){timeResolutionDiamonds56=s;}
void CTPPSDirectSimuData::setUseTimeEfficiencyCheck(bool b){useTimeEfficiencyCheck=b;}
void CTPPSDirectSimuData::setEffTimePath(std::string s){effTimePath=s;}
void CTPPSDirectSimuData::setEffTimeObject45(std::string s){effTimeObject45=s;}
void CTPPSDirectSimuData::setEffTimeObject56(std::string s){effTimeObject56=s;}



void CTPPSDirectSimuData::printInfo(std::stringstream& s) {
    s << "\n   useEmpiricalApertures = " << useEmpiricalApertures
    << "\n   empiricalAperture45 = " << empiricalAperture45
    << "\n   empiricalAperture56 = " << empiricalAperture56
    << "\n   timeResolutionDiamonds45 = " << timeResolutionDiamonds45
    << "\n   timeResolutionDiamonds56 = " << timeResolutionDiamonds56
    << "\n useTimeEfficiencyCheck= "<< useTimeEfficiencyCheck
    << "\n effTimePath= "<< effTimePath
    << "\n effTimeObject45= "<< effTimeObject45
    << "\n effTimeObject56= "<< effTimeObject56
    << std::endl;
}

std::ostream& operator<<(std::ostream& os, CTPPSDirectSimuData info) {
  std::stringstream ss;
  info.printInfo(ss);
  os << ss.str();
  return os;
}