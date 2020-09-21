#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
//#include "EvtRecEvent/EvtRecDTag.h"
//#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
// #include "SimplePIDSvc/SimplePIDSvc.h"

#include "DstEvent/TofHitStatus.h"
// #include "TFile.h"
// #include "TH1F.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "LambdacAlg/MyPid.h"
#include "TMath.h"

MyPid::MyPid(EvtRecTrack *track) { this->track = track; }
MyPid::~MyPid() {}

bool MyPid::iselectron()
{
  // pid for electron
  ParticleID *pid = ParticleID::instance();
  pid->init();
  pid->setMethod(pid->methodProbability());
  pid->setChiMinCut(4);
  pid->setRecTrack(track);
  pid->usePidSys(pid->useDedx() | pid->useTofCorr() | pid->useEmc()); // use pid sub-system
  pid->identify(pid->onlyElectron() | pid->onlyPion() | pid->onlyKaon());
  pid->calculate();
  prob_e = pid->probElectron();
  prob_mu = pid->probMuon();
  prob_pi = pid->probPion();
  prob_k = pid->probKaon();
  prob_p = pid->probProton();
  if (!(pid->IsPidInfoValid()))
  {
    return false;
  }

  if (prob_e > 0.001 && (prob_e / (prob_e + prob_pi + prob_k)) > 0.8)
  {
    return true;
  }
  return false;
}
bool MyPid::isproton()
{
  // pid for proton
  ParticleID *pid = ParticleID::instance();
  pid->init();
  pid->setMethod(pid->methodProbability());
  pid->setChiMinCut(8);
  pid->setRecTrack(track);
  pid->usePidSys(pid->useDedx() | pid->useTofCorr()); // use pid sub-system
  pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());
  pid->calculate();
  prob_e = pid->probElectron();
  prob_mu = pid->probMuon();
  prob_pi = pid->probPion();
  prob_k = pid->probKaon();
  prob_p = pid->probProton();
  if (!(pid->IsPidInfoValid()))
  {
    return false;
  }
  if (prob_p >= 0. && prob_p >= prob_k && prob_p >= prob_pi)
  {
    return true;
  }
  return false;
}

bool MyPid::ispion()
{
  // pid for kaon and pion
  ParticleID *pid = ParticleID::instance();
  pid->init();
  pid->setMethod(pid->methodProbability());
  pid->setChiMinCut(4);
  pid->setRecTrack(track);
  pid->usePidSys(pid->useDedx() | pid->useTofCorr()); // use pid sub-system
  pid->identify(pid->onlyPion() | pid->onlyKaon());
  pid->calculate();
  prob_e = pid->probElectron();
  prob_mu = pid->probMuon();
  prob_pi = pid->probPion();
  prob_k = pid->probKaon();
  prob_p = pid->probProton();
  if (!(pid->IsPidInfoValid()))
  {
    return false;
  }
  if (prob_pi >= 0.00 && prob_pi >= prob_k)
  {
    return true;
  }
  return false;
}

bool MyPid::iskaon()
{
  // pid for kaon and pion
  ParticleID *pid = ParticleID::instance();
  pid->init();
  pid->setMethod(pid->methodProbability());
  pid->setChiMinCut(4);
  pid->setRecTrack(track);
  pid->usePidSys(pid->useDedx() | pid->useTofCorr()); // use pid sub-system
  pid->identify(pid->onlyPion() | pid->onlyKaon());
  pid->calculate();
  prob_e = pid->probElectron();
  prob_mu = pid->probMuon();
  prob_pi = pid->probPion();
  prob_k = pid->probKaon();
  prob_p = pid->probProton();
  if (!(pid->IsPidInfoValid()))
  {
    return false;
  }
  if (prob_k >= 0.00 && prob_k >= prob_pi)
  {
    return true;
  }
  return false;
}

// bool SimplePIDSvc::ispion_for_Lambda_Ks()
// {

//   if (prob_pi >= 0.00 && prob_pi >= prob_k)
//     return true;
//   return false;
// }
