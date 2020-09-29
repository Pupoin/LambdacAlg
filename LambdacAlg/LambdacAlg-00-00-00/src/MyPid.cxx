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

double CalChi(ParticleID *pid)
{
  double tmp = 0;
  for (int i = 0; i < 5; ++i)
  {
    // prob_dEdx_TOF_Kpip[n_charged][i] = pid->prob(i);
    tmp += pid->chi(i);
  }
  return tmp;
}

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
    chi = CalChi(pid);
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
    chi = CalChi(pid);
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
    chi = CalChi(pid);
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
    chi = CalChi(pid);
    return true;
  }
  return false;
}

double MyPid::getChi() { return this->chi; }

// bool SimplePIDSvc::ispion_for_Lambda_Ks()
// {

//   if (prob_pi >= 0.00 && prob_pi >= prob_k)
//     return true;
//   return false;
// }
