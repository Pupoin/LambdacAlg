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

// #include "DstEvent/TofHitStatus.h"
// #include "TFile.h"
// #include "TH1F.h"
#include "TMath.h"
// #include "TofCorrection.C"
#include "EvtRecEvent/EvtRecTrack.h"
// #include <fstream>

class MyPid
{
public:
  MyPid(EvtRecTrack *track);
  ~MyPid();

  bool iselectron();
  bool isproton();
  bool ispion();
  bool iskaon();
  double getChi();
  

private:
  EvtRecTrack *track;
  double prob_e, prob_mu, prob_pi, prob_k, prob_p;
  double chi;
};