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
using CLHEP::HepLorentzVector;

class Partilce
{
public:
  Partilce(EvtRecTrack *track);
  ~Partilce();

    void setTrackIndex(int trackIndex);
    void setLorentzVector(HepLorentzVector p4);
    void setCharge(int charge);

    int getTrackIndex();
    HepLorentzVector getLorentzVector();
    int getCharge();

private:
    int trackIndex;
    HepLorentzVector p4;
    int charge;
};

  
