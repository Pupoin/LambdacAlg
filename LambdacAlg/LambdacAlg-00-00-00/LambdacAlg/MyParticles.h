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

class MyParticle
{
public:
  MyParticle();
  ~MyParticle();

  void setIndex(int Index);
  void setLorentzVector(HepLorentzVector p4);
  void setCharge(int charge);
  void setTrackParameter(WTrackParameter wtrkp)

  int getIndex();
  HepLorentzVector getLorentzVector();
  int getCharge();
  WTrackParameter getTrackParameter();

private:
    int Index;
    HepLorentzVector p4;
    int charge;
    WTrackParameter wtrkp;
};

class MyMotherParticle:public MyParticle
{
public:
  MyMotherParticle(MyParticle child1, MyParticle child2);
  ~MyMotherParticle();

  // void setChild1(MyParticle child1);
  // void setChild2(MyParticle child2);
  void setChild3(MyParticle child3);

  // MyParticle getChild1();
  // MyParticle getChild2();
  MyParticle getChild3();
private:
  MyParticle child1, child2, child3;
};

  
