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

///
///
/// class MyParticle  ///
///
class MyParticle
{
public:
  MyParticle();
  MyParticle(int index, HepLorentzVector p4, WTrackParameter *wtrkp = nullptr, int charge = 0);
  ~MyParticle();

  // void setIndex(int index);
  // void setCharge(int charge);
  // void setLorentzVector(HepLorentzVector p4);
  // void setTrackParameter(WTrackParameter *wtrkp);

  int getIndex();
  int getCharge();
  HepLorentzVector getLorentzVector();
  WTrackParameter *getTrackParameter();

private:
  int index, charge;
  HepLorentzVector p4;
  WTrackParameter *wtrkp;
};

///
///
/// class MyMotherParticle  ///
///
class MyMotherParticle : public MyParticle
{
public:
  MyMotherParticle();
  MyMotherParticle::MyMotherParticle(MyParticle child1, HepLorentzVector fitP1, MyParticle child2, HepLorentzVector fitP2,
                                     MyParticle child3 = MyParticle(), HepLorentzVector fitP3 = HepLorentzVector());
  ~MyMotherParticle();

  // void setChild1(MyParticle child1);
  // void setChild2(MyParticle child2);
  // void setChild3(MyParticle child3);

  MyParticle getChild1();
  MyParticle getChild2();
  MyParticle getChild3();

  HepLorentzVector getFitP1();
  HepLorentzVector getFitP2();
  HepLorentzVector getFitP3();

private:
  MyParticle child1, child2, child3;
  HepLorentzVector fitP1, fitP2, fitP3;
};
