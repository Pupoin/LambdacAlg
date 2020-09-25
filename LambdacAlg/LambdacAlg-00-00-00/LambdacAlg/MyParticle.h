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
  MyParticle(int index, HepLorentzVector p4, RecEmcShower *emcTrk = nullptr);
  ~MyParticle();

  // void setIndex(int index);
  // void setCharge(int charge);
  // void setLorentzVector(HepLorentzVector p4);
  // void setTrackParameter(WTrackParameter *wtrkp);

  int getIndex();
  int getCharge();
  HepLorentzVector getLorentzVector();
  WTrackParameter *getTrackParameter();
  RecEmcShower *getRecEmcShower();

private:
  int index, charge;
  HepLorentzVector p4;
  WTrackParameter *wtrkp;
  RecEmcShower *emcTrk;
};

///
///
/// class MyMotherParticleFit  ///
///
class MyMotherParticleFit : public MyParticle
{
public:
  MyMotherParticleFit();
  MyMotherParticleFit::MyMotherParticleFit(MyParticle child1, MyParticle child2,
                                           MyParticle child3 = MyParticle(),
                                           KalmanKinematicFit *kmfit=nullptr);
  ~MyMotherParticleFit();

  // void setChild1(MyParticle child1);
  // void setChild2(MyParticle child2);
  // void setChild3(MyParticle child3);

  MyParticle getChild1();
  MyParticle getChild2();
  MyParticle getChild3();

  KalmanKinematicFit *getFit();

private:
  MyParticle child1, child2, child3;
  KalmanKinematicFit *kmfit;
};
