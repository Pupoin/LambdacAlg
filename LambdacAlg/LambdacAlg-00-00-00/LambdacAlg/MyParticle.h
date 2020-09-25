#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "TMath.h"

#include "EvTimeEvent/RecEsTime.h"
#include "McTruth/McParticle.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

#include "CLHEP/Geometry/Point3D.h"
#include "DTagTool/DTagTool.h"

#include "VertexFit/Helix.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"

///
///
/// class MyParticle  ///
///
class MyParticle
{
public:
  MyParticle();
  MyParticle(int index, HepLorentzVector p4, int charge = 0, WTrackParameter wtrkp = WTrackParameter());
  MyParticle(int index, HepLorentzVector p4, RecEmcShower *emcTrk = nullptr);
  ~MyParticle();

  // void setIndex(int index);
  void setCharge(int charge);
  void setLorentzVector(HepLorentzVector p4);
  // void setTrackParameter(WTrackParameter *wtrkp);

  int getIndex();
  int getCharge();
  HepLorentzVector getLorentzVector();
  WTrackParameter *getTrackParameter();
  RecEmcShower *getRecEmcShower();
  double getMass();

private:
  int index, charge;
  HepLorentzVector p4;
  WTrackParameter wtrkp;
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
  MyMotherParticleFit(MyParticle child1, MyParticle child2,
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
