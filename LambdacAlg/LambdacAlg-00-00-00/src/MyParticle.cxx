#include "LambdacAlg/MyParticle.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "TMath.h"
#include "VertexFit/WTrackParameter.h"
#include "VertexFit/Helix.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"

using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

#pragma region for_myparticle________________________________________
///
///
/// class MyParticle  ///
///
MyParticle::MyParticle() {}
MyParticle::MyParticle(int index, HepLorentzVector p4, WTrackParameter *wtrkp, int charge = 0)
{
  this->index = index;
  this->p4 = p4;
  this->wtrkp = wtrkp;
  this->charge = charge;
}
MyParticle::MyParticle(int index, HepLorentzVector p4, RecEmcShower *emcTrk = nullptr)
{
  this->index = index;
  this->p4 = p4;
  this->emcTrk = emcTrk;
}
MyParticle::~MyParticle() {}

// ------  set -------
// void MyParticle::setIndex(int index) { this->index = index; }
// void MyParticle::setCharge(int charge) { this->charge = charge; }
// void MyParticle::setLorentzVector(HepLorentzVector p4) { this->p4 = p4; }
// void MyParticle::setTrackParameter(WTrackParameter *wtrkp) { this->wtrkp = wtrkp; }

// -------  get --------
int MyParticle::getIndex() { return index; }
int MyParticle::getCharge() { return charge; }
HepLorentzVector MyParticle::getLorentzVector() { return p4; }
WTrackParameter *MyParticle::getTrackParameter() { return wtrkp; }
RecEmcShower *MyParticle::getRecEmcShower() { return emcTrk; }

#pragma endregion

///
///
/// class MyMotherParticleFit  ///
///
#pragma region for_motherparticle____________________________________
MyMotherParticleFit::MyMotherParticleFit() {}
MyMotherParticleFit::MyMotherParticleFit(MyParticle child1, MyParticle child2, MyParticle child3 = MyParticle(),
                                         KalmanKinematicFit *kmfit);
{
  this->child1 = child1;
  this->child2 = child2;
  this->child3 = child3;
  this->kmfit = kmfit;
}
MyMotherParticleFit::~MyMotherParticleFit() {}

// void MyMotherParticleFit::setChild1(MyParticle child1) { this->child1 = child1; }
// void MyMotherParticleFit::setChild2(MyParticle child2) { this->child2 = child2; }
// void MyMotherParticleFit::setChild3(MyParticle child3) { this->child3 = child3; }

MyParticle MyMotherParticleFit::getChild1() { return child1; }
MyParticle MyMotherParticleFit::getChild2() { return child2; }
MyParticle MyMotherParticleFit::getChild3() { return child3; }

KalmanKinematicFit *MyMotherParticleFit::getFit() { return kmfit; }

#pragma endregion
