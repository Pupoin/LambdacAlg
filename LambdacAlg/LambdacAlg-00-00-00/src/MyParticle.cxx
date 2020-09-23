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

using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

#pragma region for_myparticle________________________________________
MyParticle::MyParticle(EvtRecTrack *track) {}
MyParticle::~MyParticle() {}

// ------  set -------
void MyParticle::setIndex(int Index) { this->Index = Index; }
void MyParticle::setCharge(int charge) { this->charge = charge; }
void MyParticle::setLorentzVector(HepLorentzVector p4) { this->p4 = p4; }
void MyParticle::setTrackParameter(WTrackParameter wtrkp) { this->wtrkp = wtrkp; }

// -------  get --------
int MyParticle::getIndex() { return Index; }
int MyParticle::getCharge() { return charge; }
HepLorentzVector MyParticle::getLorentzVector() { return p4; }
WTrackParameter MyParticle::getTrackParameter() { return wtrkp }

#pragma endregion

#pragma region for_motherparticle____________________________________

  MyMotherParticle::MyMotherParticle(MyParticle child1, MyParticle child2)
  {
    this->child1 = child1;
    this->child2 = child2; 
  }
  MyMotherParticle::~MyMotherParticle() { }

  // void MyMotherParticle::setChild1(MyParticle child1) { this->child1 = child1; }
  // void MyMotherParticle::setChild2(MyParticle child2) { this->child2 = child2; }
  void MyMotherParticle::setChild3(MyParticle child3) { this->child3 = child3; }

  // MyParticle MyMotherParticle::getChild1() { return child1; }
  // MyParticle MyMotherParticle::getChild2() { return child2; }
  MyParticle MyMotherParticle::getChild3() { return child3; }

#pragma endregion