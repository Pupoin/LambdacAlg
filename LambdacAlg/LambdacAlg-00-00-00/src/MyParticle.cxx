#include "/LambdacAlg/MyParticle.h"

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
WTrackParameter *MyParticle::getTrackParameter() { return wtrkp }

#pragma endregion

///
///
/// class MyMotherParticle  ///
///
#pragma region for_motherparticle____________________________________
MyMotherParticle::MyMotherParticle() {}
MyMotherParticle::MyMotherParticle(MyParticle child1, HepLorentzVector fitP1, MyParticle child2,
                                   HepLorentzVector fitP2, MyParticle child3, HepLorentzVector fitP3)
{
  this->child1 = child1;
  this->child2 = child2;
  this->child3 = child3;
  this->fitP1 = fitP1;
  this->fitP2 = fitP2;
  this->fitP3 = fitP3;
}
MyMotherParticle::~MyMotherParticle() {}

// void MyMotherParticle::setChild1(MyParticle child1) { this->child1 = child1; }
// void MyMotherParticle::setChild2(MyParticle child2) { this->child2 = child2; }
// void MyMotherParticle::setChild3(MyParticle child3) { this->child3 = child3; }

MyParticle MyMotherParticle::getChild1() { return child1; }
MyParticle MyMotherParticle::getChild2() { return child2; }
MyParticle MyMotherParticle::getChild3() { return child3; }


MyParticle MyMotherParticle::getFitP1() { return fitP1; }
MyParticle MyMotherParticle::getFitP2() { return fitP2; }
MyParticle MyMotherParticle::getFitP3() { return fitP3; }

#pragma endregion