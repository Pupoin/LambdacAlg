#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "TMath.h"
#include "LambdacAlg/Particle.h"

using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

  Partilce(EvtRecTrack *track){}
  ~Partilce(){}

// ------------------------  get ------------------------

void Partilce::setTrackIndex(int trackIndex)
{
    this->trackIndex = trackIndex;
}

void setLorentzVector(HepLorentzVector p4)
{
    this->p4 = p4;

void Partilce::setCharge(int charge)
{
    this->charge = charge;
}


// ------------------------  get ------------------------

int getTrackIndex()
{
    return trackIndex;
}

HepLorentzVector getLorentzVector()
{
    return p4;
}
int getCharge()
{
    return charge;
}