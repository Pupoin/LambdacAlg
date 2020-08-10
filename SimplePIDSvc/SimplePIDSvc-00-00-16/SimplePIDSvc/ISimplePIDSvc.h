#ifndef ISIMPLE_PID_SVC_H
#define ISIMPLE_PID_SVC_H

#include "GaudiKernel/IService.h"
#include "EvtRecEvent/EvtRecTrack.h"

/* Decaration of the interface ID */
static const InterfaceID IID_ISimplePIDSvc("ISimplePIDSvc", 1, 0);

class EvtRecDTag;

class ISimplePIDSvc : virtual public IService
{
  public :
    virtual ~ISimplePIDSvc() {}
  
    static const InterfaceID& interfaceID() { return IID_ISimplePIDSvc; }
      
    virtual void setdedxminchi(double x) = 0;
    virtual void settofminchi(double x) = 0;
    virtual void  preparePID(EvtRecTrack* track) = 0;
    //virtual bool  iselectron(bool eop=false) = 0;
    virtual bool  iselectron() = 0;
    virtual bool  isproton() = 0;
    virtual bool  ispion() = 0;
    virtual bool  ispion_for_Lambda_Ks() = 0;
    virtual bool  iskaon() = 0;
    virtual double probElectron() = 0;
    virtual double probMuon() = 0;
    virtual double probPion() = 0;
    virtual double probKaon() = 0;
    virtual double probProton() = 0;
    virtual double probelectron() = 0;
    virtual double probmuon() = 0;
    virtual double probpion() = 0;
    virtual double probkaon() = 0;
    virtual double probproton() = 0;
    virtual double getdEdxChi(int i) = 0;
    virtual double getTOFChi(int i) = 0;
    virtual double getChi2(int i) = 0;
    
};

#endif
