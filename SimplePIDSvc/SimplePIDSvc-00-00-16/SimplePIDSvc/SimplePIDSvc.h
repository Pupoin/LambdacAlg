#ifndef SIMPLE_PID_SVC_H
#define SIMPLE_PID_SVC_H

#include "GaudiKernel/Service.h"
#include "SimplePIDSvc/ISimplePIDSvc.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "TH1F.h"
#include "ParticleID/ParticleID.h"
 
class IDataProviderSvc;
template <class TYPE> class CnvFactory;
 
class SimplePIDSvc : public Service, virtual public ISimplePIDSvc
{
  friend class CnvFactory<SimplePIDSvc>;
 
  public :
    
    SimplePIDSvc(const std::string& name, ISvcLocator* svcLoc);
    virtual ~SimplePIDSvc();
  
    virtual StatusCode initialize();
    virtual StatusCode finalize();
    virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvIF);


    inline void    setdedxminchi(double x){m_dedxminchi=x;}
    inline void    settofminchi(double x){m_tofminchi=x;}
    void    preparePID(EvtRecTrack* track);
    bool    iselectron();
    bool    isproton();
    bool    ispion();
    bool    iskaon();
    bool    ispion_for_Lambda_Ks();
 

    inline double  probElectron() { return m_prob[0];};
    inline double  probMuon() { return m_prob[1];};
    inline double  probPion() { return m_prob[2];};
    inline double  probKaon() { return m_prob[3];};
    inline double  probProton() { return m_prob[4];};

    inline double  probelectron() { return mm_prob[0];};
    inline double  probmuon() { return mm_prob[1];};
    inline double  probpion() { return mm_prob[2];};
    inline double  probkaon() { return mm_prob[3];};
    inline double  probproton() { return mm_prob[4];};

    inline double  getdEdxChi(int i){return m_dedxchi[i];}
    inline double  getTOFChi(int i){return m_tofchi[i];}
    inline double  getChi2(int i){return pow(m_dedxchi[i],2)+pow(m_tofchi[i],2);}
    
    
    
    //EvtRecDTag* getBestDTag(int modeid = -1, int charm = 0);
    
 private:

    //
    bool m_tofcorrec;
    double m_dedxminchi;
    double m_tofminchi;
    int m_run;
    vector<double> m_datatof;
    vector<double> m_mctof;
    //
    double m_eop;
    
    //
        
    bool m_dedxonly[5];
    double m_p[5];
    double m_costh[5];
    double m_mass[5];
    double m_tofscale1[5];
    double m_tofscale2[5];

    double m_dedxchi[5];
    double m_tofchi[5];
    double m_tofdt[5];
    double m_tofdt1[5];
    double m_tofdt2[5];
    double m_sigma1;
    double m_sigma2;
    

    double m_prob[5];
    double mm_prob[5];
    double mmm_prob[5];
  private :
      
    IDataProviderSvc*  eventSvc_;
     
    TH1F* h_ebarlp;
    TH1F* h_ebarhp;
    TH1F* h_eendlp;
    TH1F* h_eendhp;
    TH1F* h_kbar;
    TH1F* h_kend;
    TH1F* h_pibar;
    TH1F* h_piend;
    
};



#endif
