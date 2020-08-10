#ifndef LambdacSigmapEta_Header
#define LambdacSigmapEta_Header

#include "GaudiKernel/Algorithm.h"
//you can add oher necessary header files
#include "GaudiKernel/NTuple.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
//#include "BestDTagSvc/BestDTagSvc.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "SimplePIDSvc/SimplePIDSvc.h"
#include "SimplePIDSvc/ISimplePIDSvc.h"
//
//namespace
//
class LambdacSigmapEta:public Algorithm {
	public:
		LambdacSigmapEta(const std::string& name, ISvcLocator* pSvcLocator);
		~LambdacSigmapEta();
		StatusCode initialize();
		StatusCode beginRun();   
		StatusCode execute();
		StatusCode endRun();
		StatusCode finalize();

	private:
		McDecayModeSvc*  m_svc;
		HepLorentzVector getP4(RecEmcShower* gTrk,Hep3Vector origin);
		double m_pvr0cut;
		double m_pvz0cut;
		double m_costheta;
		double m_minEnergy;
		double m_gammathCut;
		double m_gammatlCut;
		double m_maxCosThetaBarrel;
		double m_minCosThetaEndcap;
		double m_maxCosThetaEndcap;
		double m_minEndcapEnergy;
		double m_SigmaminMass;
		double m_SigmamaxMass;
        double m_SigmaMassFit;
		double m_Pi0minMass;
		double m_Pi0maxMass;
        double m_Pi0MassFit;
		double m_EtaminMass;
		double m_EtamaxMass;
        double m_EtaMassFit;
		bool   m_debug;
		bool   m_getmcinfo;
		double m_beamE;
		bool m_use_Total_TOF;
		double m_chisqMax;
		Hep3Vector m_beta;

		NTuple::Tuple* m_tuple1;
		NTuple::Item<int> m_run;

		NTuple::Item<double> m_pi0;
		NTuple::Item<double> m_eta;
		NTuple::Item<double> m_sigmap;

		NTuple::Item<int> m_event;
		NTuple::Item<int> m_rightflag;
		NTuple::Item<int> m_p4index;

//For Lc
		NTuple::Item<int> m_isig;
		NTuple::Array<int> m_charge;
		NTuple::Array<double> m_chisq_3C;
		NTuple::Array<double> m_chisq_4C;
		NTuple::Array<double> m_angle;
		NTuple::Array<double> m_deltaE_min;
   		NTuple::Array<double> m_bc;
   		NTuple::Array<double> m_sigma;
   		NTuple::Array<double> m_sigma2;

		NTuple::Matrix<double>m_Gam1_Pi0_p4_old;
		NTuple::Matrix<double>m_Gam2_Pi0_p4_old;
		NTuple::Matrix<double>m_Gam1_Pi0_p4;
		NTuple::Matrix<double>m_Gam2_Pi0_p4;

		NTuple::Matrix<double>m_Gam1_Eta_p4_old;
		NTuple::Matrix<double>m_Gam2_Eta_p4_old;
		NTuple::Matrix<double>m_Gam1_Eta_p4;
		NTuple::Matrix<double>m_Gam2_Eta_p4;

		NTuple::Matrix<double>m_gam1_pi0_Sigma_p4;
		NTuple::Matrix<double>m_gam2_pi0_Sigma_p4;
		NTuple::Matrix<double>m_p_Sigma_p4;

		NTuple::Matrix<double>m_gam1_pi0_Sigma_p4_old;
		NTuple::Matrix<double>m_gam2_pi0_Sigma_p4_old;
		NTuple::Matrix<double>m_p_Sigma_p4_old;
		NTuple::Matrix<double>m_gam1_pi0_Sigma_p4_3C;
		NTuple::Matrix<double>m_gam2_pi0_Sigma_p4_3C;
//For mc 
		NTuple::Item<int> m_mode1;
		NTuple::Item<int> m_mode2;
		NTuple::Item<int> m_mode3;

		NTuple::Item<int> m_Lmdc_P;
		NTuple::Item<int> m_Lmdc_M;

		NTuple::Item<int> m_idxmc;
		NTuple::Array<int> m_pdgid;
		NTuple::Array<int> m_motheridx;

		NTuple::Item<int> m_mcparticle_p;
		NTuple::Array<int> m_pdgid_p;
		NTuple::Array<int> m_motheridx_p;

		NTuple::Item<int> m_mcparticle_m;
		NTuple::Array<int> m_pdgid_m;
		NTuple::Array<int> m_motheridx_m;

		NTuple::Item<int> m_ndaughterAp;
		NTuple::Array<int> m_Ap_id;
		NTuple::Matrix<double> m_Ap_ptruth;

		NTuple::Item<int> m_ndaughterAm;
		NTuple::Array<int> m_Am_id;
		NTuple::Matrix<double> m_Am_ptruth;


		NTuple::Tuple* m_tuple2;                    
		NTuple::Item<int> m_run_t;
		NTuple::Item<int> m_event_t;

		NTuple::Item<int> m_mode1_t;
		NTuple::Item<int> m_mode2_t;
		NTuple::Item<int> m_mode3_t;

		NTuple::Item<int> m_Lmdc_P_t;
		NTuple::Item<int> m_Lmdc_M_t;
		NTuple::Item<int> m_Flag_t;

		NTuple::Item<int> m_ndaughterAp_t;
		NTuple::Array<int> m_Ap_id_t;
		NTuple::Matrix<double> m_Ap_ptruth_t;

		NTuple::Item<int> m_ndaughterAm_t;
		NTuple::Array<int> m_Am_id_t;
		NTuple::Matrix<double> m_Am_ptruth_t;


	protected:

};
//add your inline methods

//

#endif//LambdacSigmapEta_Header
