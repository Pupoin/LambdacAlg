#ifndef LambdacAlg_Header
#define LambdacAlg_Header

#include "GaudiKernel/Algorithm.h"
//you can add oher necessary header files
#include "GaudiKernel/NTuple.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "LambdacAlg/ReadBeamInfFromDb.h"
//#include "BestDTagSvc/BestDTagSvc.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

//
//namespace
//
class LambdacAlg:public Algorithm {
	public:
		LambdacAlg(const std::string& name, ISvcLocator* pSvcLocator);
		~LambdacAlg();
		StatusCode initialize();
		StatusCode beginRun();   
		StatusCode execute();
		StatusCode endRun();
		StatusCode finalize();
	//	bool    isGoodks(RecMdcKalTrack* pipTrk, RecMdcKalTrack* pimTrk, double& ks_1chis, double& ks_2chis, double& ks_lchue, HepLorentzVector& p4_ks_1s, double& ks_mass);

	private:
		McDecayModeSvc*  m_svc;
		HepLorentzVector getP4(RecEmcShower* gTrk,Hep3Vector origin);
		double m_costheta;
		double m_vr0cut;
		double m_vz0cut;
		double m_vz0cut1;
		double m_skim;
		double m_minEnergy;
		double m_chisqMax;

		double m_EtaMinMass;
		double m_EtaMaxMass;
		double m_Pi0MinMass;
		double m_Pi0MaxMass;
		double m_SigmaMinMass;
		double m_SigmaMaxMass;

		bool   m_debug;
		bool   m_checktotal;
		double m_maxCosThetaBarrel;
		double m_minCosThetaEndcap;
		double m_maxCosThetaEndcap;
		double m_minEndcapEnergy;
		double m_gammaAngleCut;
		double m_gammathCut;
		double m_gammatlCut;
		double m_beamE;
		Hep3Vector m_beta;
		bool m_ReadBeamEFromDB;
		bool m_usecalibBeamE;
		ReadBeamInfFromDb m_readDb;
		int m_test1C;


		NTuple::Tuple*  m_tuple2;  
		NTuple::Item<int>    m_evtNo_;
		NTuple::Item<int>    m_runNo_;
		NTuple::Item<int>    m_mode1_;
		NTuple::Item<int>    m_mode2_;
		NTuple::Item<int>    m_mode3_;
		NTuple::Item<int>    m_ndaughterAp_;
		NTuple::Array<int>   m_Ap_id_;
		NTuple::Matrix<double>   m_Ap_ptruth_;
		NTuple::Item<int>    m_ndaughterAm_;
		NTuple::Array<int>   m_Am_id_;
		NTuple::Matrix<double>   m_Am_ptruth_;


		NTuple::Tuple* m_tuple1;                    
		NTuple::Item<int> m_run;
		NTuple::Item<int> m_event;
		NTuple::Item<int> m_idxmc;
		NTuple::Item<int> m_all;		
		NTuple::Item<int> m_all_m;
		NTuple::Item<int> m_all_p;
		NTuple::Item<int> m_Kmindex;
		NTuple::Item<int> m_pbarindex;
		NTuple::Item<int> m_pindex;
		NTuple::Item<int> m_Kpindex;
		NTuple::Item<int> m_r;
		NTuple::Item<int> m_rightflag;
		NTuple::Item<int> m_signal;
		NTuple::Item<int> m_bg;
		NTuple::Item<int> m_yes;
		NTuple::Item<int> m_no;
		NTuple::Array<int> m_pdgid;
		NTuple::Array<int> m_motheridx;

		NTuple::Item<int> m_p4index;
		NTuple::Item<int> m_np;
		NTuple::Item<int> m_npbar;
		NTuple::Item<int> m_ngam;
		NTuple::Array<double>m_Kp_p4;
		NTuple::Array<double>m_Km_p4;
		NTuple::Array<double>m_p_p4;
		NTuple::Array<double>m_pbar_p4;
		NTuple::Array<double>m_sigma_p4;
		NTuple::Array<double>m_pim_p4;
		NTuple::Array<double>m_pip_p4;
		NTuple::Array<double>m_gam1a_p4;
		NTuple::Array<double>m_gam2a_p4;
		NTuple::Array<double>m_gam3a_p4;
		NTuple::Array<double>m_gam4a_p4;
		NTuple::Array<double>m_gam1b_p4;
		NTuple::Array<double>m_gam2b_p4;
		NTuple::Array<double>m_gam3b_p4;
		NTuple::Array<double>m_gam4b_p4;
		NTuple::Array<double>m_gama_p4;
		NTuple::Array<double>m_gamb_p4;
		NTuple::Array<double>m_pi0_p4;

		NTuple::Item<double> m_pi0m;
		NTuple::Item<double> m_etam;
		NTuple::Item<double> m_lambda;
		NTuple::Item<double> m_ksi;
		NTuple::Item<double> m_sigmastar;
		NTuple::Item<double> m_etaprimem;
		NTuple::Item<double> m_kstar;
		NTuple::Item<double> m_pi0m1c;
		NTuple::Item<double> m_etam1c;
		NTuple::Item<double> m_chi2;
//		NTuple::Item<double> m_chis1v;
//		NTuple::Item<double> m_chis2v;
//		NTuple::Item<double> m_lchue;
//		NTuple::Item<int> m_ks;
		NTuple::Item<double> m_Sigmam;
		NTuple::Item<double> m_numothertrackp;
		NTuple::Item<double> m_numothertrackm;
		NTuple::Item<double> m_ebeam;
		NTuple::Item<double> m_eop_pip;
		NTuple::Item<double> m_eop_pim;
		NTuple::Item<double> m_eop_p;
		NTuple::Item<double> m_eop_pbar;
		NTuple::Item<double> m_eop_Kp;
		NTuple::Item<double> m_eop_Km;
		NTuple::Item<double> m_deltaE_min;
		NTuple::Item<double> m_bc;
		NTuple::Item<int> m_mcparticle_p;
		NTuple::Array<int> m_pdgid_p;
		NTuple::Array<int> m_motheridx_p;
		NTuple::Item<int> m_mcparticle_m;
		NTuple::Array<int> m_pdgid_m;
		NTuple::Array<int> m_motheridx_m;
		NTuple::Item<int> m_ndaughterAp;
		NTuple::Item<int> m_ndaughterAm;
		NTuple::Array<int> m_Ap_id;
		NTuple::Matrix<double> m_Ap_ptruth;
		NTuple::Array<int> m_Am_id;
		NTuple::Matrix<double> m_Am_ptruth;
		NTuple::Item<double> m_mode1;
		NTuple::Item<double> m_mode2;
		NTuple::Item<double> m_mode3;







	protected:

};
//add your inline methods

//

#endif//LambdacAlg_Header
