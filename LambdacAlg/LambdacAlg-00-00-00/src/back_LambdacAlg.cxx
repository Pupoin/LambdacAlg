//ucas-PRLi (lipeirong11@mails.ucas.ac.cn)
#include "LambdacAlg/LambdacAlg.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "BestDTagSvc/BestDTagSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "DstEvent/TofHitStatus.h"

#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"


#include "McTruth/McParticle.h"
#include "EvTimeEvent/RecEsTime.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

#include "DTagTool/DTagTool.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include "CLHEP/Geometry/Point3D.h"

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include <vector>
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
typedef std::vector<double> Vdouble;
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;

int Ntotal(0),Ncut0(0),Ncut1(0),Ncut2(0),Ncut3(0),Ncut4(0),Ncut5(0),Ncut6(0),Ncut7(0),H(0),A(0);

HepLorentzVector LambdacAlg::getP4(RecEmcShower* gTrk, Hep3Vector origin)
{ 
	Hep3Vector Gm_Vec(gTrk->x(), gTrk->y(), gTrk->z());
	Hep3Vector Gm_Mom = Gm_Vec - origin;
	Gm_Mom.setMag(gTrk->energy());
	HepLorentzVector pGm(Gm_Mom, gTrk->energy());
	return pGm;
}


LambdacAlg::LambdacAlg(const std::string& name, ISvcLocator* pSvcLocator) :Algorithm(name,pSvcLocator){
	declareProperty("Vr0cut",   m_vr0cut=1.0);
	declareProperty("Vz0cut",   m_vz0cut=10.0);
	declareProperty("SkimFlag", m_skim=false);
	declareProperty("PhotonMinEnergy", m_minEnergy = 0.025);
	declareProperty("GammaAngleCut", m_gammaAngleCut=10.0);
        declareProperty("Test1C", m_test1C = 1);                                               //ec fit
	declareProperty("GammathCut",  m_gammathCut=14.0);
	declareProperty("GammatlCut",   m_gammatlCut=0.0);
	declareProperty("PhotonMaxCosThetaBarrel", m_maxCosThetaBarrel = 0.8);
	declareProperty("PhotonMinCosThetaEndcap", m_minCosThetaEndcap = 0.84);
	declareProperty("PhotonMaxCosThetaEndcap", m_maxCosThetaEndcap = 0.92);
	declareProperty("PhotonMinEndcapEnergy",   m_minEndcapEnergy   = 0.050);
	declareProperty("Debug",   m_debug   = false);
	declareProperty("BeamE",  m_beamE = 2.30 );
	declareProperty("ReadBeamEFromDB", m_ReadBeamEFromDB = false );
	declareProperty( "UseCalibBeamE",   m_usecalibBeamE = false );
}
LambdacAlg::~LambdacAlg(){
	//add your code for deconstructor

}
StatusCode LambdacAlg::initialize(){
	if (m_debug)  cout<<__LINE__<<" begin initialize "<<endl;
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"LambdacAlg::initialize()"<<endreq;
	m_beta.setX(0.011);
	m_beta.setY(0);
	m_beta.setZ(0);
	//add your code here
	StatusCode status;
	NTuplePtr nt1(ntupleSvc(), "FILE1/tree");
	if ( nt1 ) m_tuple1 = nt1;
	else {
		m_tuple1 = ntupleSvc()->book ("FILE1/tree", CLID_ColumnWiseTuple, "exam N-Tuple example");
		if ( m_tuple1 )    {
			status = m_tuple1->addItem ("run",   m_run);
			status = m_tuple1->addItem ("event", m_event);
			status = m_tuple1->addItem ("rightflag", m_rightflag);
			status = m_tuple1->addItem("indexmc",          m_idxmc, 0, 100);
			status = m_tuple1->addIndexedItem("pdgid",     m_idxmc, m_pdgid);
			status = m_tuple1->addIndexedItem("motheridx", m_idxmc, m_motheridx);
			status = m_tuple1->addItem ("p4index", m_p4index, 0, 10);
			status = m_tuple1->addIndexedItem ("K", m_p4index, m_K_p4);
			status = m_tuple1->addIndexedItem ("p", m_p4index, m_p_p4);
			status = m_tuple1->addIndexedItem ("pi", m_p4index, m_pi_p4);
			status = m_tuple1->addIndexedItem ("gam1", m_p4index ,m_gam1_p4);
			status = m_tuple1->addIndexedItem ("gam2", m_p4index ,m_gam2_p4);

			status = m_tuple1->addItem ("pi0m", m_pi0m);
			status = m_tuple1->addItem ("pi0m1c", m_pi0m1c);
			status = m_tuple1->addItem ("sigmam", m_Sigmam);
			status = m_tuple1->addItem ("num_othertrackp", m_numothertrackp);
			status = m_tuple1->addItem ("num_othertrackm", m_numothertrackm);
			status = m_tuple1->addItem ("E_beam", m_ebeam);
			status = m_tuple1->addItem ("deltaE_min", m_deltaE_min);
			status = m_tuple1->addItem ("M_BC", m_bc);
			status = m_tuple1->addItem ("chi2",  m_chi2);
			status = m_tuple1->addItem ("ndaughterAp",m_ndaughterAp_, 0, 15);
			status = m_tuple1->addIndexedItem ("Ap_id", m_ndaughterAp_, m_Ap_id_);
			status = m_tuple1->addIndexedItem ("Ap_ptruth", m_ndaughterAp_,  4, m_Ap_ptruth_);
		}
		else    {
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
			return StatusCode::FAILURE;
		}
	}
	//--------end of book--------
   //
 
   log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	return StatusCode::SUCCESS;

}
StatusCode LambdacAlg::beginRun(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"LambdacAlg::beginRun()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;

}
StatusCode LambdacAlg::execute(){
	m_rightflag=-1;
	if (m_debug)  cout<<"m_debug1 begin execute "<<endl;
	Ntotal++;
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"LambdacAlg::execute()"<<endreq;
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");

	m_run = eventHeader->runNumber();
	m_event = eventHeader->eventNumber();

if(!(m_run==35966&&(m_event==10304||m_event==77483||m_event==79947||m_event==97315)))   return StatusCode::SUCCESS;
cout<<endl<<endl<<"********************************************************************************************************"<<endl;
cout<<"run : "<<m_run<<"  event : "<<m_event<<endl;

	log << MSG::DEBUG <<"run, evtnum = " << m_run << " , " << m_event <<endreq;

	IMcDecayModeSvc* i_svc;
	StatusCode sc_DecayModeSvc = service("McDecayModeSvc", i_svc);
	if ( sc_DecayModeSvc.isFailure() )
	{
		log << MSG::FATAL << "Could not load McDecayModeSvc!" << endreq;
		return sc_DecayModeSvc;
	}
	m_svc = dynamic_cast<McDecayModeSvc*>(i_svc);

 int ndaughterAp=0; double Ap_ptruth[15][4]; int Ap_id[15];
 for ( int aa = 0; aa < 15; aa++ )
	    for ( int ll = 0; ll < 4; ll++ )
				  Ap_ptruth[aa][ll]=0;
 for ( int aa = 0; aa < 15; aa++ ) 
	    Ap_id[aa]=0;

	if (eventHeader->runNumber()<0)
	{
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");

		int numParticle = 0;
		if (!mcParticleCol)
		{
			std::cout << "Could not retrieve McParticelCol" << std::endl;
			// return StatusCode::FAILURE;
		}
		else
		{
	    Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
			for (; iter_mc != mcParticleCol->end(); iter_mc++)
			{
				            if ((*iter_mc)->primaryParticle()) continue;
				            if (!(*iter_mc)->decayFromGenerator()) continue;
					          int pdg = (*iter_mc)->particleProperty();
		                int motherpdg = ((*iter_mc)->mother()).particleProperty();
						        int mmotherpdg = (((*iter_mc)->mother()).mother()).particleProperty();
			              if (pdg==111&&motherpdg==3222&&mmotherpdg==4122)
										{
							                     const SmartRefVector<Event::McParticle>& gc = (*iter_mc)->daughterList();
																	 cout<<"gc.size="<<gc.size()<<endl;
									                 for(unsigned int ii = 0; ii < gc.size(); ii++)
																	 {
							                     if( gc[ii]->particleProperty() == -22) continue;
							                     Ap_id[ndaughterAp]=gc[ii]->particleProperty();
							                     for( int ll = 0; ll < 4; ll++ ) Ap_ptruth[ndaughterAp][ll]=gc[ii]->initialFourMomentum()[ll];     
						                       ndaughterAp++;         A++;        
		                               }// End of "gc.size() > 0" IF
				             }
										  m_ndaughterAp_= ndaughterAp;
											       for ( int aa = 0; aa < ndaughterAp; aa++ ) m_Ap_id_[aa]=Ap_id[aa];
											       for ( int aa = 0; aa < ndaughterAp; aa++ )
											       for ( int ll = 0; ll < 4; ll++ ) 
											               m_Ap_ptruth_[aa][ll]=Ap_ptruth[aa][ll];
			}
		}	
	}	
		 
	//end mc
	//TADA!!!=========================================================================
	//good track 
	if (m_debug)  cout<<__LINE__<<" begin choose good track"<<endl;

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);

	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);

	Vint iGoodtotal; iGoodtotal.clear();

	int nCharge = 0;
	Hep3Vector xorigin(0,0,0);
	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid())
	{
		double* dbv = vtxsvc->PrimaryVertex();
		double*  vv = vtxsvc->SigmaPrimaryVertex();
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
	}
	if (m_debug) cout<<__LINE__<<"choose good track"<<endl;
	for(int i = 0; i < evtRecEvent->totalCharged(); i++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();

		double pch=mdcTrk->p();
		//double x0=mdcTrk->x();
		//double y0=mdcTrk->y();
		//double z0=mdcTrk->z();
		//double phi0=mdcTrk->helix(1);
		//double xv=xorigin.x();
		//double yv=xorigin.y();
		//double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);
		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
		VFHelix helixip(point0,a,Ea);
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0=vecipa[1];

		if(fabs(Rvz0) >= m_vz0cut) continue;
		if(fabs(Rvxy0) >= m_vr0cut) continue;

		iGoodtotal.push_back(i);
		nCharge += mdcTrk->charge();
	}// Finish Good Charged Track SKction
       H++;
        int nGoodtotal = iGoodtotal.size();
        Vint iGood; iGood.clear();
	for(int i = 0; i < nGoodtotal; i++)                                              // !!!!!!!!!!!!!!!!!!trying to cut the e- out
        {
                EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGoodtotal[i];
                RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
                double energy = 0.;
                if((*itTrk)->isEmcShowerValid())
                {
                        RecEmcShower *emcTrk = (*itTrk)->emcShower();
                        energy = emcTrk->energy();
                        double m_eop_lm = energy/(mdcTrk->p());
                        if(m_eop_lm>0.8)   continue;
                        iGood.push_back(i);
                }
                else
                {      iGood.push_back(i);
                }
        }
        int nGood = iGood.size();
	if(nGood < 3) return StatusCode::SUCCESS;
	Ncut0++;

	Vint ip,iK,ipi; ip.clear(); iK.clear();ipi.clear();
	int ntrackm(0),ntrackp(0);
	for(int i = 0; i < nGood; i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
			if(mdcTrk->charge()==-1)
			{
				ntrackm++;
				ISimplePIDSvc *m_simplePIDSvc;
				Gaudi::svcLocator()->service("SimplePIDSvc",m_simplePIDSvc);
				m_simplePIDSvc->preparePID(*itTrk);
				if(m_simplePIDSvc ->ispion()) ipi.push_back(iGood[i]);     //xia meng de
                            //    A++;
			}

			if(mdcTrk->charge()==1)
			{
				ntrackp++;
				ISimplePIDSvc *m_simplePIDSvc;
				Gaudi::svcLocator()->service("SimplePIDSvc",m_simplePIDSvc);
				m_simplePIDSvc->preparePID(*itTrk);
				if(m_simplePIDSvc ->isproton()) ip.push_back(iGood[i]);
                              //       B++;
				if(m_simplePIDSvc ->iskaon()) iK.push_back(iGood[i]);        //also xia meng de
                                //       C++;
			}

		}
	int npi = ipi.size(); int np = ip.size();int nK = iK.size();
//	cout<<"npi="<<npi<<"np="<<np<<"nK="<<nK<<nGood<<endl;
	if(npi==0||np==0||nK==0)  return StatusCode::SUCCESS;
	Ncut1++;

	Vp4 ppi, pp,pK; ppi.clear(); pp.clear();pK.clear();
	for(int i = 0; i < npi; i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ipi[i];
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		mdcKalTrk->setPidType(RecMdcKalTrack::pion);                                  //still....
		HepLorentzVector pi4 = mdcKalTrk->p4(xmass[2]);
		ppi.push_back(pi4);
	}
	for(int i = 0; i < np; i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ip[i];
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		mdcKalTrk->setPidType(RecMdcKalTrack::proton);
		HepLorentzVector pp4 = mdcKalTrk->p4(xmass[4]);
		pp.push_back(pp4);
	}
	for(int i = 0; i < nK; i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iK[i];
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		mdcKalTrk->setPidType  (RecMdcKalTrack::kaon);                              //oops!!!!
		HepLorentzVector pK4 = mdcKalTrk->p4(xmass[3]);          
		pK.push_back(pK4);
	}
	//good shower 
	Vint iGam;
	iGam.clear();
	if (m_debug)  cout<<__LINE__<<"begin choose good gamma"<<endl;
	for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++)
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		Hep3Vector emcp(emcTrk->x(), emcTrk->y(), emcTrk->z());
		HepLorentzVector shP4 = getP4(emcTrk,xorigin);
		double cosThetaSh = shP4.vect().cosTheta();
		//double dthe = 200.;
		//double dphi = 200.;
		double dang = 200.;
		if (m_debug)  cout<<__LINE__<<"choose good gamma"<<endl;
		for(int j = 0; j < evtRecEvent->totalCharged(); j++)
		{
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
			if(!(*jtTrk)->isExtTrackValid()) continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();
			if(extTrk->emcVolumeNumber() == -1) continue;
			Hep3Vector extp = extTrk->emcPosition();
			double angd = extp.angle(emcp);
			//double thed = extp.theta() - emcp.theta();
			//double phid = extp.deltaPhi(emcp);
			//thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
			//phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
			if(angd < dang)
			{
				dang = angd;
				//dthe = thed;
				//dphi = phid;
			}
		}
		if(dang>=200) continue;
		double eraw = emcTrk->energy();
		double getTime = emcTrk->time();
		//dthe = dthe * 180 / (CLHEP::pi);
		//dphi = dphi * 180 / (CLHEP::pi);
		dang = dang * 180 / (CLHEP::pi);
		if(fabs(dang) < m_gammaAngleCut) continue;
		if(getTime>m_gammathCut||getTime<m_gammatlCut) continue;
		if(!((fabs(cosThetaSh) < m_maxCosThetaBarrel&&eraw > m_minEnergy)||((fabs(cosThetaSh) > m_minCosThetaEndcap)&&(fabs(cosThetaSh) < m_maxCosThetaEndcap)&&(eraw > m_minEndcapEnergy))))  continue;
		iGam.push_back(i);
	} // Finish Good Photon SKction
	int nGam = iGam.size();
	if(nGam < 2) return StatusCode::SUCCESS;
	Ncut2++;

cout<<"nGam : "<<nGam<<endl;
	Vint igam1,igam2;igam1.clear(),igam2.clear();
	Vp4 pgam1,pgam2,pgam1_1C,pgam2_1C; pgam1.clear(),pgam2.clear(),pgam1_1C.clear(),pgam2_1C.clear();
        Vdouble chi2;chi2.clear();
        KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
	// Loop each gamma pair, check ppi0 mass
	for(int i = 0; i < nGam-1; i++) 
	{
		EvtRecTrackIterator itTrki = evtRecTrkCol->begin() + iGam[i];
		RecEmcShower* emcTrki = (*itTrki)->emcShower();
		Hep3Vector emcpi(emcTrki->x(), emcTrki->y(), emcTrki->z());
		HepLorentzVector ptrki = getP4(emcTrki,xorigin);
		for(int j = i+1; j < nGam; j++) 
		{
			EvtRecTrackIterator itTrkj = evtRecTrkCol->begin() + iGam[j];
			RecEmcShower* emcTrkj = (*itTrkj)->emcShower();
			Hep3Vector emcpj(emcTrkj->x(), emcTrkj->y(), emcTrkj->z());
			HepLorentzVector ptrkj = getP4(emcTrkj,xorigin);

			HepLorentzVector p2gpi = ptrki + ptrkj;
			cout<<i<<"  "<<j<<"  ptrki : "<<ptrki<<"  ptrkj : "<<ptrkj<<"  pi0m : "<<p2gpi.m()<<endl;
			if(p2gpi.m()<0.1||p2gpi.m()>0.16) continue;
                        if(m_test1C==1) 
                   {
                                 kmfit->init();
                                 kmfit->AddTrack(0, 0.0, emcTrki);
                                 kmfit->AddTrack(1, 0.0, emcTrkj);
                                 kmfit->AddResonance(0, 0.135, 0, 1);
                                 bool oksq = kmfit->Fit();
                                 
                                 if(oksq) {
                                         chi2.push_back(kmfit->chisq());
																				 cout<<"chi2 : "<<kmfit->chisq()<<endl;
                                 //      HepLorentzVector pi0p4 = kmfit->pfit(0) + kmfit->pfit(1);
                                         pgam1_1C.push_back(kmfit->pfit(0));
                                         pgam2_1C.push_back(kmfit->pfit(1));
                                         igam1.push_back(iGam[i]);
                                         igam2.push_back(iGam[j]); 
		                         pgam1.push_back(ptrki);
			                 pgam2.push_back(ptrkj);
			cout<<i<<"  "<<j<<"  ptrki ok : "<<ptrki<<"  ptrkj ok : "<<ptrkj<<"  pi0m ok : "<<p2gpi.m()<<endl;
                                 }
                  }          
	       }
	}
	int ngam1 = igam1.size();
	int ngam2 = igam2.size();
	if(ngam1 != ngam2) return StatusCode::SUCCESS;
	Ncut3++;//Ncut3 should equal Ncut2;
	int ngam12 = ngam1;
	if(ngam12==0) return StatusCode::SUCCESS;
	Ncut4++;
	///////////////////////
	//  // get beam energy and beta 
	// ///////////////////////                           

	if(m_ReadBeamEFromDB)
	{             
		if(m_usecalibBeamE) m_readDb.setcalib(true);                      
		m_beamE=m_readDb.getbeamE(m_run,m_beamE);      
		if(m_run>0) m_beta=m_readDb.getbeta();                    
	}
	double ebeam=m_beamE;                             
	double deltaE_min = 9999;
        double chisq = 0;
	int K_index = -1 ;
	int p_index = -1 ;
	int gam1_index = -1 ;
	int gam2_index = -1 ;
	int pi_index = -1;

int a=0;	HepLorentzVector K_p4(0,0,0,0),p_p4(0,0,0,0),gam1_p4(0,0,0,0),gam2_p4(0,0,0,0),pi_p4(0,0,0,0),pgam1_1C4p(0,0,0,0),pgam2_1C4p(0,0,0,0);
	for (int i = 0; i < nK; i++)
	{
		for (int j = 0; j < np; j++)
		{
			for (int k = 0; k < ngam12; k++)
			{
				for(int l = 0; l < npi; l++)
				{
					HepLorentzVector psigma =  pp[j] + pgam1_1C4p[k] + pgam2_1C4p[k];
					if(psigma.m()<1.15||psigma.m()>1.22)continue;
					HepLorentzVector pLambda =  pK[i]+ pp[j] + pgam1_1C4p[k] + pgam2_1C4p[k] + ppi[l];
					pLambda.boost(-m_beta);
					//double deltaE = fabs(pLambda.t() - ebeam);
					double deltaE = pLambda.t() - ebeam;
					if (fabs(deltaE)<fabs(deltaE_min))
					{
a=1;
						deltaE_min=deltaE;
						K_index = iK[i];
						p_index = ip[j];
						gam1_index = igam1[k];
						gam2_index = igam2[k];
						pi_index = ipi[l];
						K_p4 = pK[i];
						p_p4 = pp[j];
						gam1_p4 = pgam1[k];
						gam2_p4 = pgam2[k];
						pi_p4 = ppi[l];
						chisq=chi2[k];
					pgam1_1C4p = pgam1_1C[k];
					pgam2_1C4p = pgam2_1C[k];
					}
				}
			}
		}
	}
//cout<<"a="<<a<<endl;
if(a==0)  return StatusCode::SUCCESS;
	for ( int jj = 0; jj < 4; jj++ ) m_K_p4[jj] = K_p4[jj];
	for ( int jj = 0; jj < 4; jj++ ) m_p_p4[jj] = p_p4[jj];
	for ( int jj = 0; jj < 4; jj++ ) m_pi_p4[jj] = pi_p4[jj];
	for ( int jj = 0; jj < 4; jj++ ) m_gam1_p4[jj] = gam1_p4[jj];
	for ( int jj = 0; jj < 4; jj++ ) m_gam2_p4[jj] = gam2_p4[jj];

	m_p4index =4;
        m_chi2 = chisq;
	m_pi0m1c = (pgam1_1C4p+pgam2_1C4p).m();
	m_pi0m = (gam1_p4+gam2_p4).m();
	m_Sigmam = (pgam1_1C4p+pgam2_1C4p +p_p4).m();
	m_numothertrackp = ntrackp-1;
	m_numothertrackm = ntrackm-1;
	m_ebeam =ebeam;
	m_deltaE_min = deltaE_min;
cout<<"m_pi0m1c : "<<m_pi0m1c<<"  m_chi2 : "<<m_chi2<<endl;
	HepLorentzVector pLambda = K_p4+p_p4+pgam1_1C4p+pgam2_1C4p+pi_p4;
	pLambda.boost(-m_beta);
	double mbc2 = ebeam*ebeam - pLambda.v().mag2();
	m_bc = mbc2 > 0 ? sqrt( mbc2 ) : -10;
	//after give value for all branch 
	m_rightflag=1;
	m_tuple1->write();

	Ncut5++;
	return StatusCode::SUCCESS;

}
StatusCode LambdacAlg::endRun(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"LambdacAlg::endRun()"<<endreq;
	//add your code here
	return StatusCode::SUCCESS;

}
StatusCode LambdacAlg::finalize(){
	MsgStream log(msgSvc(), name());
	log<<MSG::INFO<<"LambdacAlg::finalize()"<<endreq;
	//add your code here
	cout << "attention: if Ncut2!= Ncut3 ,you should check"<<endl;
	cout << "Ntotal  "<<Ntotal<<endl;
	cout << "Ncut0   "<<Ncut0<<endl;
	cout << "Ncut1   "<<Ncut1<<endl;
	cout << "Ncut2   "<<Ncut2<<endl;
	cout << "Ncut3   "<<Ncut3<<endl;
	cout << "Ncut4   "<<Ncut4<<endl;
	cout << "Ncut5   "<<Ncut5<<endl;
	cout << "H   "<<H<<endl;
	 cout<<"A"<<A<<endl;
	return StatusCode::SUCCESS;
}

//add your code here,for other member-functions

