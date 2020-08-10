#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "EventModel/EventModel.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EventModel/EventHeader.h"
//#include "EvtRecEvent/EvtRecDTag.h"
//#include "VertexFit/Helix.h"
#include "SimplePIDSvc/SimplePIDSvc.h"
#include "ParticleID/ParticleID.h"

#include "DstEvent/TofHitStatus.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include <fstream>
#include "TofCorrection.C"

		SimplePIDSvc::SimplePIDSvc(const std::string& name, ISvcLocator* svcLoc)
: Service(name, svcLoc)
{
		declareProperty("mindedxchi", m_dedxminchi=4);
		declareProperty("mintofchi", m_tofminchi=2);
		declareProperty("tofcorrection", m_tofcorrec=true);
}

SimplePIDSvc::~SimplePIDSvc()
{
}

StatusCode SimplePIDSvc::initialize()
{
		MsgStream log(messageService(), name());
		log << MSG::INFO << "@initialize()" << endreq;

		StatusCode sc = Service::initialize();

		sc = serviceLocator()->service("EventDataSvc", eventSvc_, true);
		//sc = serviceLocator()->service("VertexDbSvc", vtxSvc_, true);

		//m_dedxminchi=4;
		//m_tofminchi=2;

		//get eop files
		string filename;
		filename = string(getenv("SIMPLEPIDSVCROOT"))+"/share/eop.root";
		TFile* file=new TFile(filename.c_str(),"read");
		h_ebarlp=(TH1F*)file->Get("barlp");
		h_ebarhp=(TH1F*)file->Get("barhp");
		h_eendlp=(TH1F*)file->Get("endlp");
		h_eendhp=(TH1F*)file->Get("endhp");
		h_kbar=(TH1F*)file->Get("barkaon");
		h_kend=(TH1F*)file->Get("endkaon");
		h_pibar=(TH1F*)file->Get("barpion");
		h_piend=(TH1F*)file->Get("endpion");

		return sc;
}

StatusCode SimplePIDSvc::finalize()
{
		MsgStream log(messageService(), name());
		log << MSG::INFO << "@initialize()" << endreq;

		StatusCode sc = Service::finalize();

		delete h_ebarlp;
		delete h_ebarhp;
		delete h_eendlp;
		delete h_eendhp;
		delete h_kbar;
		delete h_kend;
		delete h_pibar;
		delete h_piend;


		return sc;
}

StatusCode SimplePIDSvc::queryInterface(const InterfaceID& riid, void** ppvIF)
{
		if ( ISimplePIDSvc::interfaceID().versionMatch(riid) ) {
				*ppvIF = dynamic_cast<ISimplePIDSvc*>(this);
		}
		else {
				return Service::queryInterface(riid, ppvIF);
		}
		addRef();
		return StatusCode::SUCCESS;
}

void SimplePIDSvc::preparePID(EvtRecTrack* track){

		//PID for proton
		ParticleID *pid = ParticleID::instance();
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(8);
		pid->setRecTrack(track);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTof()); // use PID sub-system
		//pid->identify(pid->onlyPion() | pid->onlyKaon());
		pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());
		pid->calculate();
		m_prob[0]=pid->probElectron();
		m_prob[1]=pid->probMuon();
		m_prob[2]=pid->probPion();
		m_prob[3]=pid->probKaon();
		m_prob[4]=pid->probProton();
		if(!(pid->IsPidInfoValid())) {for(int i=0; i<5; i++) m_prob[i]=-99; }

		//PID for kaon and pion
		ParticleID *PID = ParticleID::instance();
		PID->init();
		PID->setMethod(PID->methodProbability());
		PID->setChiMinCut(4);
		PID->setRecTrack(track);
		PID->usePidSys(PID->useDedx() | PID->useTof1() | PID->useTof2() | PID->useTof()); // use PID sub-system
		PID->identify(PID->onlyPion() | PID->onlyKaon());
		//PID->identify(PID->onlyPion() | PID->onlyKaon() | PID->onlyProton());
		PID->calculate();
		mm_prob[0]=PID->probElectron();
		mm_prob[1]=PID->probMuon();
		mm_prob[2]=PID->probPion();
		mm_prob[3]=PID->probKaon();
		mm_prob[4]=PID->probProton();
		if(!(PID->IsPidInfoValid())) {for(int i=0; i<5; i++) mm_prob[i]=-99; }


		//PID for electron 
		ParticleID *Epid = ParticleID::instance();
		Epid->init();
		Epid->setMethod(Epid->methodProbability());
		Epid->setChiMinCut(4);
		Epid->setRecTrack(track);
		Epid->usePidSys(Epid->useDedx() | Epid->useTof1() | Epid->useTof2() | Epid->useTof() | Epid->useEmc()); // use PID sub-system
		Epid->identify(Epid->onlyElectron() | Epid->onlyPion() | Epid->onlyKaon());
		Epid->calculate();
		mmm_prob[0]=Epid->probElectron();
		mmm_prob[1]=Epid->probMuon();
		mmm_prob[2]=Epid->probPion();
		mmm_prob[3]=Epid->probKaon();
		mmm_prob[4]=Epid->probProton();
		if(!(Epid->IsPidInfoValid())) {for(int i=0; i<5; i++) mmm_prob[i]=-99; }

}

bool SimplePIDSvc::iselectron(){

		if(mmm_prob[0]>0.001 && (mmm_prob[0]/(mmm_prob[0]+mmm_prob[2]+mmm_prob[3]))>0.8){ return true; }
		return false;
}
bool SimplePIDSvc::isproton(){
		if(m_prob[4]>=0 && m_prob[4]>=m_prob[3] && m_prob[4]>=m_prob[2]){ return true; }
		return false;
}

bool SimplePIDSvc::ispion(){

		if(mm_prob[2]>=0.00 && mm_prob[2]>=mm_prob[3]) return true;
//if(mm_prob[3]>=0.00 )  return true;
		return false;
}

bool SimplePIDSvc::iskaon(){

		if(mm_prob[3]>=0.00 && mm_prob[3]>=mm_prob[2]){ return true; }
//		if(mm_prob[3]>=0.00 ){ return true; }                      //loose kaon
		return false;
}

//bool SimplePIDSvc::ispion_for_Lambda_Ks(){
//
//		if(mm_prob[2]>=0.001 || mm_prob[2]>=mm_prob[3]) return true;
//		return false;
//}

bool SimplePIDSvc::ispion_for_Lambda_Ks(){

		if(mm_prob[2]>=0.00 && mm_prob[2]>=mm_prob[3]) return true;
		return false;
}

