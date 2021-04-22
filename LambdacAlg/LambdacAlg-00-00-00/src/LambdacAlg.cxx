//ucas-PRLi (lipeirong11@mails.ucas.ac.cn)
#include "LambdacAlg/LambdacAlg.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "BestDTagSvc/BestDTagSvc.h"
// #include "SimplePIDSvc/ISimplePIDSvc.h"
#include "LambdacAlg/MyPid.h"

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
using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
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
const double xmass[5] = {0.000511, 0.105658, 0.139571, 0.493677, 0.938272};
typedef std::vector<double> Vdouble;
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;

int Ntotal(0), Ncut0(0), Ncut1(0), Ncut2(0), Ncut3(0), Ncut4(0), Ncut5(0), Ncut6(0), Ncut7(0), H(0), A(0), all(0), al(0), Proton(0), all_m(0), all_p(0);

int Npp(0), Npm(0);

HepLorentzVector LambdacAlg::getP4(RecEmcShower *gTrk, Hep3Vector origin)
{
  Hep3Vector Gm_Vec(gTrk->x(), gTrk->y(), gTrk->z());
  Hep3Vector Gm_Mom = Gm_Vec - origin;
  Gm_Mom.setMag(gTrk->energy());
  HepLorentzVector pGm(Gm_Mom, gTrk->energy());
  return pGm;
}

LambdacAlg::LambdacAlg(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
  declareProperty("Vr0cut", m_vr0cut = 1.0);
  declareProperty("Vz0cut", m_vz0cut = 10.0);
  declareProperty("Vz0cut1", m_vz0cut1 = 20.0);
  declareProperty("SkimFlag", m_skim = false);
  declareProperty("PhotonMinEnergy", m_minEnergy = 0.025);
  declareProperty("GammaAngleCut", m_gammaAngleCut = 10.0);
  declareProperty("Test1C", m_test1C = 1); //ec fit
  declareProperty("GammathCut", m_gammathCut = 14.0);
  declareProperty("GammatlCut", m_gammatlCut = 0.0);
  declareProperty("PhotonMaxCosThetaBarrel", m_maxCosThetaBarrel = 0.8);
  declareProperty("PhotonMinCosThetaEndcap", m_minCosThetaEndcap = 0.86);
  declareProperty("PhotonMaxCosThetaEndcap", m_maxCosThetaEndcap = 0.92);
  declareProperty("PhotonMinEndcapEnergy", m_minEndcapEnergy = 0.050);
  declareProperty("Debug", m_debug = false);
  declareProperty("BeamE", m_beamE = 2.2997650);
  declareProperty("ReadBeamEFromDB", m_ReadBeamEFromDB = false);
  declareProperty("UseCalibBeamE", m_usecalibBeamE = false);
  declareProperty("CheckTotal", m_checktotal = 0);

  //
  declareProperty("ChisqMax", m_chisqMax = 200);
  declareProperty("Costheta", m_costheta = 0.93);
  declareProperty("EtaMinMass", m_EtaMinMass = 0.50);
  declareProperty("EtaMaxMass", m_EtaMaxMass = 0.56);
  declareProperty("Pi0MinMass", m_Pi0MinMass = 0.115);
  declareProperty("Pi0MaxMass", m_Pi0MaxMass = 0.15);
  declareProperty("SigmaMinMass", m_SigmaMinMass = 1.174);
  declareProperty("SigmaMaxMass", m_SigmaMaxMass = 1.2);
  declareProperty("EtaPrimeMinMass", m_EtaPrimeMinMass = 0.946);
  declareProperty("EtaPrimeMaxMass", m_EtaPrimeMaxMass = 0.968);

  declareProperty("OmegaMinMass", m_OmegaMinMass = 0.76);
  declareProperty("OmegaMaxMass", m_OmegaMaxMass = 0.8);
}
LambdacAlg::~LambdacAlg()
{
  //add your code for deconstructor
}
StatusCode LambdacAlg::initialize()
{
  if (m_debug)
    cout << __LINE__ << " begin initialize " << endl;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "LambdacAlg::initialize()" << endreq;
  m_beta.setX(0.011);
  m_beta.setY(0);
  m_beta.setZ(0);
  //add your code here
  StatusCode status;

  if (m_checktotal)
  {
    NTuplePtr nt2(ntupleSvc(), "FILE1/tree_truth");
    if (nt2)
      m_tuple2 = nt2;
    else
    {
      m_tuple2 = ntupleSvc()->book("FILE1/tree_truth", CLID_ColumnWiseTuple, "tree_truth N-Tuple example");
      if (m_tuple2)
      {
        status = m_tuple2->addItem("runNo", m_runNo_);
        status = m_tuple2->addItem("evtNo", m_evtNo_);
        status = m_tuple2->addItem("mode1", m_mode1_);
        status = m_tuple2->addItem("mode2", m_mode2_);
        status = m_tuple2->addItem("mode3", m_mode3_);
        status = m_tuple2->addItem("ndaughterAp", m_ndaughterAp_, 0, 15);
        status = m_tuple2->addIndexedItem("Ap_id", m_ndaughterAp_, m_Ap_id_);
        status = m_tuple2->addIndexedItem("Ap_ptruth", m_ndaughterAp_, 4, m_Ap_ptruth_);
        status = m_tuple2->addItem("ndaughterAm", m_ndaughterAm_, 0, 15);
        status = m_tuple2->addIndexedItem("Am_id", m_ndaughterAm_, m_Am_id_);
        status = m_tuple2->addIndexedItem("Am_ptruth", m_ndaughterAm_, 4, m_Am_ptruth_);
      }
      else
      {
        log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple2) << endmsg;
        return StatusCode::FAILURE;
      }
    }
  }

  NTuplePtr nt1(ntupleSvc(), "FILE1/tree");
  if (nt1)
    m_tuple1 = nt1;
  else
  {
    m_tuple1 = ntupleSvc()->book("FILE1/tree", CLID_ColumnWiseTuple, "exam N-Tuple example");
    if (m_tuple1)
    {
      status = m_tuple1->addItem("run", m_run);
      status = m_tuple1->addItem("event", m_event);
      status = m_tuple1->addItem("rightflag", m_rightflag);
      status = m_tuple1->addItem("signal", m_signal);
      status = m_tuple1->addItem("bg", m_bg);
      status = m_tuple1->addItem("yes", m_yes);
      status = m_tuple1->addItem("no", m_no);
      status = m_tuple1->addItem("indexmc", m_idxmc, 0, 100);
      status = m_tuple1->addIndexedItem("pdgid", m_idxmc, m_pdgid);
      status = m_tuple1->addIndexedItem("motheridx", m_idxmc, m_motheridx);
      status = m_tuple1->addItem("p4index", m_p4index, 0, 10);
      status = m_tuple1->addIndexedItem("Km", m_p4index, m_Km_p4);
      status = m_tuple1->addIndexedItem("Kp", m_p4index, m_Kp_p4);
      status = m_tuple1->addIndexedItem("p", m_p4index, m_p_p4);
      status = m_tuple1->addIndexedItem("pbar", m_p4index, m_pbar_p4);
      status = m_tuple1->addIndexedItem("pim", m_p4index, m_pim_p4);
      status = m_tuple1->addIndexedItem("pip", m_p4index, m_pip_p4);
      status = m_tuple1->addIndexedItem("gam1a", m_p4index, m_gam1a_p4);
      status = m_tuple1->addIndexedItem("gam2a", m_p4index, m_gam2a_p4);
      status = m_tuple1->addIndexedItem("gam3a", m_p4index, m_gam3a_p4);
      status = m_tuple1->addIndexedItem("gam4a", m_p4index, m_gam4a_p4);
      status = m_tuple1->addIndexedItem("pi0", m_p4index, m_pi0_p4);
      status = m_tuple1->addIndexedItem("gam1b", m_p4index, m_gam1b_p4);
      status = m_tuple1->addIndexedItem("gam2b", m_p4index, m_gam2b_p4);
      status = m_tuple1->addIndexedItem("gam3b", m_p4index, m_gam3b_p4);
      status = m_tuple1->addIndexedItem("gam4b", m_p4index, m_gam4b_p4);
      status = m_tuple1->addIndexedItem("gama", m_p4index, m_gama_p4);
      status = m_tuple1->addIndexedItem("gamb", m_p4index, m_gamb_p4);

      status = m_tuple1->addItem("pi0m", m_pi0m);
      status = m_tuple1->addItem("etam", m_etam);
      status = m_tuple1->addItem("ksi", m_ksi);
      status = m_tuple1->addItem("lambda", m_lambda);
      status = m_tuple1->addItem("etaprimem", m_etaprimem);
      status = m_tuple1->addItem("kstar", m_kstar);
      status = m_tuple1->addItem("sigmastar", m_sigmastar);
      status = m_tuple1->addItem("all", m_all);
      status = m_tuple1->addItem("Kmindex", m_Kmindex);
      status = m_tuple1->addItem("pbarindex", m_pbarindex);
      status = m_tuple1->addItem("pindex", m_pindex);
      status = m_tuple1->addItem("Kpindex", m_Kpindex);
      status = m_tuple1->addItem("r", m_r);
      status = m_tuple1->addItem("pi0m1c", m_pi0m1c);
      status = m_tuple1->addItem("etam1c", m_etam1c);
      status = m_tuple1->addItem("sigmam", m_Sigmam);
      status = m_tuple1->addItem("num_othertrackp", m_numothertrackp);
      status = m_tuple1->addItem("num_othertrackm", m_numothertrackm);
      status = m_tuple1->addItem("E_beam", m_ebeam);
      status = m_tuple1->addItem("deltaE_min", m_deltaE_min);
      status = m_tuple1->addItem("M_BC", m_bc);
      status = m_tuple1->addItem("np", m_np);
      status = m_tuple1->addItem("npbar", m_npbar);
      status = m_tuple1->addItem("eop_pim", m_eop_pim);
      status = m_tuple1->addItem("eop_pip", m_eop_pip);
      status = m_tuple1->addItem("eop_p", m_eop_p);
      status = m_tuple1->addItem("eop_pbar", m_eop_pbar);
      status = m_tuple1->addItem("eop_Km", m_eop_Km);
      status = m_tuple1->addItem("eop_Kp", m_eop_Kp);
      status = m_tuple1->addItem("ngam", m_ngam);
      status = m_tuple1->addItem("chi2", m_chi2);

      status = m_tuple1->addItem("indexmc_p", m_mcparticle_p, 0, 100);
      status = m_tuple1->addIndexedItem("pdgid_p", m_mcparticle_p, m_pdgid_p);
      status = m_tuple1->addIndexedItem("motheridx_p", m_mcparticle_p, m_motheridx_p);
      status = m_tuple1->addItem("indexmc_m", m_mcparticle_m, 0, 100);
      status = m_tuple1->addIndexedItem("pdgid_m", m_mcparticle_m, m_pdgid_m);
      status = m_tuple1->addIndexedItem("motheridx_m", m_mcparticle_m, m_motheridx_m);
      status = m_tuple1->addItem("ndaughterAp", m_ndaughterAp, 0, 15);
      status = m_tuple1->addIndexedItem("Ap_id", m_ndaughterAp, m_Ap_id);
      status = m_tuple1->addIndexedItem("Ap_ptruth", m_ndaughterAp, 4, m_Ap_ptruth);
      status = m_tuple1->addItem("ndaughterAm", m_ndaughterAm, 0, 15);
      status = m_tuple1->addIndexedItem("Am_id", m_ndaughterAm, m_Am_id);
      status = m_tuple1->addIndexedItem("Am_ptruth", m_ndaughterAm, 4, m_Am_ptruth);
      status = m_tuple1->addItem("mode1", m_mode1);
      status = m_tuple1->addItem("mode2", m_mode2);
      status = m_tuple1->addItem("mode3", m_mode3);
    }
    else
    {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
      return StatusCode::FAILURE;
    }
  }

  //--------end of book--------
  //

  log << MSG::INFO << "successfully return from initialize()" << endmsg;
  return StatusCode::SUCCESS;
}
StatusCode LambdacAlg::beginRun()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "LambdacAlg::beginRun()" << endreq;
  //add your code here
  return StatusCode::SUCCESS;
}
StatusCode LambdacAlg::execute()
{
  m_rightflag = -1;
  int signal = -999;
  int bg = -1;
  int yes = -1;
  int no = -1;
  Ntotal++;
  if (m_debug)
    cout << "\033[31m" << "m_debug3 begin execute, Ntotal(event): " << Ntotal << "\033[0m" << endl;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "LambdacAlg::execute()" << endreq;
  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");

  m_run = eventHeader->runNumber();
  m_event = eventHeader->eventNumber();
  int runNo = eventHeader->runNumber();
  int eventNo = eventHeader->eventNumber();
  int mm_mode1 = ((eventHeader->flag1() / 1000000)) % 1000;
  int mm_mode2 = (eventHeader->flag1() / 1000) % 1000;
  int mm_mode3 = eventHeader->flag1() % 1000;

  //if(!(m_run==35966&&(m_event==10304||m_event==77483||m_event==79947||m_event==97315)))   return StatusCode::SUCCESS;
  //cout<<endl<<endl<<"********************************************************************************************************"<<endl;

  log << MSG::DEBUG << "run, evtnum = " << m_run << " , " << m_event << endreq;

  IMcDecayModeSvc *i_svc;
  StatusCode sc_DecayModeSvc = service("McDecayModeSvc", i_svc);
  if (sc_DecayModeSvc.isFailure())
  {
    log << MSG::FATAL << "Could not load McDecayModeSvc!" << endreq;
    return sc_DecayModeSvc;
  }
  m_svc = dynamic_cast<McDecayModeSvc *>(i_svc);

  int gamp = 0;
  int gampmi = 0;
  int game = 0;
  int gamemi = 0;
  int pim = 0;
  int piplus = 0;
  int pimmi = 0;
  int piplusmi = 0;
  int pi0 = 0;
  int pi0mi = 0;
  int proton = 0;
  int pm = 0;
  int sigm = 0;
  int sigmm = 0;
  int eta = 0;

  int numParticle_p = 0;
  int numParticle_m = 0;
  int numParticle = 0;
  int M_pdgid_p[100];
  int M_motheridx_p[100];

  int M_pdgid_m[100];
  int M_motheridx_m[100];

  int M_pdgid[100];
  int M_motheridx[100];

  int ndaughterAp = 0;
  double Ap_ptruth[15][4];
  int Ap_id[15];
  for (int aa = 0; aa < 15; aa++)
    for (int ll = 0; ll < 4; ll++)
      Ap_ptruth[aa][ll] = 0;
  for (int aa = 0; aa < 15; aa++)
    Ap_id[aa] = 0;

  int ndaughterAm = 0;
  double Am_ptruth[15][4];
  int Am_id[15];
  for (int aa = 0; aa < 15; aa++)
    for (int ll = 0; ll < 4; ll++)
      Am_ptruth[aa][ll] = 0;
  for (int aa = 0; aa < 15; aa++)
    Am_id[aa] = 0;

  Vint me;
  me.clear();
  if (eventHeader->runNumber() < 0)
  {
    SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
    Vint pdgid;
    ;
    Vint motherindex;
    pdgid.clear();
    motherindex.clear();

    int numParticle = 0;
    if (!mcParticleCol)
    {
      std::cout << "Could not retrieve McParticelCol" << std::endl;
      // return StatusCode::FAILURE;
    }
    /*else
    {
      Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
      for (; iter_mc != mcParticleCol->end(); iter_mc++)
      {

        if ((*iter_mc)->primaryParticle())
          continue;
        if (!(*iter_mc)->decayFromGenerator())
          continue;
        int pdg = (*iter_mc)->particleProperty();
        int motherpdg = ((*iter_mc)->mother()).particleProperty();
        int mmotherpdg = (((*iter_mc)->mother()).mother()).particleProperty();
        if ((*iter_mc)->particleProperty() == 4122)
        {
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();
            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          } // End of "gc.size() > 0" IF
        }
        if ((*iter_mc)->particleProperty() == -4122)
        {
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();
            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }

        //mode2
        if (mm_mode2 == 67 && pdg == 3222 && motherpdg == 4122)
        { //Sgm+ ->p pi0
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();
            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          } // End of "gc.size() > 0" IF
        }

        if (mm_mode2 == 67 && pdg == 331 && motherpdg == 4122)
        { //eta' ->pip pim eta
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();
            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          } // End of "gc.size() > 0" IF
        }

        if ((mm_mode2 == 67) && pdg == 111 && motherpdg == 3222 && mmotherpdg == 4122)
        { ////Sgm+ ->p pi0, pi0->gam gam //store 2 gamma
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();
            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          } // End of "gc.size() > 0" IF
        }

        if (mm_mode2 == 67 && pdg == 221 && motherpdg == 331 && mmotherpdg == 4122)
        { //eta' ->pip pim eta, eta->gam gam // 2 gamma
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();
            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          } // End of "gc.size() > 0" IF
        }

        //mode3
        if (mm_mode3 == 67 && pdg == -3222 && motherpdg == -4122)
        { //Sgm- ->pbar pi0
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();
            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }
        if (mm_mode3 == 67 && pdg == 331 && motherpdg == -4122)
        { //eta' ->pip pim eta
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();
            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }

        if ((mm_mode3 == 67) && pdg == 111 && motherpdg == -3222 && mmotherpdg == -4122)
        { ////Sgm- ->pbar pi0, pi0->gam gam //store 2 gamma
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();
            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }
        if (mm_mode3 == 67 && pdg == 221 && motherpdg == 331 && mmotherpdg == -4122)
        { //eta' ->pip pim eta, eat->gam gam //two gamma
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();
            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }
      }

      //trace Lambda_C+

      for (iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++)
      {
        if ((*iter_mc)->primaryParticle())
          continue;
        if (!(*iter_mc)->decayFromGenerator())
          continue;

        if ((*iter_mc)->particleProperty() == 4122)
        {
          int mode = m_svc->extract(*iter_mc, pdgid, motherindex);
          numParticle = pdgid.size();
          for (int i = 0; i != pdgid.size(); i++)
          {
            M_pdgid[i] = pdgid[i];
            if (m_debug)
              cout << "M_pdgid[i]=" << M_pdgid[i] << endl;
            M_motheridx[i] = motherindex[i];
            if (m_debug)
              cout << "M_motheridx[i]=" << M_motheridx[i] << endl;
          }
        }
        pdgid.clear();
        motherindex.clear();

        //trace Lambda_C-
        if ((*iter_mc)->particleProperty() == -4122)
        {
          int mode = m_svc->extract(*iter_mc, pdgid, motherindex);
          numParticle_m = pdgid.size();
          for (int i = 0; i != pdgid.size(); i++)
          {
            M_pdgid_m[i] = pdgid[i];
            if (m_debug)
              cout << "M_pdgid_m[i]=" << M_pdgid_m[i] << endl;
            M_motheridx_m[i] = motherindex[i];
            if (m_debug)
              cout << "M_motheridx_m[i]=" << M_motheridx_m[i] << endl;
          }
        }
      }

      //if(mm_mode2==67&&ndaughterAp==11&&Ap_id[0]==3222&&Ap_id[1]==331&&Ap_id[2]==2212&&Ap_id[3]==111&&Ap_id[4]==211&&Ap_id[5]==-211&&Ap_id[6]==221&&Ap_id[7]==22&&Ap_id[8]==22&&Ap_id[9]==22&&Ap_id[10]==22)
      //{
      //	signal=1;
      //	all++;
      //}
      //
      //					 if (mm_mode3==67&&ndaughterAm==11&&Am_id[0]==-3222&&Am_id[1]==331&&Am_id[2]==-2212&&Am_id[3]==111&&Am_id[4]==211&&Am_id[5]==-211&&Am_id[6]==221&&Am_id[7]==22&&Am_id[7]==22&&Am_id[9]==22&&Am_id[10]==22)
      //						{
      //							signal=1;
      //							all++;
      //							}

      if (mm_mode3 == 0)
      {
        signal = 1;
        all++;
      }
      if (mm_mode2 == 0)
      {
        signal = 1;
        H++;
      }
    }*/
    else
    {
      Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
      //
      for (; iter_mc != mcParticleCol->end(); iter_mc++)
      {
        if ((*iter_mc)->primaryParticle())
          continue;
        if (!(*iter_mc)->decayFromGenerator())
          continue;

        int pdg = (*iter_mc)->particleProperty();
        int motherpdg = ((*iter_mc)->mother()).particleProperty();
        int mmotherpdg = (((*iter_mc)->mother()).mother()).particleProperty();
        if (m_debug)
          cout << __LINE__ << "mcParticleCol pdg: " << pdg << ", motherpdg: " << motherpdg
               << ", mmotherpdg:" << mmotherpdg << endl;

        // truth for lambda_c+ to ... ------------------------------------
        if ((*iter_mc)->particleProperty() == 4122)
        { // lambdac+ -> sigma+ omega, store sigma+ omega
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ << " ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;

          } // End of "gc.size() > 0" IF
        }

        if (pdg == 3222 && motherpdg == 4122)
        { // Sgm+ ->p pi0, store p pi0
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ << " ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          } // End of "gc.size() > 0" IF
        }
        if (pdg == 111 && motherpdg == 3222 && mmotherpdg == 4122)
        { ////Sgm+ ->p pi0, pi0->gam gam //store 2 gamma
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ << " ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          } // End of "gc.size() > 0" IF
        }
        if (pdg == 223 && motherpdg == 4122)
        { //// omega ->pi+ pi- eta, store pi+ pi- pi0
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ << " ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          }
        }

        if (pdg == 111 && motherpdg == 223 && mmotherpdg == 4122)
        { // pi0 ->gam gam //store 2 gamma
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ << " ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          } // End of "gc.size() > 0" IF
        }

        // truth for lambda_c- to ... ------------------------------------
        if ((*iter_mc)->particleProperty() == -4122)
        { // lambdac- -> sigma- oemga, store sigma- omega
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ << " ndaughterAm[" << ndaughterAm << "]: " << Am_id[ndaughterAm] << endl;

            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }
        // mode3
        if (pdg == -3222 && motherpdg == -4122)
        { // Sgm- ->pbar pi0, store pbar pi0
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ << " ndaughterAm[" << ndaughterAm << "]: " << Am_id[ndaughterAm] << endl;

            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }
        if (pdg == 111 && motherpdg == -3222 && mmotherpdg == -4122)
        { ////Sgm- ->pbar pi0, pi0->gam gam //store 2 gamma
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ << " ndaughterAm[" << ndaughterAm << "]: " << Am_id[ndaughterAm] << endl;

            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }
        if (pdg == 223 && motherpdg == -4122)
        { //// omega ->pi+ pi- pi0, store pi+ pi- pi0
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ << " ndaughterAm[" << ndaughterAm << "]: " << Am_id[ndaughterAm] << endl;

            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          }
        }

        if (pdg == 111 && motherpdg == 223 && mmotherpdg == -4122)
        { // pi0 ->gam gam //store 2 gamma
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ << " ndaughterAm[" << ndaughterAm << "]: " << Am_id[ndaughterAm] << endl;

            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }
        if (m_debug)
          cout << endl;
      }

      if (m_debug)
      {
        cout << endl;
        for (int i = 0; i < ndaughterAp; i++)
        {
          cout << "Ap_id[" << i << "]: " << Ap_id[i] << ", ";
        }
        cout << endl;

        for (int i = 0; i < ndaughterAm; i++)
        {
          cout << "Am_id[" << i << "]: " << Am_id[i] << ", ";
        }
        cout << "\n" << endl;
      }

      // trace both
      for (iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++)
      {
        if ((*iter_mc)->primaryParticle())
          continue;
        if (!(*iter_mc)->decayFromGenerator())
          continue;

        if ((*iter_mc)->particleProperty() == 4122 || (*iter_mc)->particleProperty() == -4122)
        {
          int mode = m_svc->extract(*iter_mc, pdgid, motherindex);
          numParticle = pdgid.size();

          if (m_debug)
          {
            cout << endl;
            cout << __LINE__ << " +||- iter_mc: " << (*iter_mc)->particleProperty() << ", numParticle: " << numParticle
                 << endl;
          }

          for (int i = 0; i != pdgid.size(); i++)
          {
            M_pdgid[i] = pdgid[i];
            M_motheridx[i] = motherindex[i];

            if (m_debug)
            {
              cout << __LINE__ << "M_pdgid[" << i << "]=" << M_pdgid[i] << endl;
              cout << __LINE__ << "M_motheridx[" << i << "]=" << M_motheridx[i] << endl;
            }
          }
        }

        pdgid.clear();
        motherindex.clear();

        // trace Lambda_C+
        if ((*iter_mc)->particleProperty() == 4122)
        {
          int mode = m_svc->extract(*iter_mc, pdgid, motherindex);
          numParticle_p = pdgid.size();

          if (m_debug)
          {
            cout << endl;
            cout << __LINE__ << __LINE__ << "+iter_mc: " << (*iter_mc)->particleProperty()
                 << ", numParticle: " << numParticle << endl;
          }

          for (int i = 0; i != pdgid.size(); i++)
          {
            M_pdgid_p[i] = pdgid[i];
            M_motheridx_p[i] = motherindex[i];

            if (m_debug)
            {
              cout << __LINE__ << "M_pdgid[" << i << "]=" << M_pdgid[i] << endl;
              cout << __LINE__ << "M_motheridx[" << i << "]=" << M_motheridx[i] << endl;
            }
          }
        }

        pdgid.clear();
        motherindex.clear();
        // trace Lambda_C-
        if ((*iter_mc)->particleProperty() == -4122)
        {
          int mode = m_svc->extract(*iter_mc, pdgid, motherindex);
          numParticle_m = pdgid.size();

          if (m_debug)
          {
            cout << endl;
            cout << __LINE__ << "-iter_mc: " << (*iter_mc)->particleProperty() << ", numParticle: " << numParticle
                 << endl;
          }

          for (int i = 0; i != pdgid.size(); i++)
          {
            M_pdgid_m[i] = pdgid[i];
            M_motheridx_m[i] = motherindex[i];

            if (m_debug)
            {
              cout << __LINE__ << "M_pdgid[" << i << "]=" << M_pdgid[i] << endl;
              cout << __LINE__ << "M_motheridx[" << i << "]=" << M_motheridx[i] << endl;
            }
          }
        }
      }

      if (ndaughterAp == 11 && Ap_id[0] == 3222 && Ap_id[1] == 223 && Ap_id[2] == 2212 && Ap_id[3] == 111 &&
          Ap_id[4] == 211 && Ap_id[5] == -211 && Ap_id[6] == 111 && Ap_id[7] == 22 && Ap_id[8] == 22 &&
          Ap_id[9] == 22 && Ap_id[10] == 22)
      {
        signal = 1;
        all++;
        all_p++;
      }
      if (ndaughterAm == 11 && Am_id[0] == -3222 && Am_id[1] == 223 && Am_id[2] == -2212 && Am_id[3] == 111 &&
          Am_id[4] == 211 && Am_id[5] == -211 && Am_id[6] == 111 && Am_id[7] == 22 && Am_id[8] == 22 &&
          Am_id[9] == 22 && Am_id[10] == 22)
      {
        signal = -1;
        all++;
        all_m++;
      }

      if (m_debug)
      {
        cout << __LINE__ << __LINE__ << " all: " << all << ", +signal: " << all_p << ", -signal: " << all_m << endl;
        cout << endl;
      }
    }
  }

  

  if (m_checktotal)
  {
    m_runNo_ = runNo;
    m_evtNo_ = eventNo;
    m_mode1_ = mm_mode1;
    m_mode2_ = mm_mode2;
    m_mode3_ = mm_mode3;
    m_ndaughterAp_ = ndaughterAp;
    for (int aa = 0; aa < ndaughterAp; aa++)
      m_Ap_id_[aa] = Ap_id[aa];
    for (int aa = 0; aa < ndaughterAp; aa++)
      for (int ll = 0; ll < 4; ll++)
      {
        m_Ap_ptruth_[aa][ll] = Ap_ptruth[aa][ll];

        if (m_debug)
          cout << __LINE__ << "    " << Ap_id[aa] << ", " << Ap_ptruth[aa][0] << " " << Ap_ptruth[aa][1] << " " << Ap_ptruth[aa][2] << " " << Ap_ptruth[aa][3] << " " << endl;
      }

    m_ndaughterAm_ = ndaughterAm;
    for (int aa = 0; aa < ndaughterAm; aa++)
      m_Am_id_[aa] = Am_id[aa];
    for (int aa = 0; aa < ndaughterAm; aa++)
      for (int ll = 0; ll < 4; ll++)
      {
        m_Am_ptruth_[aa][ll] = Am_ptruth[aa][ll];
        if (m_debug)
          cout << __LINE__ << "    " << Am_id[aa] << ", " << Am_ptruth[aa][0] << " " << Am_ptruth[aa][1] << " " << Am_ptruth[aa][2] << " " << Am_ptruth[aa][3] << " " << endl;
      }
    m_tuple2->write();
  }
  //end mc
  //TADA!!!=========================================================================
  //good track
  if (m_debug)
    cout << __LINE__ << " begin choose good track" << endl;

  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);

  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);

  Vint iGoodtotal, iGoodforp, othertrack;
  iGoodforp.clear(), iGoodtotal.clear(), othertrack.clear();

  int nCharge = 0;
  Hep3Vector xorigin(0, 0, 0);
  IVertexDbSvc *vtxsvc;
  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if (vtxsvc->isVertexValid())
  {
    double *dbv = vtxsvc->PrimaryVertex();
    double *vv = vtxsvc->SigmaPrimaryVertex();
    xorigin.setX(dbv[0]);
    xorigin.setY(dbv[1]);
    xorigin.setZ(dbv[2]);
  }
  if (m_debug)
    cout << __LINE__ << "choose good track" << endl;
  for (int i = 0; i < evtRecEvent->totalCharged(); i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
    if (!(*itTrk)->isMdcTrackValid())
      continue;
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();

    double pch = mdcTrk->p();
    //double x0=mdcTrk->x();
    //double y0=mdcTrk->y();
    //double z0=mdcTrk->z();
    //double phi0=mdcTrk->helix(1);
    //double xv=xorigin.x();
    //double yv=xorigin.y();
    //double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);
    HepVector a = mdcTrk->helix();
    HepSymMatrix Ea = mdcTrk->err();
    HepPoint3D point0(0., 0., 0.); // the initial point for MDC recosntruction
    HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
    VFHelix helixip(point0, a, Ea);
    helixip.pivot(IP);
    HepVector vecipa = helixip.a();
    double Rvxy0 = fabs(vecipa[0]); //the nearest distance to IP in xy plane
    double Rvz0 = vecipa[3];        //the nearest distance to IP in z direction
    double Rvphi0 = vecipa[1];
    //A++;
    double costheta = cos(mdcTrk->theta());
    if (fabs(costheta) >= m_costheta)
      continue;
    if (fabs(Rvz0) >= m_vz0cut)
      continue;
    if (fabs(Rvxy0) >= m_vr0cut)
      continue;

    iGoodforp.push_back(i);
    iGoodtotal.push_back(i);
    nCharge += mdcTrk->charge();
  } // Finish Good Charged Track SKction
  if (m_debug)
    cout << __LINE__ << endl;
  int nGoodforp = iGoodforp.size();
  if (nGoodforp < 1)
    return StatusCode::SUCCESS;
  //  int nothertrack = othertrack.size();
  Vint ip, ipbar;
  ip.clear();
  ipbar.clear();
  for (int i = 0; i < nGoodforp; i++)
  {
    cout <<__LINE__<< " " <<  iGoodforp[i] << endl;
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGoodforp[i];
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
    if (mdcTrk->charge() == 1)
    {
      // ISimplePIDSvc *m_simplePIDSvc;
      // Gaudi::svcLocator()->service("SimplePIDSvc", m_simplePIDSvc);
      // m_simplePIDSvc->preparePID(*itTrk);
      MyPid *m_simplePIDSvc = new MyPid(*itTrk);
      if (m_simplePIDSvc->isproton())
      {        
        ip.push_back(iGoodforp[i]);
        cout <<__LINE__<< " " <<  iGoodforp[i] << endl;
      }    
    }
    if (mdcTrk->charge() == -1)
    {
      // ISimplePIDSvc *m_simplePIDSvc;
      // Gaudi::svcLocator()->service("SimplePIDSvc", m_simplePIDSvc);
      // m_simplePIDSvc->preparePID(*itTrk);
      MyPid *m_simplePIDSvc = new MyPid(*itTrk);
      if (m_simplePIDSvc->isproton())
      {        
        ipbar.push_back(iGoodforp[i]);
        cout <<__LINE__<< " " <<  iGoodforp[i] << endl;
      }    
    }
  }

  int nGoodtotal = iGoodtotal.size();
  if (nGoodtotal < 3)
    return StatusCode::SUCCESS;
  if (fabs(signal) == 1)
    Ncut0++;
  if (m_debug)
    cout << __LINE__ << endl;

  Vint iKp, ipip, iKm, ipim;
  iKp.clear();
  ipip.clear();
  iKm.clear();
  ipim.clear();
  int ntrackm(0), ntrackp(0);
  for (int i = 0; i < nGoodtotal; i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGoodtotal[i];
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
    if (mdcTrk->charge() == 1)
    {
      ntrackm++;
      // ISimplePIDSvc *m_simplePIDSvc;
      // Gaudi::svcLocator()->service("SimplePIDSvc", m_simplePIDSvc);
      // m_simplePIDSvc->preparePID(*itTrk);
      MyPid *m_simplePIDSvc = new MyPid(*itTrk);
      if (m_simplePIDSvc->ispion())
        ipip.push_back(iGoodtotal[i]);
      if (m_simplePIDSvc->iskaon())
        iKp.push_back(iGoodtotal[i]);
    }

    if (mdcTrk->charge() == -1)
    {
      ntrackp++;
      // ISimplePIDSvc *m_simplePIDSvc;
      // Gaudi::svcLocator()->service("SimplePIDSvc", m_simplePIDSvc);
      // m_simplePIDSvc->preparePID(*itTrk);
      MyPid *m_simplePIDSvc = new MyPid(*itTrk);
      if (m_simplePIDSvc->iskaon())
        iKm.push_back(iGoodtotal[i]);
      if (m_simplePIDSvc->ispion())
        ipim.push_back(iGoodtotal[i]);
    }
  }
  if (m_debug)
    cout << __LINE__ << endl;
  int npim = ipim.size();
  int nKp = iKp.size();
  int np = ip.size();
  int npip = ipip.size();
  int nKm = iKm.size();
  int npbar = ipbar.size();
  m_np = np;
  m_npbar = npbar;
  //	cout<<"npi="<<npi<<"np="<<np<<"nK="<<nK<<nGood<<endl;
  //	if (np>0&&nKp>0&&npim>0)     Proton++;
  if ((npim == 0 || npip == 0 || np == 0) && (npim == 0 || npip == 0 || npbar == 0))
    return StatusCode::SUCCESS;
  //if (np==0&&npbar==0) return StatusCode::SUCCESS;

  if (m_debug)
    cout << __LINE__ << endl;
  if (fabs(signal) == 1)
    Ncut1++;

  Vp4 ppim, ppip, pp, ppbar, pKp, pKm;
  ppim.clear();
  ppip.clear();
  ppbar.clear();
  pp.clear();
  pKp.clear();
  pKm.clear();
  for (int i = 0; i < npim; i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ipim[i];
    RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
    mdcKalTrk->setPidType(RecMdcKalTrack::pion); //still....
    HepLorentzVector pim4 = mdcKalTrk->p4(xmass[2]);
    double pim3 = pim4.rho();
    ppim.push_back(pim4);
    double pimenergy = 0.;
    double eop_pim = 0.;
    if ((*itTrk)->isEmcShowerValid())
    {
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      pimenergy = emcTrk->energy();
      eop_pim = pimenergy / (pim3);
    }

    m_eop_pim = eop_pim;
  }
  if (m_debug)
    cout << __LINE__ << endl;
  for (int i = 0; i < npip; i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ipip[i];
    RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
    mdcKalTrk->setPidType(RecMdcKalTrack::pion); //still....
    HepLorentzVector pip4 = mdcKalTrk->p4(xmass[2]);
    double pip3 = pip4.rho();
    ppip.push_back(pip4);
    double pipenergy = 0.;
    double eop_pip = 0.;
    if ((*itTrk)->isEmcShowerValid())
    {
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      pipenergy = emcTrk->energy();
      eop_pip = pipenergy / (pip3);
    }

    m_eop_pip = eop_pip;
  }
  if (m_debug)
    cout << __LINE__ << endl;
  for (int i = 0; i < np; i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ip[i];
    RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
    mdcKalTrk->setPidType(RecMdcKalTrack::proton);
    HepLorentzVector pp4 = mdcKalTrk->p4(xmass[4]);
    double p3 = pp4.rho();
    pp.push_back(pp4);
    double penergy = 0.;
    double eop_p = 0.;
    if ((*itTrk)->isEmcShowerValid())
    {
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      penergy = emcTrk->energy();
      eop_p = penergy / (p3);
    }

    m_eop_p = eop_p;
  }
  if (m_debug)
    cout << __LINE__ << endl;
  for (int i = 0; i < npbar; i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ipbar[i];
    RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
    mdcKalTrk->setPidType(RecMdcKalTrack::proton);
    HepLorentzVector ppbar4 = mdcKalTrk->p4(xmass[4]);
    double pbar3 = ppbar4.rho();
    ppbar.push_back(ppbar4);
    double pbarenergy = 0.;
    double eop_pbar = 0.;
    if ((*itTrk)->isEmcShowerValid())
    {
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      pbarenergy = emcTrk->energy();
      eop_pbar = pbarenergy / (pbar3);
    }

    m_eop_pbar = eop_pbar;
  }
  if (m_debug)
    cout << __LINE__ << endl;
  for (int i = 0; i < nKm; i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iKm[i];
    RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
    mdcKalTrk->setPidType(RecMdcKalTrack::kaon); //oops!!!!
    HepLorentzVector pKm4 = mdcKalTrk->p4(xmass[3]);
    double Km3 = pKm4.rho();
    pKm.push_back(pKm4);
    double Kmenergy = 0.;
    double eop_Km = 0.;
    if ((*itTrk)->isEmcShowerValid())
    {
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      Kmenergy = emcTrk->energy();
      eop_Km = Kmenergy / (Km3);
    }

    m_eop_Km = eop_Km;
  }
  if (m_debug)
    cout << __LINE__ << endl;
  for (int i = 0; i < nKp; i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iKp[i];
    RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
    mdcKalTrk->setPidType(RecMdcKalTrack::kaon); //oops!!!!
    HepLorentzVector pKp4 = mdcKalTrk->p4(xmass[3]);
    double Kp3 = pKp4.rho();
    pKp.push_back(pKp4);
    double Kpenergy = 0.;
    double eop_Kp = 0.;
    if ((*itTrk)->isEmcShowerValid())
    {
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      Kpenergy = emcTrk->energy();
      eop_Kp = Kpenergy / (Kp3);
    }

    m_eop_Kp = eop_Kp;
  }
  if (m_debug)
    cout << __LINE__ << endl;

  //good shower
  Vint iGam;
  iGam.clear();
  if (m_debug)
    cout << __LINE__ << "begin choose good gamma" << endl;
  for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
    if (!(*itTrk)->isEmcShowerValid())
      continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    Hep3Vector emcp(emcTrk->x(), emcTrk->y(), emcTrk->z());
    HepLorentzVector shP4 = getP4(emcTrk, xorigin);
    double cosThetaSh = shP4.vect().cosTheta();
    //double dthe = 200.;
    //double dphi = 200.;
    double dang = 200.;
    if (m_debug)
      cout << __LINE__ << "choose good gamma" << endl;
    for (int j = 0; j < evtRecEvent->totalCharged(); j++)
    {
      EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
      if (!(*jtTrk)->isExtTrackValid())
        continue;
      RecExtTrack *extTrk = (*jtTrk)->extTrack();
      if (extTrk->emcVolumeNumber() == -1)
        continue;
      Hep3Vector extp = extTrk->emcPosition();
      double angd = extp.angle(emcp);
      //double thed = extp.theta() - emcp.theta();
      //double phid = extp.deltaPhi(emcp);
      //thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      //phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      if (angd < dang)
      {
        dang = angd;
        //dthe = thed;
        //dphi = phid;
      }
    }
    if (dang >= 200)
      continue;
    double eraw = emcTrk->energy();
    double getTime = emcTrk->time();
    //dthe = dthe * 180 / (CLHEP::pi);
    //dphi = dphi * 180 / (CLHEP::pi);
    dang = dang * 180 / (CLHEP::pi);
    if (fabs(dang) < m_gammaAngleCut)
      continue;
    if (getTime > m_gammathCut || getTime < m_gammatlCut)
      continue;
    if (!((fabs(cosThetaSh) < m_maxCosThetaBarrel && eraw > m_minEnergy) || ((fabs(cosThetaSh) > m_minCosThetaEndcap) && (fabs(cosThetaSh) < m_maxCosThetaEndcap) && (eraw > m_minEndcapEnergy))))
      continue;
    iGam.push_back(i);
  } // Finish Good Photon SKction
  if (m_debug)
    cout << __LINE__ << endl;
  int nGam = iGam.size();
  if (nGam < 4)
    return StatusCode::SUCCESS;
  if (fabs(signal) == 1)
    Ncut2++;

  if (m_debug)
    cout << __LINE__ << endl;
  Vint igam1, igam2, igam3, igam4;
  igam1.clear(), igam2.clear(), igam3.clear(), igam4.clear();
  Vp4 pgam1, pgam2, pgam3, pgam4, pgam1_1C, pgam2_1C, pgam3_1C, pgam4_1C;
  pgam1.clear(), pgam2.clear(), pgam3.clear(), pgam4.clear(), pgam1_1C.clear(), pgam2_1C.clear(), pgam3_1C.clear(), pgam4_1C.clear();
  Vdouble chi2;
  chi2.clear();
  KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
  // Loop each gamma pair, check ppi0 mass
  for (int i = 0; i < nGam - 1; i++)
  {
    EvtRecTrackIterator itTrki = evtRecTrkCol->begin() + iGam[i];
    RecEmcShower *emcTrki = (*itTrki)->emcShower();
    Hep3Vector emcposi(emcTrki->x(), emcTrki->y(), emcTrki->z());
    HepLorentzVector ptrki = getP4(emcTrki, xorigin);
    for (int j = i + 1; j < nGam; j++)
    {
      EvtRecTrackIterator itTrkj = evtRecTrkCol->begin() + iGam[j];
      RecEmcShower *emcTrkj = (*itTrkj)->emcShower();
      Hep3Vector emcposj(emcTrkj->x(), emcTrkj->y(), emcTrkj->z());
      HepLorentzVector ptrkj = getP4(emcTrkj, xorigin);

      HepLorentzVector p2geta = ptrki + ptrkj;
      if (p2geta.m() < 0.115 || p2geta.m() > 0.15)
        continue;
      if (m_test1C == 1)
      {
        kmfit->init();
        kmfit->setChisqCut(m_chisqMax);
        kmfit->AddTrack(0, 0.0, emcTrki);
        kmfit->AddTrack(1, 0.0, emcTrkj);
        kmfit->AddResonance(0, 0.1349770, 0, 1);
        bool oksq = kmfit->Fit();
        if (oksq)
        {
          igam1.push_back(iGam[i]);
          pgam1.push_back(ptrki);
          pgam1_1C.push_back(kmfit->pfit(0));
          igam2.push_back(iGam[j]);
          pgam2.push_back(ptrkj);
          pgam2_1C.push_back(kmfit->pfit(1));
        }
      }
    }
  }

  for (int k = 0; k < nGam - 1; k++)
  {
    EvtRecTrackIterator itTrkk = evtRecTrkCol->begin() + iGam[k];
    RecEmcShower *emcTrkk = (*itTrkk)->emcShower();
    Hep3Vector emcposk(emcTrkk->x(), emcTrkk->y(), emcTrkk->z());
    HepLorentzVector ptrkk = getP4(emcTrkk, xorigin);
    for (int l = k + 1; l < nGam; l++)
    {
      EvtRecTrackIterator itTrkl = evtRecTrkCol->begin() + iGam[l];
      RecEmcShower *emcTrkl = (*itTrkl)->emcShower();
      Hep3Vector emcposl(emcTrkl->x(), emcTrkl->y(), emcTrkl->z());
      HepLorentzVector ptrkl = getP4(emcTrkl, xorigin);

      HepLorentzVector p2gpi = ptrkk + ptrkl;
      if (p2gpi.m() < 0.115 || p2gpi.m() > 0.15)
        continue;
      if (m_test1C == 1)
      {
        kmfit->init();
        kmfit->setChisqCut(m_chisqMax);
        kmfit->AddTrack(0, 0.0, emcTrkk);
        kmfit->AddTrack(1, 0.0, emcTrkl);
        kmfit->AddResonance(0, 0.1349770, 0, 1);
        bool oksq = kmfit->Fit();

        if (oksq)
        {
          chi2.push_back(kmfit->chisq());
          //      HepLorentzVector pi0p4 = kmfit->pfit(0) + kmfit->pfit(1);
          pgam3_1C.push_back(kmfit->pfit(0));
          pgam4_1C.push_back(kmfit->pfit(1));
          igam3.push_back(iGam[k]);
          igam4.push_back(iGam[l]);
          pgam3.push_back(ptrkk);
          pgam4.push_back(ptrkl);
        }
      }
    }
  }
  if (m_debug)
    cout << __LINE__ << endl;
  int ngam1 = igam1.size();
  int ngam2 = igam2.size();
  int ngam3 = igam3.size();
  int ngam4 = igam4.size();
  if ((ngam3 != ngam4) || (ngam1 != ngam2))
    return StatusCode::SUCCESS;
  int ngam12 = ngam1;
  int ngam34 = ngam3;
  if (ngam12 == 0)
    return StatusCode::SUCCESS;
  if (fabs(signal) == 1)
    Ncut3++; //Ncut3 should equal Ncut2;
  if (ngam34 == 0)
    return StatusCode::SUCCESS;
  //cout<<"a="<<ngam<<endl;
  if (fabs(signal) == 1)
    Ncut4++;

  if (m_debug)
    cout << __LINE__ << endl;
  ///////////////////////
  //  // get beam energy and beta
  // ///////////////////////

  if (m_ReadBeamEFromDB)
  {
    if (m_usecalibBeamE)
      m_readDb.setcalib(true);
    m_beamE = m_readDb.getbeamE(m_run, m_beamE);
    if (m_run > 0)
      m_beta = m_readDb.getbeta();
  }
  double ebeam = m_beamE;
  double deltaE_mina = 9999;
  double deltaE_minb = 9999;
  double deltaE_minc = 9999;
  double deltaE_mind = 9999;
  double chisq = 0;
  int Km_index = -1;
  int Kp_index = -1;
  int p_index = -1;
  int pbar_index = -1;
  int gam1_index = -1;
  int gam3_index = -1;
  int gam1_indexm = -1;
  int gam3_indexm = -1;
  int gam2_index = -1;
  int gam4_index = -1;
  int gam2_indexm = -1;
  int gam4_indexm = -1;
  int pim_index = -1;
  int pip_index = -1;
  int a = 0;
  int b = 0;
  int c = 0;
  int d = 0;
  int num_cana = 0;
  int num_canb = 0;
  HepLorentzVector pKm_p4(0, 0, 0, 0), p_p4(0, 0, 0, 0), gam1a_p4(0, 0, 0, 0), gam2a_p4(0, 0, 0, 0), gam1b_p4(0, 0, 0, 0), gam2b_p4(0, 0, 0, 0), gam3a_p4(0, 0, 0, 0), gam4a_p4(0, 0, 0, 0), gam3b_p4(0, 0, 0, 0), gam4b_p4(0, 0, 0, 0), pipa_p4(0, 0, 0, 0), pipb_p4(0, 0, 0, 0), pgam1a_1C4p(0, 0, 0, 0), pgam2a_1C4p(0, 0, 0, 0), pgam1b_1C4p(0, 0, 0, 0), pgam2b_1C4p(0, 0, 0, 0), pgam3a_1C4p(0, 0, 0, 0), pgam4a_1C4p(0, 0, 0, 0), pgam3b_1C4p(0, 0, 0, 0), pgam4b_1C4p(0, 0, 0, 0), pKp_p4(0, 0, 0, 0), pbar_p4(0, 0, 0, 0), pima_p4(0, 0, 0, 0), pimb_p4(0, 0, 0, 0);

  for (int i = 0; i < npim; i++)
  {
    for (int j = 0; j < npbar; j++)
    {
      for (int k = 0; k < ngam34; k++)
      {
        for (int l = 0; l < npip; l++)
        {
          for (int m = 0; m < ngam12; m++)
          {
            HepLorentzVector psigma = ppbar[j] + pgam3_1C[k] + pgam4_1C[k];
            HepLorentzVector etap = ppim[i] + pgam1_1C[m] + pgam2_1C[m] + ppip[l];

            if (igam1[m] == igam3[k] || igam1[m] == igam4[k])
              continue;
            if (igam2[m] == igam3[k] || igam2[m] == igam4[k])
              continue;
            //					HepLorentzVector plambuda =  ppbar[j] + ppip[l];
            //					HepLorentzVector kshort =  pKm[i] + ppip[l];
            //					if(plambuda.m()<1.12&&plambuda.m()>1.11)continue;
            if(m_debug)
              cout<< __LINE__ << " " << " k " <<  k << " m "<< m << " psigma.m() "<<  psigma.m() << " etap.m() " << etap.m()<< endl;
            if (psigma.m() < 1.174 || psigma.m() > 1.2)
              continue;
            if (etap.m() < 0.76 || etap.m() > 0.8)
              continue;
            //					if(kshort.m()>0.48&&kshort.m()<0.52)continue;
            if (ipim[i] == ipbar[j])
              continue;
            if(m_debug)
              cout<< __LINE__ << " 0000000000"<< " " << " k " <<  k << " m "<< m << " psigma.m() "<<  psigma.m() << " etap.m() " << etap.m()<< endl;
            HepLorentzVector pLambda = ppim[i] + ppbar[j] + pgam3_1C[k] + pgam4_1C[k] + ppip[l] + pgam1_1C[m] + pgam2_1C[m];
            pLambda.boost(-m_beta);
            //double deltaE = fabs(pLambda.t() - ebeam);
            double deltaEb = pLambda.t() - ebeam;
            d++;
            if (m_debug)
              cout << __LINE__ << " deltaE_minb " << deltaE_minb <<  " deltaEb:" << deltaEb  << endl;
            if (fabs(deltaEb) < fabs(deltaE_minb))
            {
              b = 1;
              deltaE_minb = deltaEb;
              pipb_p4 = ppip[l];
              pbar_p4 = ppbar[j];
              gam3b_p4 = pgam3[k];
              gam4b_p4 = pgam4[k];
              gam1b_p4 = pgam1[m];
              gam2b_p4 = pgam2[m];
              pimb_p4 = ppim[i];
              pgam3b_1C4p = pgam3_1C[k];
              pgam4b_1C4p = pgam4_1C[k];
              pgam1b_1C4p = pgam1_1C[m];
              pgam2b_1C4p = pgam2_1C[m];
              num_canb = d;
            }
          }
        }
      }
    }
  }
  if (m_debug)
    cout << __LINE__ << endl;

  for (int i = 0; i < npip; i++)
  {
    for (int j = 0; j < np; j++)
    {
      for (int k = 0; k < ngam34; k++)
      {
        for (int l = 0; l < npim; l++)
        {
          for (int m = 0; m < ngam12; m++)
          {
            HepLorentzVector psigma = pp[j] + pgam3_1C[k] + pgam4_1C[k];
            HepLorentzVector etap = ppip[i] + pgam1_1C[m] + pgam2_1C[m] + ppim[l];
            if (igam1[m] == igam3[k] || igam1[m] == igam4[k])
              continue;
            if (igam2[m] == igam3[k] || igam2[m] == igam4[k])
              continue;
            //				HepLorentzVector plambuda =  pp[j] + ppim[l];
            //					HepLorentzVector kshort =  pKp[i] + ppim[l];

            //					if(plambuda.m()<1.12&&plambuda.m()>1.11)continue;
            if(m_debug)
              cout<< __LINE__ << " " << " k " <<  k << " m "<< m << " psigma.m() "<<  psigma.m() << " etap.m() " << etap.m()<< endl;
            if (psigma.m() < 1.174 || psigma.m() > 1.2)
              continue;
            if (etap.m() < 0.76 || etap.m() > 0.8)
              continue;
            //					if(kshort.m()>0.48&&kshort.m()<0.52)continue;
            if (ipip[i] == ip[j])
              continue;
            if(m_debug)
              cout<< __LINE__ << " 0000000000"<< " " << " k " <<  k << " m "<< m << " psigma.m() "<<  psigma.m() << " etap.m() " << etap.m()<< endl;
            HepLorentzVector pLambda = ppip[i] + pp[j] + pgam3_1C[k] + pgam4_1C[k] + ppim[l] + pgam1_1C[m] + pgam2_1C[m];
            pLambda.boost(-m_beta);
            //double deltaE = fabs(pLambda.t() - ebeam);
            double deltaEa = pLambda.t() - ebeam;
            c++;
            if (m_debug)
              cout << __LINE__ << " deltaE_mina " << deltaE_mina <<  " deltaEa:" << deltaE_mina  << endl;
            if (fabs(deltaEa) < fabs(deltaE_mina))
            {
              a = 1;
              deltaE_mina = deltaEa;
              pipa_p4 = ppip[i];
              p_p4 = pp[j];
              gam3a_p4 = pgam3[k];
              gam4a_p4 = pgam4[k];
              gam1a_p4 = pgam1[m];
              gam2a_p4 = pgam2[m];
              pima_p4 = ppim[l];
              pgam3a_1C4p = pgam3_1C[k];
              pgam4a_1C4p = pgam4_1C[k];
              pgam1a_1C4p = pgam1_1C[m];
              pgam2a_1C4p = pgam2_1C[m];
              num_cana = c;
            }
          }
        }
      }
    }
  }

  if (m_debug)
    cout << __LINE__ << endl;
  //cout<<"a="<<a<<endl;
  if (a == 0 && b == 0)
    return StatusCode::SUCCESS;
  if (m_debug)
    cout << __LINE__ << " start write" << endl;
  if (b == 1)
  {

    m_mode1 = mm_mode1;
    m_mode2 = mm_mode2;
    m_mode3 = mm_mode3;
    m_idxmc = numParticle;
    for (int i = 0; i < numParticle; i++)
    {
      m_pdgid[i] = M_pdgid[i];
      m_motheridx[i] = M_motheridx[i];
    }

    m_ndaughterAp = ndaughterAp;
    for (int aa = 0; aa < ndaughterAp; aa++)
      m_Ap_id[aa] = Ap_id[aa];
    for (int aa = 0; aa < ndaughterAp; aa++)
      for (int ll = 0; ll < 4; ll++)
        m_Ap_ptruth[aa][ll] = Ap_ptruth[aa][ll];

    m_ndaughterAm = ndaughterAm;
    for (int aa = 0; aa < ndaughterAm; aa++)
      m_Am_id[aa] = Am_id[aa];
    for (int aa = 0; aa < ndaughterAm; aa++)
      for (int ll = 0; ll < 4; ll++)
        m_Am_ptruth[aa][ll] = Am_ptruth[aa][ll];

    m_mcparticle_p = numParticle_p;
    m_mcparticle_m = numParticle_m;
    for (int i = 0; i < numParticle_p; i++)
    {
      m_pdgid_p[i] = M_pdgid_p[i];
      m_motheridx_p[i] = M_motheridx_p[i];
    }
    for (int i = 0; i < numParticle_m; i++)
    {
      m_pdgid_m[i] = M_pdgid_m[i];
      m_motheridx_m[i] = M_motheridx_m[i];
    }

    m_all = all;
    m_run = eventHeader->runNumber();
    m_event = eventHeader->eventNumber();
    m_signal = num_canb;
    //  m_signal =signal;
    m_bg = bg;
    for (int jj = 0; jj < 4; jj++)
      m_pbar_p4[jj] = pbar_p4[jj];
    for (int jj = 0; jj < 4; jj++)
      m_pim_p4[jj] = pimb_p4[jj];
    for (int jj = 0; jj < 4; jj++)
      m_pip_p4[jj] = pipb_p4[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam3b_p4[jj] = gam3b_p4[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam4b_p4[jj] = gam4b_p4[jj];

    m_pbarindex = pbar_index;
    m_p4index = 4;
    m_pi0m1c = (pgam3b_1C4p + pgam4b_1C4p).m();
    m_etam1c = (pgam1b_1C4p + pgam2b_1C4p).m();
    m_numothertrackp = ntrackp - 1;
    m_numothertrackm = ntrackm - 1;
    m_pi0m = (gam3b_p4 + gam4b_p4).m();
    m_etam = (gam1b_p4 + gam2b_p4).m();
    m_etaprimem = (pgam1b_1C4p + pgam2b_1C4p + pimb_p4 + pipb_p4).m();
    m_ebeam = ebeam;
    m_Sigmam = (pgam3b_1C4p + pgam4b_1C4p + pbar_p4).m();
    //	cout<<"deltaE="<<deltaE_minb<<endl;
    m_deltaE_min = deltaE_minb;
    //	cout<<"deltaE_min="<<m_deltaE_min<<endl;
    HepLorentzVector pLambda = pipb_p4 + pbar_p4 + pgam3b_1C4p + pgam4b_1C4p + pimb_p4 + pgam1b_1C4p + pgam2b_1C4p;

    pLambda.boost(-m_beta);
    if(m_debug) cout<< __LINE__ << " pLambda " << pLambda.px()<< " " << pLambda.py() << " "<< pLambda.pz() << " "<< pLambda.e()<< endl;
    double mbc2 = ebeam * ebeam - pLambda.v().mag2();
    m_bc = mbc2 > 0 ? sqrt(mbc2) : -10;
    m_rightflag = 2;
    m_np = np;
    m_npbar = npbar;
    if(m_debug)
    {
      cout << __LINE__ << "my raw" << endl;      
      HepLorentzVector etap4=(gam1b_p4 + gam2b_p4); cout << __LINE__ << " etap4 " << etap4.px()<< " " << etap4.py() << " "<< etap4.pz() << " "<< etap4.e() << " " << etap4.m() << endl;
      HepLorentzVector pi0p4=(gam3b_p4 + gam4b_p4); cout << __LINE__ << " pi0p4 " << pi0p4.px()<< " " << pi0p4.py() << " "<< pi0p4.pz() << " "<< pi0p4.e() << " " << pi0p4.m() << endl;
      cout << __LINE__ << " pimb_p4 " << pimb_p4.px()<< " " << pimb_p4.py() << " "<< pimb_p4.pz() << " "<< pimb_p4.e() << " " << pimb_p4.m() << endl;
      cout << __LINE__ << " pipb_p4 " << pipb_p4.px()<< " " << pipb_p4.py() << " "<< pipb_p4.pz() << " "<< pipb_p4.e() << " " << pipb_p4.m() << endl;
     
      cout << __LINE__ << "after fit" << endl;      
      HepLorentzVector etap41c=(pgam1b_1C4p + pgam2b_1C4p); cout << __LINE__ << " etap41c " << etap41c.px()<< " " << etap41c.py() << " "<< etap41c.pz() << " "<< etap41c.e() << " " << etap41c.m() << endl;
      HepLorentzVector pi0p41c=(pgam3b_1C4p + pgam4b_1C4p); cout << __LINE__ << " pi0p41c " << pi0p41c.px()<< " " << pi0p41c.py() << " "<< pi0p41c.pz() << " "<< pi0p41c.e() << " " << pi0p41c.m() << endl;
      HepLorentzVector etaprimep41c=(pgam1b_1C4p + pgam2b_1C4p + pimb_p4 + pipb_p4); cout << __LINE__ << " etaprimep41c" << etaprimep41c.px()<< " " << etaprimep41c.py() << " "<< etaprimep41c.pz() << " "<< etaprimep41c.e() << " " << etaprimep41c.m() << endl;
      HepLorentzVector sigmaa1c=(pgam3b_1C4p + pgam4b_1C4p + pbar_p4); cout << __LINE__ << " sigmaa1c" << sigmaa1c.px()<< " " << sigmaa1c.py() << " "<< sigmaa1c.pz() << " "<< sigmaa1c.e() << " " << sigmaa1c.m()<< endl;

      cout << __LINE__ << " b==1, -" << " m_pi0m1c  " << m_pi0m1c << "  m_etam1c " << m_etam1c << "  m_Sigmam " << m_Sigmam << "  m_etaprimem " << m_etaprimem << " m_deltaE_min  " << m_deltaE_min << "  m_bc " << m_bc << endl;
    }
    m_tuple1->write();
  }
  if (a == 1)
  {

    m_mode1 = mm_mode1;
    m_mode2 = mm_mode2;
    m_mode3 = mm_mode3;
    m_idxmc = numParticle;
    for (int i = 0; i < numParticle; i++)
    {
      m_pdgid[i] = M_pdgid[i];
      m_motheridx[i] = M_motheridx[i];
    }

    m_ndaughterAp = ndaughterAp;
    for (int aa = 0; aa < ndaughterAp; aa++)
      m_Ap_id[aa] = Ap_id[aa];
    for (int aa = 0; aa < ndaughterAp; aa++)
      for (int ll = 0; ll < 4; ll++)
        m_Ap_ptruth[aa][ll] = Ap_ptruth[aa][ll];

    m_ndaughterAm = ndaughterAm;
    for (int aa = 0; aa < ndaughterAm; aa++)
      m_Am_id[aa] = Am_id[aa];
    for (int aa = 0; aa < ndaughterAm; aa++)
      for (int ll = 0; ll < 4; ll++)
        m_Am_ptruth[aa][ll] = Am_ptruth[aa][ll];

    m_mcparticle_p = numParticle_p;
    m_mcparticle_m = numParticle_m;
    for (int i = 0; i < numParticle_p; i++)
    {
      m_pdgid_p[i] = M_pdgid_p[i];
      m_motheridx_p[i] = M_motheridx_p[i];
    }
    for (int i = 0; i < numParticle_m; i++)
    {
      m_pdgid_m[i] = M_pdgid_m[i];
      m_motheridx_m[i] = M_motheridx_m[i];
    }
    m_all = all;
    m_run = eventHeader->runNumber();
    m_event = eventHeader->eventNumber();
    // m_signal =signal;
    m_signal = num_cana;
    m_bg = bg;
    for (int jj = 0; jj < 4; jj++)
      m_p_p4[jj] = p_p4[jj];
    for (int jj = 0; jj < 4; jj++)
      m_pim_p4[jj] = pima_p4[jj];
    for (int jj = 0; jj < 4; jj++)
      m_pip_p4[jj] = pipa_p4[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam3a_p4[jj] = gam3a_p4[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam4a_p4[jj] = gam4a_p4[jj];

    m_pindex = p_index;
    m_p4index = 4;
    m_chi2 = chisq;
    m_pi0m1c = (pgam3a_1C4p + pgam4a_1C4p).m();
    m_etam1c = (pgam1a_1C4p + pgam2a_1C4p).m();
    m_numothertrackp = ntrackp - 1;
    m_numothertrackm = ntrackm - 1;
    m_pi0m = (gam3a_p4 + gam4a_p4).m();
    m_etam = (gam1a_p4 + gam2a_p4).m();
    m_etaprimem = (pgam1a_1C4p + pgam2a_1C4p + pima_p4 + pipa_p4).m();
    m_ebeam = ebeam;
    m_Sigmam = (pgam3a_1C4p + pgam4a_1C4p + p_p4).m();
    m_deltaE_min = deltaE_mina;
    //	cout<<"deltaE="<<deltaE_mina<<endl;
    HepLorentzVector pLambda = pipa_p4 + p_p4 + pgam3a_1C4p + pgam4a_1C4p + pima_p4 + pgam1a_1C4p + pgam2a_1C4p;
    pLambda.boost(-m_beta);
    if(m_debug) cout<< __LINE__ << " pLambda " << pLambda.px()<< " " << pLambda.py() << " "<< pLambda.pz() << " "<< pLambda.e()<< endl;
    double mbc2 = ebeam * ebeam - pLambda.v().mag2();
    m_bc = mbc2 > 0 ? sqrt(mbc2) : -10;
    m_rightflag = 1;
    m_np = np;
    m_npbar = npbar;
    if(m_debug)
    {
      cout << __LINE__ << "my raw" << endl;            
      HepLorentzVector etap4=(gam1a_p4 + gam2a_p4); cout << __LINE__ << " etap4 " << etap4.px()<< " " << etap4.py() << " "<< etap4.pz() << " "<< etap4.e() << " " << etap4.m() << endl;
      HepLorentzVector pi0p4=(gam3a_p4 + gam4a_p4); cout << __LINE__ << " pi0p4 " << pi0p4.px()<< " " << pi0p4.py() << " "<< pi0p4.pz() << " "<< pi0p4.e() << " " << pi0p4.m() << endl;
      cout << __LINE__ << " pima_p4 " << pima_p4.px()<< " " << pima_p4.py() << " "<< pima_p4.pz() << " "<< pima_p4.e() << " " << pima_p4.m() << endl;
      cout << __LINE__ << " pipa_p4 " << pipa_p4.px()<< " " << pipa_p4.py() << " "<< pipa_p4.pz() << " "<< pipa_p4.e() << " " << pipa_p4.m() << endl;
      
      cout << __LINE__ << "after fit" << endl ;     
      HepLorentzVector etap41c=(pgam1a_1C4p + pgam2a_1C4p); cout << __LINE__ << " etap41c " << etap41c.px()<< " " << etap41c.py() << " "<< etap41c.pz() << " "<< etap41c.e() << " " << etap41c.m() << endl;
      HepLorentzVector pi0p41c=(pgam3a_1C4p + pgam4a_1C4p); cout << __LINE__ << " pi0p41c " << pi0p41c.px()<< " " << pi0p41c.py() << " "<< pi0p41c.pz() << " "<< pi0p41c.e() << " " << pi0p41c.m() << endl;
      HepLorentzVector etaprimep41c=(pgam1a_1C4p + pgam2a_1C4p + pima_p4 + pipa_p4); cout << __LINE__ << " etaprimep41c" << etaprimep41c.px()<< " " << etaprimep41c.py() << " "<< etaprimep41c.pz() << " "<< etaprimep41c.e() << " " << etaprimep41c.m() << endl;
      HepLorentzVector sigmaa1c=(pgam3a_1C4p + pgam4a_1C4p + p_p4); cout << __LINE__ << " sigmaa1c" << sigmaa1c.px()<< " " << sigmaa1c.py() << " "<< sigmaa1c.pz() << " "<< sigmaa1c.e() << " " << sigmaa1c.m() << endl;

      cout << __LINE__ << " a==1, -" << " m_pi0m1c  " << m_pi0m1c << "  m_etam1c " << m_etam1c << "  m_Sigmam " << m_Sigmam << "  m_etaprimem " << m_etaprimem << " m_deltaE_min  " << m_deltaE_min << "  m_bc " << m_bc << endl;
    }
    m_tuple1->write();
  }
  Ncut5++;
  if(m_debug)
  {  
    cout << "Ntotal  " << Ntotal << endl;
    cout << "Ncut0   " << Ncut0 << endl;
    cout << "Ncut1   " << Ncut1 << endl;
    cout << "Ncut2   " << Ncut2 << endl;
    cout << "Ncut3   " << Ncut3 << endl;
    cout << "Ncut4   " << Ncut4 << endl;
    cout << "Ncut5   " << Ncut5 << endl;
    cout << "all       " << all << endl;
    cout << "H       " << H << endl;
  }
  return StatusCode::SUCCESS;
}
StatusCode LambdacAlg::endRun()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "LambdacAlg::endRun()" << endreq;
  //add your code here
  return StatusCode::SUCCESS;
}
StatusCode LambdacAlg::finalize()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "LambdacAlg::finalize()" << endreq;
  //add your code here
  cout << "attention: if Ncut2!= Ncut3 ,you should check" << endl;
  cout << "Ntotal  " << Ntotal << endl;
  cout << "Ncut0   " << Ncut0 << endl;
  cout << "Ncut1   " << Ncut1 << endl;
  cout << "Ncut2   " << Ncut2 << endl;
  cout << "Ncut3   " << Ncut3 << endl;
  cout << "Ncut4   " << Ncut4 << endl;
  cout << "Ncut5   " << Ncut5 << endl;
  cout << "all       " << all << endl;
  cout << "H       " << H << endl;
  //	cout << "Proton       " <<Proton << endl;
  //	cout << "al       " <<al << endl;
  //	cout << "yes     " << yes << endl;
  //	cout << "Npp>0:  " << Npp << endl;
  //	cout << "Npm>0:  " << Npm << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  cout << "-------------------------------           -------------------------------" << endl;
  cout << "--------------------                                ---------------------" << endl;
  cout << "-------------  omega modify zhou to me code, v100 ------------------" << endl;
  cout << "--------------------                               ----------------------" << endl;
  cout << "------------------------------           --------------------------------" << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  return StatusCode::SUCCESS;
}

//add your code here,for other member-functions
