// ucas-PRLi (lipeirong11@mails.ucas.ac.cn)
#include "LambdacAlg/LambdacAlg.h"
#include "BestDTagSvc/BestDTagSvc.h"
#include "McDecayModeSvc/McDecayModeSvc.h"

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"

#include "DstEvent/TofHitStatus.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "TMath.h"

#include "EvTimeEvent/RecEsTime.h"
#include "McTruth/McParticle.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;

#include "CLHEP/Geometry/Point3D.h"
#include "DTagTool/DTagTool.h"
#include "ParticleID/ParticleID.h"
// #include "SimplePIDSvc/ISimplePIDSvc.h"
#include "LambdacAlg/MyParticle.h"
#include "LambdacAlg/MyPid.h"

#include "VertexFit/Helix.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"

#define RED "\033[31m"     /* Red */
#define GREEN "\033[32m"   /* Green */
#define YELLOW "\033[33m"  /* Yellow */
#define BLUE "\033[34m"    /* Blue */
#define ENDCOLOR "\033[0m" /* ENDCOLOR */

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include <vector>
const double xmass[5] = {0.000511, 0.105658, 0.139571, 0.493677, 0.938272};
typedef std::vector<double> Vdouble;
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;

int Ntotal(0), Ncut0(0), Ncut1(0), Ncut2(0), Ncut3(0), Ncut4(0), Ncut5(0), Ncut6(0), Ncut7(0), H(0), A(0), all(0),
    all_m(0), all_p(0), al(0);

int N1(0), N2(0), N3(0), N4(0);

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
  declareProperty("Costheta", m_costheta = 0.93);
  declareProperty("Vr0cut", m_vr0cut = 1.0);
  declareProperty("Vz0cut", m_vz0cut = 10.0);
  declareProperty("Vz0cut1", m_vz0cut1 = 20.0);
  declareProperty("ChisqMax", m_chisqMax = 200);

  declareProperty("SkimFlag", m_skim = false);
  declareProperty("PhotonMinEnergy", m_minEnergy = 0.025);
  declareProperty("GammaAngleCut", m_gammaAngleCut = 10.0);
  declareProperty("Test1C", m_test1C = 1); // ec fit
  declareProperty("GammathCut", m_gammathCut = 14.0);
  declareProperty("GammatlCut", m_gammatlCut = 0.0);
  declareProperty("PhotonMaxCosThetaBarrel", m_maxCosThetaBarrel = 0.8);
  declareProperty("PhotonMinCosThetaEndcap", m_minCosThetaEndcap = 0.86);
  declareProperty("PhotonMaxCosThetaEndcap", m_maxCosThetaEndcap = 0.92);
  declareProperty("PhotonMinEndcapEnergy", m_minEndcapEnergy = 0.050);

  declareProperty("EtaMinMass", m_EtaMinMass = 0.50);
  declareProperty("EtaMaxMass", m_EtaMaxMass = 0.56);
  declareProperty("Pi0MinMass", m_Pi0MinMass = 0.115);
  declareProperty("Pi0MaxMass", m_Pi0MaxMass = 0.15);
  declareProperty("SigmaMinMass", m_SigmaMinMass = 1.174);
  declareProperty("SigmaMaxMass", m_SigmaMaxMass = 1.2);

  declareProperty("EtaPrimeMinMass", m_EtaPrimeMinMass = 0.946);
  declareProperty("EtaPrimeMaxMass", m_EtaPrimeMaxMass = 0.968);

  declareProperty("Debug", m_debug = false);
  declareProperty("Isqqbar", m_isqqbar = false);

  declareProperty("BeamE", m_beamE = 2.313);
  declareProperty("ReadBeamEFromDB", m_ReadBeamEFromDB = false);
  declareProperty("UseCalibBeamE", m_usecalibBeamE = false);
  declareProperty("CheckTotal", m_checktotal = true);
}
LambdacAlg::~LambdacAlg()
{
  // add your code for deconstructor
}
StatusCode LambdacAlg::initialize()
{
  if (m_debug)
    cout << __LINE__ << " debug begin initialize " << endl;
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "LambdacAlg::initialize()" << endreq;
  m_beta.setX(0.011);
  m_beta.setY(0);
  m_beta.setZ(0);
  // add your code here
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
      // status = m_tuple1->addItem("yes", m_yes);
      // status = m_tuple1->addItem("no", m_no);

      status = m_tuple1->addItem("indexmc", m_idxmc, 0, 100);
      status = m_tuple1->addIndexedItem("pdgid", m_idxmc, m_pdgid);
      status = m_tuple1->addIndexedItem("motheridx", m_idxmc, m_motheridx);
      status = m_tuple1->addItem("p4index", m_p4index, 0, 10);
      // status = m_tuple1->addIndexedItem("Km", m_p4index, m_Km_p4);
      // status = m_tuple1->addIndexedItem("Kp", m_p4index, m_Kp_p4);
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
      status = m_tuple1->addIndexedItem("gamb", m_p4index, m_gamb_p4);
      status = m_tuple1->addIndexedItem("gama", m_p4index, m_gama_p4);
      // 1c __________________________________________________________________
      status = m_tuple1->addIndexedItem("gam1_1c", m_p4index, m_gam1_p4_1c);
      status = m_tuple1->addIndexedItem("gam2_1c", m_p4index, m_gam2_p4_1c);
      status = m_tuple1->addIndexedItem("gam3_1c", m_p4index, m_gam3_p4_1c);
      status = m_tuple1->addIndexedItem("gam4_1c", m_p4index, m_gam4_p4_1c);
      status = m_tuple1->addItem("sigmam1c", m_Sigmam1c);
      status = m_tuple1->addItem("etaprimem1c", m_etaprimem1c);

      // R1c _________________________________________________________________
      status = m_tuple1->addIndexedItem("pall", m_p4index, m_pall_p4);
      // for pi
      status = m_tuple1->addIndexedItem("p_pim", m_p4index, m_pim_p4_r1c);
      status = m_tuple1->addIndexedItem("p_pip", m_p4index, m_pip_p4_r1c);

      status = m_tuple1->addIndexedItem("gam1", m_p4index, m_gam1_p4);
      status = m_tuple1->addIndexedItem("gam2", m_p4index, m_gam2_p4);
      status = m_tuple1->addIndexedItem("gam3", m_p4index, m_gam3_p4);
      status = m_tuple1->addIndexedItem("gam4", m_p4index, m_gam4_p4);

      status = m_tuple1->addIndexedItem("pall_r1c", m_p4index, m_pall_p4_r1c);
      status = m_tuple1->addIndexedItem("gam1_r1c", m_p4index, m_gam1_p4_r1c);
      status = m_tuple1->addIndexedItem("gam2_r1c", m_p4index, m_gam2_p4_r1c);
      status = m_tuple1->addIndexedItem("gam3_r1c", m_p4index, m_gam3_p4_r1c);
      status = m_tuple1->addIndexedItem("gam4_r1c", m_p4index, m_gam4_p4_r1c);
      status = m_tuple1->addItem("chi2_min_r1c", m_chi2_min_r1c);

      status = m_tuple1->addItem("pi0mR1c", m_pi0mR1c);
      status = m_tuple1->addItem("etamR1c", m_etamR1c);
      status = m_tuple1->addItem("sigmamR1c", m_SigmamR1c);
      status = m_tuple1->addItem("pcharge", m_pcharge);
      // 2c _________________________________________________________________



      // ____________________________________________________________________
      status = m_tuple1->addItem("pi0m", m_pi0m);
      status = m_tuple1->addItem("etam", m_etam);
      status = m_tuple1->addItem("lambda", m_lambda);
      status = m_tuple1->addItem("ksi", m_ksi);
      status = m_tuple1->addItem("kstar", m_kstar);
      status = m_tuple1->addItem("sigmastar", m_sigmastar);
      status = m_tuple1->addItem("etaprimem", m_etaprimem);
      status = m_tuple1->addItem("all", m_all);
      status = m_tuple1->addItem("all_m", m_all_m);
      status = m_tuple1->addItem("all_p", m_all_p);
      status = m_tuple1->addItem("Kmindex", m_Kmindex);
      status = m_tuple1->addItem("pbarindex", m_pbarindex);
      status = m_tuple1->addItem("pindex", m_pindex);
      status = m_tuple1->addItem("Kpindex", m_Kpindex);
      status = m_tuple1->addItem("r", m_r);
      status = m_tuple1->addItem("pi0m1c", m_pi0m1c);
      status = m_tuple1->addItem("etam1c", m_etam1c);

      status = m_tuple1->addItem("sigmam", m_Sigmam);
      // status = m_tuple1->addItem("num_othertrackp", m_numothertrackp);
      // status = m_tuple1->addItem("num_othertrackm", m_numothertrackm);
      status = m_tuple1->addItem("E_beam", m_ebeam);
      status = m_tuple1->addItem("deltaE_min", m_deltaE_min);
      status = m_tuple1->addItem("M_BC", m_bc);
      status = m_tuple1->addItem("np", m_np);
      status = m_tuple1->addItem("npbar", m_npbar);
      // status = m_tuple1->addItem("eop_pim", m_eop_pim);
      // status = m_tuple1->addItem("eop_pip", m_eop_pip);
      // status = m_tuple1->addItem("eop_p", m_eop_p);
      // status = m_tuple1->addItem("eop_pbar", m_eop_pbar);
      // status = m_tuple1->addItem("eop_Km", m_eop_Km);
      // status = m_tuple1->addItem("eop_Kp", m_eop_Kp);
      // status = m_tuple1->addItem("ngam", m_ngam);
      // status = m_tuple1->addItem("chi2", m_chi2);

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
  // add your code here
  return StatusCode::SUCCESS;
}
StatusCode LambdacAlg::execute()
{
  m_rightflag = -1;
  int signal = -9999;
  int bg = -1;
  int yes = -1;
  int no = -1;
  Ntotal++;

  if (m_debug)
    cout << RED << "m_debug3 begin execute, Ntotal(event): " << Ntotal << "\033[0m" << endl;

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

  // if ((m_run == 35412 || m_run == 35413 || m_run == 35414 || m_run == 37373 || m_run == 37733 || m_run == 37734 ||
  //      m_run == 37736 || m_run == 37738 || m_run == 37740 || m_run == 37741))
  //   return StatusCode::SUCCESS;
  // if(!(m_run==35966&&(m_event==10304||m_event==77483||m_event==79947||m_event==97315)))   return StatusCode::SUCCESS;
  // cout<<endl<<endl<<"********************************************************************************************************"<<endl;

  log << MSG::DEBUG << "run, evtnum = " << m_run << " , " << m_event << endreq;

  IMcDecayModeSvc *i_svc;
  StatusCode sc_DecayModeSvc = service("McDecayModeSvc", i_svc);
  if (sc_DecayModeSvc.isFailure())
  {
    log << MSG::FATAL << "Could not load McDecayModeSvc!" << endreq;
    return sc_DecayModeSvc;
  }
  m_svc = dynamic_cast<McDecayModeSvc *>(i_svc);

  // int protons = 0;
  // int pms = 0;
  // int pi0s = 0;
  // int pi0mis = 0;
  // int gams = 0;
  // int gamms = 0;
  // int gam = 0;
  // int gamm = 0;
  // int pim = 0;
  // int piplus = 0;
  // int kplus = 0;
  // int km = 0;
  // int pi0 = 0;
  // int pi0mi = 0;
  // int proton = 0;
  // int pm = 0;
  // int sigm = 0;
  // int sigmm = 0;
  // int kstarp = 0;
  // int kstarm = 0;

  int numParticle_p = 0;
  int numParticle_m = 0;
  int numParticle = 0;
  int M_pdgid_p[100];
  int M_motheridx_p[2000];

  int M_pdgid_m[100];
  int M_motheridx_m[2000];

  int M_pdgid[100];
  int M_motheridx[2000];

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
  // truth
  if (m_debug)
    std::cout << "mode1=" << mm_mode1 << ", mode2=" << mm_mode2 << ", mode3=" << mm_mode3 << ", runNo" << runNo
              << ", eventNo " << eventNo << std::endl;

#pragma region for_mc______________________________________________________________
  if (eventHeader->runNumber() < 0)
  {
    SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
    Vint pdgid;

    Vint motherindex;
    pdgid.clear();
    motherindex.clear();

    if (!mcParticleCol)
    {
      std::cout << __LINE__ <<  "Could not retrieve McParticelCol" << std::endl;
      // return StatusCode::FAILURE;
    }
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
          cout << __LINE__ <<  "mcParticleCol pdg: " << pdg << ", motherpdg: " << motherpdg << ", mmotherpdg:" << mmotherpdg << endl;

        // truth for lambda_c+ to ... ------------------------------------
        if ((*iter_mc)->particleProperty() == 4122)
        {// lambdac+ -> sigma+ eta', store sigma+ eta'
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ <<  " ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

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
              cout << __LINE__ <<  " ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

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
              cout << __LINE__ <<  " ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          } // End of "gc.size() > 0" IF
        }
        if (pdg == 331 && motherpdg == 4122)
        { //// eta'->pi+ pi- eta, store pi+ pi- eta
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ <<  " ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          }
        }

        if (pdg == 221 && motherpdg == 331 && mmotherpdg == 4122)
        { // eta->gam gam //store 2 gamma
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ <<  " ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          } // End of "gc.size() > 0" IF
        }      

        // truth for lambda_c- to ... ------------------------------------
        if ((*iter_mc)->particleProperty() == -4122)
        {// lambdac- -> sigma- eta', store sigma- eta'
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();

            if (m_debug)
              cout << __LINE__ <<  " ndaughterAm[" << ndaughterAm << "]: " << Am_id[ndaughterAm] << endl;

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
              cout << __LINE__ <<  " ndaughterAm[" << ndaughterAm << "]: " << Am_id[ndaughterAm] << endl;

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
              cout << __LINE__ <<  " ndaughterAm[" << ndaughterAm << "]: " << Am_id[ndaughterAm] << endl;

            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }
        if (pdg == 331 && motherpdg == -4122)
        { //// eta'->pi+ pi- eta, store pi+ pi- eta
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

        if (pdg == 221 && motherpdg == 331 && mmotherpdg == -4122)
        { // eta->gam gam //store 2 gamma
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
            cout << __LINE__ << __LINE__ << "+iter_mc: " << (*iter_mc)->particleProperty() << ", numParticle: " << numParticle
                 << endl;
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
            cout << __LINE__ << "-iter_mc: " << (*iter_mc)->particleProperty() << ", numParticle: " << numParticle << endl;
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

      if (ndaughterAp == 11 && Ap_id[0] == 3222 && Ap_id[1] == 331 && Ap_id[2] == 2212 &&
          Ap_id[3] == 111 && Ap_id[4] == 211 && Ap_id[5] == -211 && Ap_id[6] == 221 &&
          Ap_id[7] == 22 && Ap_id[8]==22 && Ap_id[9] == 22 && Ap_id[10] == 22 )
      {
        signal = 1;
        all++;
        all_p++;
      }
      if (ndaughterAm == 11 && Am_id[0] == -3222 && Am_id[1] == 331 && Am_id[2] == -2212 &&
          Am_id[3] == 111 && Am_id[4] == 211 && Am_id[5] == -211 && Am_id[6] == 221 &&
          Am_id[7] == 22 && Am_id[8]==22 && Am_id[9] == 22 && Am_id[10] == 22 )
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

  // write truth
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
        m_Ap_ptruth_[aa][ll] = Ap_ptruth[aa][ll];

    m_ndaughterAm_ = ndaughterAm;
    for (int aa = 0; aa < ndaughterAm; aa++)
      m_Am_id_[aa] = Am_id[aa];

    for (int aa = 0; aa < ndaughterAm; aa++)
      for (int ll = 0; ll < 4; ll++)
        m_Am_ptruth_[aa][ll] = Am_ptruth[aa][ll];

    m_tuple2->write();
  }

  // end mc;
#pragma endregion

#pragma region for_track______________________________________________________________
  // good track  =============
  if (m_debug)
    cout << "\n" << __LINE__ << " begin choose good track" << endl;

  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);

  Vint goodTrack;
  goodTrack.clear();

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
    cout << __LINE__ << " choose good track" << endl;
  for (int i = 0; i < evtRecEvent->totalCharged(); i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
    if (!(*itTrk)->isMdcTrackValid())
      continue;
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();

    double pch = mdcTrk->p();
    HepVector a = mdcTrk->helix();
    HepSymMatrix Ea = mdcTrk->err();
    HepPoint3D point0(0., 0., 0.); // the initial point for MDC recosntruction
    HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
    VFHelix helixip(point0, a, Ea);
    helixip.pivot(IP);
    HepVector vecipa = helixip.a();
    double Rvxy0 = fabs(vecipa[0]); // the nearest distance to IP in xy plane
    double Rvz0 = vecipa[3];        // the nearest distance to IP in z direction
    double Rvphi0 = vecipa[1];
    double costheta = cos(mdcTrk->theta());

    if (costheta >= m_costheta)
      continue;
    if (fabs(Rvz0) >= m_vz0cut1)
      continue;
    if (fabs(Rvz0) >= m_vz0cut)
      continue;
    if (fabs(Rvxy0) >= m_vr0cut)
      continue;

    goodTrack.push_back(i);
    nCharge += mdcTrk->charge();
  }
  // Finish Good Charged Track SKction
  if (goodTrack.size() < 3)
  {
    cout << __LINE__ << "return StatusCode::SUCCESS; goodTrack.size() < 3" << endl;
    return StatusCode::SUCCESS;
  }
  if (abs(signal) == 1)
    Ncut0++;

  Vint trackProntonP, trackProntonPbar;
  std::vector<MyParticle> proton, piPlus, piMin;
  proton.clear();
  piPlus.clear();
  piMin.clear();
  trackProntonP.clear();
  trackProntonPbar.clear();
  for (int i = 0; i < goodTrack.size(); i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + goodTrack[i];
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();

    RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
    MyPid *pid = new MyPid(*itTrk);
    
    if (pid->isproton())
    {
      if (mdcTrk->charge() == 1)
        trackProntonP.push_back(goodTrack[i]);
      if (mdcTrk->charge() == -1)
        trackProntonPbar.push_back(goodTrack[i]);

      mdcKalTrk->setPidType(RecMdcKalTrack::proton);
      HepLorentzVector p4 = mdcKalTrk->p4(xmass[4]);
      WTrackParameter wtrkp(xmass[4], mdcKalTrk->getZHelixP(), mdcKalTrk->getZErrorP());

      if (m_debug)
        cout << __LINE__ << " i " << i << " p4.m() " << p4.m() << endl;
      MyParticle tmp(goodTrack[i], p4, mdcTrk->charge(), wtrkp);
      // tmp.chi = pid->getChi();
      proton.push_back(tmp);
    }

    if (pid->ispion())
    {
      cout << __LINE__ << endl;
      // RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
      mdcKalTrk->setPidType(RecMdcKalTrack::pion);
      HepLorentzVector p4 = mdcKalTrk->p4(xmass[2]);
      WTrackParameter wtrkp(xmass[2], mdcKalTrk->getZHelixP(), mdcKalTrk->getZErrorP());

      if (mdcTrk->charge() == 1)
      {
        MyParticle tmp(goodTrack[i], p4, mdcTrk->charge(), wtrkp);
        tmp.chi = pid->getChi();
        piPlus.push_back(tmp);
        // piPlus.push_back(goodTrack[i]);
      }
        
      if (mdcTrk->charge() == -1)
      {
        MyParticle tmp(goodTrack[i], p4, mdcTrk->charge(), wtrkp);
        tmp.chi = pid->getChi();
        piMin.push_back(tmp);
        // piMin.push_back(goodTrack[i]);
      }
    }
  }

  if(m_debug)
  {
    cout << __LINE__ << " piMin.size() " << piMin.size() << " piPlus.size() " << piPlus.size() << " proton.size() " <<  proton.size() << endl;
  }
  if(piMin.size()<1 ||  piPlus.size()<1 || proton.size() <1)
  {   
    if(m_debug)
    {
      // cout << piMin.size() << " " << piPlus.size() << " " <<  proton.size() << endl;
      cout << __LINE__ << "return StatusCode::SUCCESS; piMin.size()<1 ||  piPlus.size()<1 || proton.size() <1" << endl;
    }
    return StatusCode::SUCCESS;
  }
  if (abs(signal) == 1)
    Ncut1++;

  int np = trackProntonP.size();
  int npbar = trackProntonPbar.size();
  // m_np = np;
  // m_npbar = npbar;

  // if (m_isqqbar)
  // {
    //if (np == 0 && npbar == 0)
      //return StatusCode::SUCCESS;
  // }


#pragma endregion

#pragma region for_shower______________________________________________________________
  // good shower ======================
  std::vector<MyParticle> emcGamma;
  emcGamma.clear();
  if (m_debug)
    cout << __LINE__ << " begin choose good gamma" << endl;
  for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
    if (!(*itTrk)->isEmcShowerValid())
      continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    Hep3Vector emcp(emcTrk->x(), emcTrk->y(), emcTrk->z());
    HepLorentzVector shP4 = getP4(emcTrk, xorigin);
    double cosThetaSh = shP4.vect().cosTheta();
    // double dthe = 200.;
    // double dphi = 200.;
    double dang = 200.;
    if (m_debug)
      cout << __LINE__ << " choose good gamma" << endl;

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
      // double thed = extp.theta() - emcp.theta();
      // double phid = extp.deltaPhi(emcp);
      // thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      // phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
      if (angd < dang)
      {
        dang = angd;
        // dthe = thed;
        // dphi = phid;
      }
    }
    if (dang >= 200)
      continue;

    double eraw = emcTrk->energy();
    double getTime = emcTrk->time();
    // dthe = dthe * 180 / (CLHEP::pi);
    // dphi = dphi * 180 / (CLHEP::pi);
    dang = dang * 180 / (CLHEP::pi);
    if (fabs(dang) < m_gammaAngleCut)
      continue;
    if (getTime > m_gammathCut || getTime < m_gammatlCut)
      continue;
    if (!((fabs(cosThetaSh) < m_maxCosThetaBarrel && eraw > m_minEnergy) ||
          ((fabs(cosThetaSh) > m_minCosThetaEndcap) && (fabs(cosThetaSh) < m_maxCosThetaEndcap) &&
           (eraw > m_minEndcapEnergy))))
      continue;
    if (m_debug)
      cout << __LINE__ << " i " << i << " shP4.m() " << shP4.e() << endl;

    MyParticle tmp(i, shP4, emcTrk);
    emcGamma.push_back(tmp);
  }
  // Finish Good Photon Slection
#pragma endregion

  if (emcGamma.size() < 4)
  {
    cout << __LINE__ << "return StatusCode::SUCCESS;  emcGamma.size() < 4" << endl;
    return StatusCode::SUCCESS;
  }
  
  if (abs(signal) == 1)
  {
    if (m_debug)
      cout << __LINE__ << " " << "emcGamma.size(): " << emcGamma.size() << endl;
    Ncut2++;
  }

  if (m_debug)
    cout << __LINE__ << endl;

#pragma region loop_gamma_check_eta_pi0______________________________________________________________

  std::vector<MyMotherParticleFit> eta, eta_1c;
  eta.clear();
  eta_1c.clear();

  KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
  // Loop each gamma pair, check eta mass  ----------------------------
  for (int i = 0; i < emcGamma.size() - 1; i++)
  {
    RecEmcShower *emcTrki = emcGamma[i].getRecEmcShower();
    HepLorentzVector ptrki = emcGamma[i].getLorentzVector();

    for (int j = i + 1; j < emcGamma.size(); j++)
    {
      RecEmcShower *emcTrkj = emcGamma[j].getRecEmcShower();
      HepLorentzVector ptrkj = emcGamma[j].getLorentzVector();

      HepLorentzVector p2geta = ptrki + ptrkj;
      if (m_debug)
        cout << __LINE__ << " i,j  " << i << "," << j << " p2geta.m()  " << p2geta.m() << endl;

      if (p2geta.m() < m_EtaMinMass || p2geta.m() > m_EtaMaxMass)
        continue;
      if (m_debug)
        cout << __LINE__ << " 00000000 " << " i,j  " << i << "," << j << " p2geta.m()  " << p2geta.m() << endl;

        MyMotherParticleFit tmp(emcGamma[i], emcGamma[j]);
        // tmp.setLorentzVector(p2geta);
        eta.push_back(tmp);
      }
    }
  }
  // Loop each gamma pair, check pi0 mass ---------------------
  std::vector<MyMotherParticleFit> pi0, pi0_1c;
  pi0.clear();
  pi0_1c.clear();
  for (int k = 0; k < emcGamma.size() - 1; k++)
  {
    RecEmcShower *emcTrkk = emcGamma[k].getRecEmcShower();
    HepLorentzVector ptrkk = emcGamma[k].getLorentzVector();
    for (int l = k + 1; l < emcGamma.size(); l++)
    {
      RecEmcShower *emcTrkl = emcGamma[l].getRecEmcShower();
      HepLorentzVector ptrkl = emcGamma[l].getLorentzVector();

      HepLorentzVector p2gpi = ptrkk + ptrkl;
      if (m_debug)
        cout << __LINE__ << " k,l " << k << "," << l << " p2gpi.m() " << p2gpi.m() << endl;

      if (p2gpi.m() < m_Pi0MinMass || p2gpi.m() > m_Pi0MaxMass)
        continue;
      if (m_debug)
        cout << __LINE__ << " 00000000 "<< " k,l " << k << "," << l << " p2gpi.m() " << p2gpi.m() << endl;

        MyMotherParticleFit tmp(emcGamma[k], emcGamma[l]);
        // tmp.setLorentzVector(p2gpi);
        pi0.push_back(tmp);

      }
    }
  }
  if(m_debug)
  {
    cout << __LINE__ <<" eta.size() " << eta.size()  << endl;
    cout << __LINE__ <<" pi0.size() " << pi0.size()  << endl;
  }

  // if (eta_1c.size() != eta.size() || pi0_1c.size() != pi0.size() )
  // {
  //   if(m_debug)
  //     cout << __LINE__ << "return StatusCode::SUCCESS; eta_1c.size() != eta.size() || pi0_1c.size() != pi0.size()" << endl;
  //   return StatusCode::SUCCESS;
  // }

  if (abs(signal) == 1)
    Ncut3++; // Ncut3 should equal Ncut2;


  // if (pi0.size() == 0 || eta.size() == 0 || eta_1c.size() == 0|| pi0_1c.size()==0)
  // {
  //   if(m_debug)    
  //     cout << __LINE__ << "return StatusCode::SUCCESS; pi0.size() == 0 || eta.size() == 0 || eta_1c.size() == 0|| pi0_1c.size()==0" << endl;
  //   return StatusCode::SUCCESS;
  // }
  if (abs(signal) == 1)
    Ncut4++; // Ncut4 should equal Ncut5;

  
  // Loop each pi+ pi- gamma pair, check eta' mass ---------------------
  std::vector<MyMotherParticleFit> pipm_afterSelectEtaPrime, eta_afterSelectEtaPrime;
  pipm_afterSelectEtaPrime.clear();
  eta_afterSelectEtaPrime.clear();

  for(int k = 0; k < eta.size(); k++)
  {
    for (int l = 0; l < pi0.size(); l++)
    {        
      for (int i = 0; i < piMin.size(); i++)
      {
        for(int j = 0; j< piPlus.size(); j++)
        {
          // cout << __LINE__ << endl;
          // cout << eta_1c[k].getChild1().getIndex() << " " <<  pi0_1c[l].getChild1().getIndex() << "    "
          //      << eta_1c[k].getChild1().getIndex() << " " << pi0_1c[l].getChild2().getIndex()<< endl;

          if (eta[k].getChild1().getIndex() == pi0[l].getChild1().getIndex() ||
              eta[k].getChild1().getIndex() == pi0[l].getChild2().getIndex())
            continue;
          // cout << __LINE__ << endl;
          // cout << eta_1c[k].getChild2().getIndex() << " " <<  pi0_1c[l].getChild1().getIndex() << "    "
          //      << eta_1c[k].getChild2().getIndex() << " " << pi0_1c[l].getChild2().getIndex()<< endl;

          if (eta[k].getChild2().getIndex() == pi0[l].getChild1().getIndex() ||
              eta[k].getChild2().getIndex() == pi0[l].getChild2().getIndex())
          continue;

          HepLorentzVector tmp = piMin[i].getLorentzVector() + piPlus[j].getLorentzVector() + eta_1c[k].getMotherLorentzVector(2);
          
          if(m_debug)
            cout << __LINE__ << " eta prime m(): " << tmp.m() << endl;

          if (tmp.m() > m_EtaPrimeMaxMass || tmp.m() < m_EtaPrimeMinMass )
            continue;
            
          if(m_debug)
            cout << __LINE__ <<  "00000000" << " eta prime m(): " << tmp.m() << endl;

          // MyMotherParticleFit tmp(piMin[i], piPlus[j], eta[k]);
          pipm_afterSelectEtaPrime.push_back( MyMotherParticleFit(piMin[i], piPlus[j]));
          eta_afterSelectEtaPrime.push_back(eta[k]);
        }
      }
    }
  }

  if(pipm_afterSelectEtaPrime.size()!=eta_afterSelectEtaPrime.size() )
  {
    if(m_debug)
      cout << __LINE__ << "return StatusCode::SUCCESS; pipm_afterSelectEtaPrime.size()!=eta_afterSelectEtaPrime.size() " << endl;
    return StatusCode::SUCCESS;
  }

  if (abs(signal) == 1)
    Ncut5++; // Ncut4 should equal Ncut5;

  if( pipm_afterSelectEtaPrime.size()==0 || eta_afterSelectEtaPrime.size()==0)
  {
    if(m_debug)
      cout << __LINE__ << "return StatusCode::SUCCESS; pipm_afterSelectEtaPrime.size()==0 || eta_afterSelectEtaPrime.size()==0" << endl;
    return StatusCode::SUCCESS;
  }

#pragma endregion


  if (m_debug)
    cout << __LINE__ << endl;

  // get beam energy and beta
  if (m_ReadBeamEFromDB)
  {
    if (m_usecalibBeamE)
      m_readDb.setcalib(true);
    m_beamE = m_readDb.getbeamE(m_run, m_beamE);
    if (m_run > 0)
      m_beta = m_readDb.getbeta();
    if (m_debug)
      cout << __LINE__ << "beam from db:" << m_beamE << ", mbeta: " << m_beta << endl;
  }

  double ebeam = m_beamE;
  HepLorentzVector HepCMS(0.011 * 2 * ebeam, 0., 0., 2 * ebeam);

#pragma region save_lambda_c - _and_lambda_c + ________________________________________________________________________

  // for data from raw
  HepLorentzVector eta_pg1(0, 0, 0, 0), eta_pg2(0, 0, 0, 0), pi0_pg3(0, 0, 0, 0), pi0_pg4(0, 0, 0, 0), p_p4(0, 0, 0, 0), pip_p4(0, 0, 0, 0), pim_p4(0, 0, 0, 0);
  // // for best chi2
  HepLorentzVector p_p4_r1c(0, 0, 0, 0), pi0g1_p4_r1c(0, 0, 0, 0), etag1_p4_r1c(0, 0, 0, 0), pi0g2_p4_r1c(0, 0, 0, 0),
      etag2_p4_r1c(0, 0, 0, 0), pip_p4_r1c(0,0,0,0), pim_p4_r1c(0,0,0,0);
  // for best delta E
  // HepLorentzVector pi0g1_p4_1c(0, 0, 0, 0), etag1_p4_1c(0, 0, 0, 0),
  //     pi0g2_p4_1c(0, 0, 0, 0), etag2_p4_1c(0, 0, 0, 0);

  double minChi2_r1c = 999999999, minChi2_r2c = 99999999999, deltaE_minb = 99999999;
  double chisq = 0;
  int flag = 0, pcharge=0;

  for (int i = 0; i < proton.size(); i++)
  {
    for (int j = 0; j < pi0.size(); j++)
    {
      for (int k = 0; k < eta_afterSelectEtaPrime.size(); k++)
      {
        HepLorentzVector psigma = proton[i].getLorentzVector() + pi0[j].getMotherLorentzVector(2);
        if (m_debug)
          cout << __LINE__ << " i:"<< i << "  j:" << j << " k :"  <<  k << "  psigma.m():" << psigma.m() << endl;
          
        if (psigma.m() < m_SigmaMinMass || psigma.m() > m_SigmaMaxMass)
          continue;

        if (m_debug)
          cout << __LINE__ << " i:"<< i << "  j:" << j << " k :"  <<  k << "  psigma.m():" << psigma.m() << endl;

        // HepLorentzVector pLambda = proton[i].getLorentzVector() + pi0[j].getMotherLorentzVector(2) + 
        //                            eta_afterSelectEtaPrime[k].getMotherLorentzVector(2) + 
        //                            pipm_afterSelectEtaPrime[k].getMotherLorentzVector(2);
        // pLambda.boost(-m_beta);
        // double deltaEb = pLambda.t() - ebeam;

        // //   // ______________________________  minimum delta E ____________________________________
        // if (m_debug)
        //   cout << "fabs(deltaEb): " << fabs(deltaEb) << ", fabs(deltaE_minb): " << fabs(deltaE_minb) << endl;
        // if (fabs(deltaEb) < fabs(deltaE_minb))
        // {
        //   flag = 1;
        //   deltaE_minb = deltaEb;
        //   pcharge = proton[i].getCharge();
        //   // chisq = chi2[k];
        //   // m_p_p4_r1c = kmfit1->pfit(0);
        //   pi0g1_p4_1c = pi0_1c[j].getChild1().getLorentzVector();
        //   pi0g2_p4_1c = pi0_1c[j].getChild2().getLorentzVector();
        //   etag1_p4_1c = eta_1c_afterSelectEtaPrime[k].getChild1().getLorentzVector();
        //   etag2_p4_1c = eta_1c_afterSelectEtaPrime[k].getChild2().getLorentzVector();

        //   p_p4 = proton[i].getLorentzVector();

        //   eta_pg1 = eta_afterSelectEtaPrime[k].getChild1().getLorentzVector();
        //   eta_pg2 = eta_afterSelectEtaPrime[k].getChild2().getLorentzVector();
        //   pi0_pg3 = pi0[j].getChild1().getLorentzVector();
        //   pi0_pg4 = pi0[j].getChild2().getLorentzVector();

        //   pim_p4 = pipm[k].getChild1().getLorentzVector();
        //   pip_p4 = pipm[k].getChild2().getLorentzVector();
        // }
        
        // _______________________________________________  1C  ______________________________________________
        KalmanKinematicFit *kmfit2 = KalmanKinematicFit::instance();
        kmfit2->init();
        kmfit2->setChisqCut(1e3);
        kmfit2->setIterNumber(10);

        // kmfit2->AddTrack(0, proton[i].getTrackParameter());
        // kmfit2->AddTrack(1, pi0[k].getFit()->wVirtualTrack(0));
        // kmfit2->AddTrack(2, eta[j].getFit()->wVirtualTrack(0));

        kmfit2->AddTrack(0, proton[i].getTrackParameter());
        kmfit2->AddTrack(1, 0.0, pi0[j].getChild1().getRecEmcShower());
        kmfit2->AddTrack(2, 0.0, pi0[j].getChild2().getRecEmcShower());
        kmfit2->AddTrack(3, 0.0, eta_afterSelectEtaPrime[k].getChild1().getRecEmcShower());
        kmfit2->AddTrack(4, 0.0, eta_afterSelectEtaPrime[k].getChild2().getRecEmcShower());
        kmfit2->AddTrack(5, 0.0, pipm_afterSelectEtaPrime[k].getChild1().getTrackParameter());
        kmfit2->AddTrack(6, 0.0, pipm_afterSelectEtaPrime[k].getChild2().getTrackParameter());

        kmfit2->AddMissTrack(7, 2.28646);

        kmfit2->AddResonance(0, 0.547862, 3, 4);
        kmfit2->AddResonance(1, 0.1349770, 1, 2);
        kmfit2->AddFourMomentum(2, HepCMS);

        // MyMotherParticleFit tmp2;
        bool okvs1 = kmfit2->Fit();
        cout << __LINE__  << endl;
        if (okvs1)
        {
          kmfit2->BuildVirtualParticle(0);
          // LcWTrk_1C = kmfit2->wVirtualTrack(0);
          
          // ______________________________  1c minimum chi2 ____________________________________
          if (m_debug) 
            cout << __LINE__ << " minChi2: " << minChi2_r1c << " chi2:" << kmfit2->chisq() << endl;
          if (kmfit2->chisq() < minChi2_r1c)
          {
            flag = 1;
            pcharge = proton[i].getCharge();

            minChi2_r1c = kmfit2->chisq();
            
            p_p4_r1c = kmfit2->pfit(0);
            pi0g1_p4_r1c = kmfit2->pfit(1);
            pi0g2_p4_r1c = kmfit2->pfit(2);
            etag1_p4_r1c = kmfit2->pfit(3);
            etag2_p4_r1c = kmfit2->pfit(4);
            pim_p4_r1c = kmfit2->pfit(5);
            pip_p4_r1c = kmfit2->pfit(6);

            p_p4 = proton[i].getLorentzVector();
            pim_p4 = pipm_afterSelectEtaPrime[k].getChild1().getLorentzVector();
            pip_p4 = pipm_afterSelectEtaPrime[k].getChild2().getLorentzVector();
            eta_pg1 = eta[k].getChild1().getLorentzVector();
            eta_pg2 = eta[k].getChild2().getLorentzVector();
            pi0_pg3 = pi0[j].getChild1().getLorentzVector();
            pi0_pg4 = pi0[j].getChild2().getLorentzVector();

          }
        }
        if (m_debug &&  flag==0 )
          cout << __LINE__ <<  " fit wrong"  << endl;

      }
    }
  }
#pragma endregion

  if (m_debug)
    cout << __LINE__ << endl;

#pragma region write__________________________________________________________________

  if (flag == 1)
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

    // .....
    m_signal = signal;
    // m_bg = bg;
    m_ebeam = ebeam;
    
    //  _____________  1c  _______________________
    for (int jj = 0; jj < 4; jj++)
      m_pall_p4_r1c[jj] = p_p4_r1c[jj];
    for (int jj = 0; jj < 4; jj++)
      m_pim_p4_r1c[jj] = pim_p4_r1c[jj];
    for (int jj = 0; jj < 4; jj++)
      m_pip_p4_r1c[jj] = pip_p4_r1c[jj];
    //   1,2 -> eta             3,4 -> pi       
    for (int jj = 0; jj < 4; jj++)
      m_gam1_p4_1c[jj] = etag1_p4_r1c[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam2_p4_1c[jj] = etag2_p4_r1c[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam3_p4_1c[jj] = pi0g1_p4_r1c[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam4_p4_1c[jj] = pi0g2_p4_r1c[jj];

    // ___________ from raw ___________________
    // proton, pi+, pi-, from raw
    for (int jj = 0; jj < 4; jj++)
      m_pall_p4[jj] = p_p4[jj];
    for (int jj = 0; jj < 4; jj++)
      m_pim_p4[jj] = pim_p4[jj];
    for (int jj = 0; jj < 4; jj++)
      m_pip_p4[jj] = pip_p4[jj];
    // for four gammas, from raw
    for (int jj = 0; jj < 4; jj++)
      m_gam1_p4[jj] = eta_pg1[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam2_p4[jj] = eta_pg2[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam3_p4[jj] = pi0_pg3[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam4_p4[jj] = pi0_pg4[jj];

    // m_pbarindex = pbar_index;
    m_p4index = 4;
    m_pcharge = pcharge;
    m_chi2_min_r1c = minChi2_r1c;

    m_pi0m = (pi0_pg3 + pi0_pg4).m();
    m_etam = (eta_pg1 + eta_pg2).m();
    m_Sigmam = (p_p4 + pi0_pg3 + pi0_pg4).m();
    m_etaprimem = (pim_p4 + pip_p4 + eta_pg1 + eta_pg2).m();

    if(m_debug)
      cout << __LINE__ << " m_pi0m " << m_pi0m << " m_etam " << m_etam << " m_Sigmam " << m_Sigmam
           << " m_etaprimem " << m_etaprimem << endl;

    m_pi0mR1c = (pi0g1_p4_r1c + pi0g2_p4_r1c).m();
    m_etamR1c = (etag1_p4_r1c + etag2_p4_r1c).m();
    m_SigmamR1c = (p_p4_r1c + pi0g1_p4_r1c + pi0g2_p4_r1c).m();
    m_etaprimemR1c = (pim_p4_r1c + pip_p4_r1c + etag1_p4_r1c + etag2_p4_r1c).m();
    
    if(m_debug)
      cout << __LINE__ << " m_pi0m1c " << m_pi0m1c << " m_etam1c " << m_etam1c << " m_Sigmam1c " <<  m_Sigmam1c 
           << " m_etaprimem1c " <<  m_etaprimem1c << endl;

    HepLorentzVector pLambda = p_p4_r1c + pim_p4_r1c + pip_p4_r1c + etag1_p4_r1c + etag2_p4_r1c + pi0g1_p4_r1c + pi0g2_p4_r1c;
    if(m_debug)
      cout << __LINE__ << " pLambda.m() " << pLambda.m() << endl;
    pLambda.boost(-m_beta);
    m_deltaE_min = pLambda.t() - ebeam;

    double mbc2 = ebeam * ebeam - pLambda.v().mag2();
    m_bc = mbc2 > 0 ? sqrt(mbc2) : -10;

    cout << __LINE__  << " m_bc " << m_bc  << " m_deltaE_min " << m_deltaE_min<< endl;

    // m_rightflag = 2;
    m_np = np;
    m_npbar = npbar;

    m_tuple1->write();
    if (m_debug)
      cout << __LINE__ << " write() " << endl;
    
    Ncut7++;
  }

#pragma endregion

  Ncut6++;
  cout << endl;
  cout << "attention: if Ncut2!= Ncut3 , Ncut4 != Ncut5;  you should check" << endl;
  cout << "Ntotal  " << Ntotal << endl;
  cout << "Ncut0   " << Ncut0 << endl;
  cout << "Ncut1   " << Ncut1 << endl;
  cout << "Ncut2   " << Ncut2 << endl;
  cout << "Ncut3   " << Ncut3 << endl;
  cout << "Ncut4   " << Ncut4 << endl;
  cout << "Ncut5   " << Ncut5 << endl;
  cout << "Ncut6   " << Ncut6 << endl;
  cout << "Ncut7   " << Ncut7 << endl;
  cout << "all       " << all << endl;

  return StatusCode::SUCCESS;
}
StatusCode LambdacAlg::endRun()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "LambdacAlg::endRun()" << endreq;
  // add your code here
  return StatusCode::SUCCESS;
}
StatusCode LambdacAlg::finalize()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "LambdacAlg::finalize()" << endreq;
  // add your code here
  cout << endl;
  cout << "attention: if Ncut2!= Ncut3 , Ncut4 != Ncut5;  you should check" << endl;
  cout << "Ntotal  " << Ntotal << endl;
  cout << "Ncut0   " << Ncut0 << endl;
  cout << "Ncut1   " << Ncut1 << endl;
  cout << "Ncut2   " << Ncut2 << endl;
  cout << "Ncut3   " << Ncut3 << endl;
  cout << "Ncut4   " << Ncut4 << endl;
  cout << "Ncut5   " << Ncut5 << endl;
  cout << "Ncut6   " << Ncut6 << endl; 
  cout << "all       " << all << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  cout << "-------------------------------           -------------------------------" << endl;
  cout << "--------------------                                ---------------------" << endl;
  cout << "-----------------  s i g m a    e t a '  1 c   v1.1 ------------------" << endl;
  cout << "--------------------                               ----------------------" << endl;
  cout << "------------------------------           --------------------------------" << endl;
  cout << "-------------------------------------------------------------------------" << endl;

  //	cout << "al       " << al << endl;
  //	cout << "yes     " << yes << endl;
  //	cout << "Npp>0:  " << Npp << endl;
  //	cout << "Npm>0:  " << Npm << endl;

  return StatusCode::SUCCESS;
}

// add your code here,for other member-functions
