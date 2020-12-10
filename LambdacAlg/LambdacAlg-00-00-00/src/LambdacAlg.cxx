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
const double xmass[5] = {0.000511, 0.105658, 0.139570, 0.493677, 0.938272};
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

      // 1c _________________________________________________________________
      status = m_tuple1->addIndexedItem("pall", m_p4index, m_pall_p4);
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
      status = m_tuple1->addItem("pi0mORetam1mR1c", m_pi0mORetam1mR1c);
      status = m_tuple1->addItem("sigmamR1c", m_SigmamR1c);
      // 2c _________________________________________________________________



      // ____________________________________________________________________
      status = m_tuple1->addItem("pi0m", m_pi0m);
      status = m_tuple1->addItem("pi0mORetam1m", m_pi0mORetam1m);

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
      std::cout << "Could not retrieve McParticelCol" << std::endl;
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
          cout << "mcParticleCol pdg: " << pdg << ", motherpdg: " << motherpdg << ", mmotherpdg:" << mmotherpdg << endl;

        // truth for lambda_c+ to ... ------------------------------------
        if ((*iter_mc)->particleProperty() == 4122)
        {
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << "4122 ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;

          } // End of "gc.size() > 0" IF
        }

        // mode2 mm_mode2 == 28 &&
        if (pdg == 3222 && motherpdg == 4122)
        { // Sgm+ ->p pi0
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << "4122 3222  ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

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
              cout << "4122 3222 111 ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          } // End of "gc.size() > 0" IF
        }
        if (pdg == 111 && motherpdg == 4122)
        { //// eta->gam gam
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << "4122 111 ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

            for (int ll = 0; ll < 4; ll++)
              Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAp++;
          }
        }

        // truth for lambda_c- to ... ------------------------------------
        if ((*iter_mc)->particleProperty() == -4122)
        {
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();

            if (m_debug)
              cout << "-4122 ndaughterAm[" << ndaughterAm << "]: " << Am_id[ndaughterAm] << endl;

            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }
        // mode3
        if (pdg == -3222 && motherpdg == -4122)
        { // Sgm- ->pbar pi0
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();

            if (m_debug)
              cout << "-4122 111 ndaughterAm[" << ndaughterAm << "]: " << Ap_id[ndaughterAm] << endl;

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
              cout << "4122 111 ndaughterAm[" << ndaughterAm << "]: " << Ap_id[ndaughterAm] << endl;

            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }
        if (pdg == 111 && motherpdg == -4122)
        { //// eta->gam gam
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();

            if (m_debug)
              cout << "4122 111 ndaughterAm[" << ndaughterAm << "]: " << Ap_id[ndaughterAm] << endl;

            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          }
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
              cout << "M_pdgid[" << i << "]=" << M_pdgid[i] << endl;
              cout << "M_motheridx[" << i << "]=" << M_motheridx[i] << endl;
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
            cout << __LINE__ << "+iter_mc: " << (*iter_mc)->particleProperty() << ", numParticle: " << numParticle
                 << endl;
          }

          for (int i = 0; i != pdgid.size(); i++)
          {
            M_pdgid_p[i] = pdgid[i];
            M_motheridx_p[i] = motherindex[i];

            if (m_debug)
            {
              cout << "M_pdgid[" << i << "]=" << M_pdgid[i] << endl;
              cout << "M_motheridx[" << i << "]=" << M_motheridx[i] << endl;
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
            cout << "-iter_mc: " << (*iter_mc)->particleProperty() << ", numParticle: " << numParticle << endl;
          }

          for (int i = 0; i != pdgid.size(); i++)
          {
            M_pdgid_m[i] = pdgid[i];
            M_motheridx_m[i] = motherindex[i];

            if (m_debug)
            {
              cout << "M_pdgid[" << i << "]=" << M_pdgid[i] << endl;
              cout << "M_motheridx[" << i << "]=" << M_motheridx[i] << endl;
            }
          }
        }
      }

      if (ndaughterAp == 8 && Ap_id[0] == 3222 && Ap_id[1] == 111 && Ap_id[2] == 2212 && Ap_id[3] == 111 &&
          Ap_id[4] == 22 && Ap_id[5] == 22 && Ap_id[6] == 22 && Ap_id[7] == 22)
      {
        signal = 1;
        all++;
        all_p++;
      }
      if (ndaughterAm == 8 && Am_id[0] == -3222 && Am_id[1] == 111 && Am_id[2] == -2212 && Am_id[3] == 111 &&
          Am_id[4] == 22 && Am_id[5] == 22 && Am_id[6] == 22 && Am_id[7] == 22)
      {
        signal = -1;
        all++;
        all_m++;
      }

      if (m_debug)
      {
        cout << __LINE__ << " all: " << all << ", +signal: " << all_p << ", -signal: " << all_m << endl;
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
  if (goodTrack.size() < 1)
    return StatusCode::SUCCESS;

  Vint trackProntonP, trackProntonPbar;
  std::vector<MyParticle> proton;
  proton.clear();
  trackProntonP.clear();
  trackProntonPbar.clear();
  for (int i = 0; i < goodTrack.size(); i++)
  {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + goodTrack[i];
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
    MyPid *pid = new MyPid(*itTrk);
    if (pid->isproton())
    {
      if (mdcTrk->charge() == 1)
        trackProntonP.push_back(goodTrack[i]);
      if (mdcTrk->charge() == -1)
        trackProntonPbar.push_back(goodTrack[i]);

      RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
      mdcKalTrk->setPidType(RecMdcKalTrack::proton);
      HepLorentzVector p4 = mdcKalTrk->p4(xmass[4]);

      // EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ip[i];
      // RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
      // mdcKalTrk->setPidType(RecMdcKalTrack::proton);
      // HepLorentzVector pp4 = mdcKalTrk->p4(xmass[4]);

      WTrackParameter wtrkp(xmass[4], mdcKalTrk->getZHelixP(), mdcKalTrk->getZErrorP());
      if (m_debug)
        cout << __LINE__ << " i " << i << " p4.m() " << p4.m() << endl;
      MyParticle tmp(goodTrack[i], p4, mdcTrk->charge(), wtrkp);
      tmp.chi = pid->getChi();
      proton.push_back(tmp);
    }
  }

  if (abs(signal) == 1)
    Ncut0++;

  int np = trackProntonP.size();
  int npbar = trackProntonPbar.size();
  // m_np = np;
  // m_npbar = npbar;

  // if (m_isqqbar)
  // {
    // if (np == 0 && npbar == 0)
    //   return StatusCode::SUCCESS;
  // }

  if (abs(signal) == 1)
    Ncut1++;

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

  if (m_debug)
    cout << __LINE__ << endl;
  if (emcGamma.size() < 4)
    return StatusCode::SUCCESS;
  if (abs(signal) == 1)
  {
    if (m_debug)
      cout << __LINE__ << " " << "emcGamma.size(): " << emcGamma.size() << endl;
    Ncut2++;
  }

  if (m_debug)
    cout << __LINE__ << endl;

#pragma region loop_gamma_check_2_pi0______________________________________________________________

  // Vint igam1, igam2, igam3, igam4;
  // igam1.clear(), igam2.clear(), igam3.clear(), igam4.clear();

  // Vp4 pgam1, pgam2, pgam3, pgam4, pgam1_1C, pgam2_1C, pgam3_1C, pgam4_1C;
  // pgam1.clear(), pgam2.clear(), pgam3.clear(), pgam4.clear(), pgam1_1C.clear(), pgam2_1C.clear(), pgam3_1C.clear(),
  //     pgam4_1C.clear();

  // Vdouble chi2;
  // chi2.clear();

  // std::vector<MyMotherParticleFit> eta;
  // eta.clear();

  // // KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
  // // Loop each gamma pair, check eta mass  ----------------------------
  // for (int i = 0; i < emcGamma.size() - 1; i++)
  // {
  //   RecEmcShower *emcTrki = emcGamma[i].getRecEmcShower();
  //   HepLorentzVector ptrki = emcGamma[i].getLorentzVector();

  //   for (int j = i + 1; j < emcGamma.size(); j++)
  //   {
  //     RecEmcShower *emcTrkj = emcGamma[j].getRecEmcShower();
  //     HepLorentzVector ptrkj = emcGamma[j].getLorentzVector();

  //     HepLorentzVector p2geta = ptrki + ptrkj;
  //     // if (m_debug)
  //     //   cout << __LINE__ << " i,j  " << i << "," << j << " p2geta.m()  " << p2geta.m() << endl;

  //     if (p2geta.m() < m_EtaMinMass || p2geta.m() > m_EtaMaxMass)
  //       continue;
  //     if (m_debug)
  //       cout << __LINE__ << " 00000000 " << " p2geta.m()  " << p2geta.m() << endl;
  //     // if (m_test1C == 1)
  //     // {
  //     //   kmfit->init();
  //     //   kmfit->setChisqCut(m_chisqMax);

  //     //   kmfit->AddTrack(0, 0.0, emcTrki);
  //     //   kmfit->AddTrack(1, 0.0, emcTrkj);
  //     //   // 0.547862
  //     //   kmfit->AddResonance(0, 0.547862, 0, 1);
  //     //   bool oksq = kmfit->Fit();
  //     //   if (oksq)
  //     //   {
  //     //     kmfit->BuildVirtualParticle(0);
  //     MyMotherParticleFit tmp(emcGamma[i], emcGamma[j]);
  //     tmp.setLorentzVector(p2geta);
  //     eta.push_back(tmp);
  //       // }
  //     // }
  //   }
  // }

  // Loop each gamma pair, check pi0 mass ---------------------
  std::vector<MyMotherParticleFit> pi0;
  pi0.clear();
  // KalmanKinematicFit *kmfit_pi = KalmanKinematicFit::instance();
  for (int k = 0; k < emcGamma.size() - 1; k++)
  {
    RecEmcShower *emcTrkk = emcGamma[k].getRecEmcShower();
    HepLorentzVector ptrkk = emcGamma[k].getLorentzVector();
    for (int l = k + 1; l < emcGamma.size(); l++)
    {
      RecEmcShower *emcTrkl = emcGamma[l].getRecEmcShower();
      HepLorentzVector ptrkl = emcGamma[l].getLorentzVector();

      HepLorentzVector p2gpi = ptrkk + ptrkl;
      // if (m_debug)
      //   cout << __LINE__ << " k,l " << k << "," << l << " p2gpi.m() " << p2gpi.m() << endl;

      if (p2gpi.m() < m_Pi0MinMass || p2gpi.m() > m_Pi0MaxMass)
        continue;
      if (m_debug)
        cout << __LINE__ << " 00000000 " << " p2gpi.m() " << p2gpi.m() << endl;

      // if (m_test1C == 1)
      // {
      //   kmfit_pi->init();
      //   kmfit_pi->setChisqCut(m_chisqMax);

      //   kmfit_pi->AddTrack(0, 0.0, emcTrkk);
      //   kmfit_pi->AddTrack(1, 0.0, emcTrkl);
      //   // 0.1349770
      //   kmfit_pi->AddResonance(0, 0.1349770, 0, 1);
      //   bool oksq = kmfit_pi->Fit();

      //   if (oksq)
      //   {
      //     // double chisq_1C = kmfit_pi->chisq(0);
      //     kmfit_pi->BuildVirtualParticle(0);
      //     // WTrackParameter Pi0WTrk = kmfit_pi->wVirtualTrack(0);

      MyMotherParticleFit tmp(emcGamma[k], emcGamma[l]);
      tmp.setLorentzVector(p2gpi);
      pi0.push_back(tmp);
        // }
      // }
    }
  }
  cout << __LINE__ << " proton.size() "<< proton.size() << " pi0.size() " << pi0.size() /*<< " eta.size() "<< eta.size()*/ << endl;
  // int ngam1 = igam1.size();
  // int ngam2 = igam2.size();
  // int ngam3 = igam3.size();
  // int ngam4 = igam4.size();
  // if ((ngam3 != ngam4) || (ngam1 != ngam2))
  //   return StatusCode::SUCCESS;
  // int ngam12 = ngam1;
  // int ngam34 = ngam3;

  // if (m_debug)
  //   cout << __LINE__ << " "
  //        << "2*ngam12: " << 2 * ngam12 << ", 2*ngam34: " << 2 * ngam34 << endl;
#pragma endregion

  if (abs(signal) == 1)
    Ncut3++; // Ncut3 should equal Ncut2;

  if (pi0.size() < 2 )
  {
    cout << __LINE__ << "return StatusCode::SUCCESS; pi0.size() < 2" << endl;
    return StatusCode::SUCCESS;
  }
  if (m_debug)
    cout << __LINE__ << endl;

  if (abs(signal) == 1)
    Ncut4++;

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

  HepLorentzVector eta_pg1(0, 0, 0, 0), eta_pg2(0, 0, 0, 0), pi_pg3(0, 0, 0, 0), pi_pg4(0, 0, 0, 0), p_p4(0,0,0,0),
                   m_p_p4_r1c(0, 0, 0, 0), m_pi0g1_p4_r1c(0, 0, 0, 0), m_etag1_p4_r1c(0, 0, 0, 0),
                   m_pi0g2_p4_r1c(0, 0, 0, 0), m_etag2_p4_r1c(0, 0, 0, 0);
  double minChi2_r1c = 999999999, minChi2_r2c = 99999999999;
  int flag = 0;

  for (int i = 0; i < proton.size(); i++)
  {
    for (int j = 0; j < pi0.size(); j++)
    {
      // k is the pi from sigma
      for (int k = 0; k < pi0.size(); k++)
      {
        // if (pi0[j].getChild1().getIndex() == pi0[k].getChild1().getIndex() ||
        //     pi0[j].getChild1().getIndex() == pi0[k].getChild2().getIndex())
        //   continue;
        // if (pi0[j].getChild2().getIndex() == pi0[k].getChild1().getIndex() ||
        //     pi0[j].getChild2().getIndex() == pi0[k].getChild2().getIndex())
        //   continue;

        // if (m_debug)
        //   cout << __LINE__ << proton[i].getLorentzVector()[0] << " " << proton[i].getLorentzVector()[1] << " "
        //        << proton[i].getLorentzVector()[2] << " " << proton[i].getLorentzVector()[3] << " " << endl;
        if(j==k)
            continue;

        // MyMotherParticleFit tmp0(proton[i], pi0[k]);
        // tmp0.setLorentzVector(psigma);
        // sigma = tmp0;

        // MyMotherParticleFit tmp(proton[i], pi0[k], eta[j]);

        // HepLorentzVector pLambda = proton[i].getLorentzVector() + pi0[k].getLorentzVector() + eta[j].getLorentzVector();
        // pLambda.boost(-m_beta);
        // tmp.setLorentzVector(pLambda);
        // lambdac = tmp;
        
        // 1C
        KalmanKinematicFit *kmfit1 = KalmanKinematicFit::instance();
        kmfit1->init();
        kmfit1->setChisqCut(1e3);
        kmfit1->setIterNumber(10);

        // kmfit1->AddTrack(0, proton[i].getTrackParameter());
        // kmfit1->AddTrack(1, pi0[k].getFit()->wVirtualTrack(0));
        // kmfit1->AddTrack(2, eta[j].getFit()->wVirtualTrack(0));

        kmfit1->AddTrack(0, proton[i].getTrackParameter());
        kmfit1->AddTrack(1, 0.0, pi0[k].getChild1().getRecEmcShower());
        kmfit1->AddTrack(2, 0.0, pi0[k].getChild2().getRecEmcShower());
        kmfit1->AddTrack(3, 0.0, pi0[j].getChild1().getRecEmcShower());
        kmfit1->AddTrack(4, 0.0, pi0[j].getChild2().getRecEmcShower());
        kmfit1->AddMissTrack(5, 2.28646);

        kmfit1->AddResonance(0, 0.1349770, 3, 4);
        kmfit1->AddResonance(1, 0.1349770, 1, 2);
        kmfit1->AddFourMomentum(2, HepCMS);

        // MyMotherParticleFit tmp2;
        bool okvs1 = kmfit1->Fit();
        cout << __LINE__  << endl;
        if (okvs1)
        {
          kmfit1->BuildVirtualParticle(0);
          // LcWTrk_1C = kmfit1->wVirtualTrack(0);

        HepLorentzVector psigma = kmfit1->pfit(0) + kmfit1->pfit(1) + kmfit1->pfit(2);
        if (m_debug)
          cout << __LINE__ << " psigma.m() " << psigma.m() << endl;
        if (psigma.m() < m_SigmaMinMass || psigma.m() > m_SigmaMaxMass)
          continue;
        if (m_debug)
          cout << __LINE__ << "  psigma.m():" << psigma.m() << endl;

          
          // ______________________________  1c minimum chi2 ____________________________________
          if (m_debug) 
            cout << __LINE__ << " minChi2: " << minChi2_r1c << " chi2:" << kmfit1->chisq() << endl;
          if (kmfit1->chisq() < minChi2_r1c)
          {
            flag = 1;
            minChi2_r1c = kmfit1->chisq();
            
            m_p_p4_r1c = kmfit1->pfit(0);
            m_pi0g1_p4_r1c = kmfit1->pfit(1);
            m_pi0g2_p4_r1c = kmfit1->pfit(2);
            m_etag1_p4_r1c = kmfit1->pfit(3);
            m_etag2_p4_r1c = kmfit1->pfit(4);

            p_p4 = proton[i].getLorentzVector();
            eta_pg1 = pi0[j].getChild1().getLorentzVector();
            eta_pg2 = pi0[j].getChild2().getLorentzVector();
            pi_pg3 = pi0[k].getChild1().getLorentzVector();
            pi_pg4 = pi0[k].getChild2().getLorentzVector();

          }
        }
        if (m_debug &&  flag==0 )
          cout << __LINE__ <<  " fit wrong"  << endl;
        // // 2c
        // kmfit->init();
        // kmfit->setChisqCut(1e3);
        // kmfit->setIterNumber(10);

        // kmfit->AddTrack(0, proton[i].getTrackParameter());
        // kmfit->AddTrack(1, pi0[k].getFit()->wVirtualTrack(0));
        // kmfit->AddTrack(2, eta[j].getFit()->wVirtualTrack(0));
        // kmfit->AddMissTrack(3, 2.28646);

        // kmfit->AddResonance(0, 1.18937, 0, 1);
        // kmfit->AddFourMomentum(1, HepCMS);


        // bool okvs1 = kmfit->Fit();
        // if (okvs1)
        // {
        //   kmfit->BuildVirtualParticle(0);
        //   // LcWTrk_1C = kmfit->wVirtualTrack(0);
        //   MyMotherParticleFit tmp3(proton[i], pi0[k], eta[j], kmfit);
        //   recoilLambdac_r2c = tmp3;
        // }

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
    
    //   1,2 -> eta/pi             3,4 -> pi       
    for (int jj = 0; jj < 4; jj++)
       m_pall_p4_r1c[jj] = m_p_p4_r1c[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam1_p4_r1c[jj] = m_etag1_p4_r1c[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam2_p4_r1c[jj] = m_etag2_p4_r1c[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam3_p4_r1c[jj] = m_pi0g1_p4_r1c[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam4_p4_r1c[jj] = m_pi0g2_p4_r1c[jj];

    for (int jj = 0; jj < 4; jj++)
      m_pall_p4[jj] = p_p4[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam1_p4[jj] = eta_pg1[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam2_p4[jj] = eta_pg2[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam3_p4[jj] = pi_pg3[jj];
    for (int jj = 0; jj < 4; jj++)
      m_gam4_p4[jj] = pi_pg4[jj];

    // m_pbarindex = pbar_index;
    m_p4index = 4;
    m_chi2_min_r1c = minChi2_r1c;

    m_pi0m = (pi_pg3 + pi_pg4).m();
    m_pi0mORetam1m = (eta_pg1 + eta_pg2).m();
    m_Sigmam = (p_p4 + pi_pg3 + pi_pg4).m();

    cout << __LINE__ << " m_pi0m " << m_pi0m << " m_pi0mORetam1m " << m_pi0mORetam1m << " m_chi2_min_r1c " << m_chi2_min_r1c<< endl;


    m_pi0mR1c = (m_pi0g1_p4_r1c + m_pi0g2_p4_r1c).m();
    m_pi0mORetam1mR1c = (m_etag1_p4_r1c + m_etag2_p4_r1c).m();
    m_SigmamR1c = (m_p_p4_r1c + m_pi0g1_p4_r1c + m_pi0g2_p4_r1c).m();

    cout << __LINE__ << " m_pi0mR1c " << m_pi0mR1c << " m_pi0mORetam1mR1c " << m_pi0mORetam1mR1c << " m_SigmamR1c " <<  m_SigmamR1c << endl;


    HepLorentzVector pLambda = m_p_p4_r1c + m_pi0g1_p4_r1c + m_pi0g2_p4_r1c + m_etag1_p4_r1c + m_etag2_p4_r1c;
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
    
    Ncut6++;
  }

#pragma endregion

  Ncut5++;

  cout << "attention: if Ncut2!= Ncut3 ,you should check" << endl;
  cout << "Ntotal  " << Ntotal << endl;
  cout << "Ncut0   " << Ncut0 << endl;
  cout << "Ncut1   " << Ncut1 << endl;
  cout << "Ncut2   " << Ncut2 << endl;
  cout << "Ncut3   " << Ncut3 << endl;
  cout << "Ncut4   " << Ncut4 << endl;
  cout << "Ncut5   " << Ncut5 << endl;
  cout << "Ncut6   " << Ncut6 << endl;
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
  cout << "attention: if Ncut2!= Ncut3 ,you should check" << endl;
  cout << "Ntotal  " << Ntotal << endl;
  cout << "Ncut0   " << Ncut0 << endl;
  cout << "Ncut1   " << Ncut1 << endl;
  cout << "Ncut2   " << Ncut2 << endl;
  cout << "Ncut3   " << Ncut3 << endl;
  cout << "Ncut4   " << Ncut4 << endl;
  cout << "Ncut5   " << Ncut5 << endl;
  cout << "Ncut6   " << Ncut6 << endl; 
  cout << "all       " << all << endl;
  //	cout << "al       " << al << endl;
  //	cout << "yes     " << yes << endl;
  //	cout << "Npp>0:  " << Npp << endl;
  //	cout << "Npm>0:  " << Npm << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  cout << "-------------------------------           -------------------------------" << endl;
  cout << "--------------------                                ---------------------" << endl;
  cout << "--------------  for sigma pi0_r3c   v1.2 -----------------" << endl;
  cout << "--------------------                               ----------------------" << endl;
  cout << "------------------------------           --------------------------------" << endl;
  cout << "-------------------------------------------------------------------------" << endl;
  return StatusCode::SUCCESS;
}

// add your code here,for other member-functions
