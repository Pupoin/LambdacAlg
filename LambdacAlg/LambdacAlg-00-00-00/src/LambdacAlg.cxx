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

      // Proton
      status = m_tuple1->addItem("NumOfProton", m_nProton, 0, 200);
      status = m_tuple1->addIndexedItem("Proton_ID", m_nProton, m_Proton_ID);
      status = m_tuple1->addIndexedItem("Proton_Charge", m_nProton, m_Proton_Charge);
      status = m_tuple1->addIndexedItem("Proton_P4", m_nProton, 4, m_Proton_P4);
      // status = m_tuple1->addIndexedItem("Proton_PID", m_nProton, 5, m_Proton_PID);


      // Pi0
      status = m_tuple1->addItem("NumOfPi0", m_nPi0, 0, 200);
      status = m_tuple1->addIndexedItem("Pi0_Gam1_ID", m_nPi0, m_Pi0_Gam1_ID);
      status = m_tuple1->addIndexedItem("Pi0_Gam2_ID", m_nPi0, m_Pi0_Gam2_ID);
      status = m_tuple1->addIndexedItem("Pi0_m", m_nPi0, m_Pi0);
      status = m_tuple1->addIndexedItem("Pi0_P4", m_nPi0, 4, m_Pi0_P4);

      // eta
      status = m_tuple1->addItem("NumOfEta", m_nEta, 0, 200);
      status = m_tuple1->addIndexedItem("Eta_Gam1_ID", m_nEta, m_Eta_Gam1_ID);
      status = m_tuple1->addIndexedItem("Eta_Gam2_ID", m_nEta, m_Eta_Gam2_ID);
      status = m_tuple1->addIndexedItem("Eta_m", m_nEta, m_Eta);
      status = m_tuple1->addIndexedItem("Eta_P4", m_nEta, 4, m_Eta_P4);

      // Sigma+
      status = m_tuple1->addItem("NumOfSigmap", m_nSigmap, 0, 200);
      status = m_tuple1->addIndexedItem("Sigmap_Proton_ID", m_nSigmap, m_Sigmap_Proton_ID);
      status = m_tuple1->addIndexedItem("Sigmap_Pi0_ID", m_nSigmap, m_Sigmap_Pi0_ID);
      status = m_tuple1->addIndexedItem("Sigmap_m", m_nSigmap, m_Sigmap);

      // Lc
      status = m_tuple1->addItem("NumOfLc", m_nLc, 0, 1000);
      status = m_tuple1->addIndexedItem("Lc_Charge", m_nLc, m_Lc_Charge);
      // status = m_tuple1->addIndexedItem("Lc_Sigmap_ID", m_nLc, m_Lc_Sigmap_ID);
      // status = m_tuple1->addIndexedItem("Lc_Ks_ID", m_nLc, m_Lc_Ks_ID);
      status = m_tuple1->addIndexedItem("Lc_Mass", m_nLc, m_Lc_Mass); // no Recoil 1C update, no Beam constrain
      status = m_tuple1->addIndexedItem("Lc_MBC", m_nLc, m_Lc_MBC);   // no Recoil 1C update, Beam constrain
      status = m_tuple1->addIndexedItem("Lc_De", m_nLc, m_Lc_De);
      // Recoil 1C
      // status = m_tuple1->addIndexedItem("Lc_Chisq_1c", m_nLc, m_Lc_Chisq_1c);
      // status = m_tuple1->addIndexedItem("Lc_Proton_P4_1c", m_nLc, 4, m_Lc_Proton_P4_1c);
      // status = m_tuple1->addIndexedItem("Lc_Pi0_P4_1c", m_nLc, 4, m_Lc_Pi0_P4_1c);
      // status = m_tuple1->addIndexedItem("Lc_Ks_P4_1c", m_nLc, 4, m_Lc_Ks_P4_1c);
      // status = m_tuple1->addIndexedItem("Lc_P4_1c", m_nLc, 4, m_Lc_P4_1c);
      // status = m_tuple1->addIndexedItem("Lc_Mass_1c", m_nLc, m_Lc_Mass_1c); // Recoil 1C update, no Beam constrain
      // status = m_tuple1->addIndexedItem("Lc_MBC_1c", m_nLc, m_Lc_MBC_1c);   // Recoil 1C update, Beam constrain
      // status = m_tuple1->addIndexedItem("Lc_De_1c", m_nLc, m_Lc_De_1c);

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
      status = m_tuple1->addIndexedItem("gamb", m_p4index, m_gamb_p4);
      status = m_tuple1->addIndexedItem("gama", m_p4index, m_gama_p4);

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
        if (pdg == 221 && motherpdg == 4122)
        { //// eta->gam gam
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Ap_id[ndaughterAp] = gc[ii]->particleProperty();

            if (m_debug)
              cout << "4122 221 ndaughterAp[" << ndaughterAp << "]: " << Ap_id[ndaughterAp] << endl;

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
              cout << "-4122 221 ndaughterAm[" << ndaughterAm << "]: " << Ap_id[ndaughterAm] << endl;

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
              cout << "4122 221 ndaughterAm[" << ndaughterAm << "]: " << Ap_id[ndaughterAm] << endl;

            for (int ll = 0; ll < 4; ll++)
              Am_ptruth[ndaughterAm][ll] = gc[ii]->initialFourMomentum()[ll];
            ndaughterAm++;
          } // End of "gc.size() > 0" IF
        }
        if (pdg == 221 && motherpdg == -4122)
        { //// eta->gam gam
          const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList();
          for (unsigned int ii = 0; ii < gc.size(); ii++)
          {
            if (gc[ii]->particleProperty() == -22)
              continue;
            Am_id[ndaughterAm] = gc[ii]->particleProperty();

            if (m_debug)
              cout << "4122 221 ndaughterAm[" << ndaughterAm << "]: " << Ap_id[ndaughterAm] << endl;

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

      if (ndaughterAp == 8 && Ap_id[0] == 3222 && Ap_id[1] == 221 && Ap_id[2] == 2212 && Ap_id[3] == 111 &&
          Ap_id[4] == 22 && Ap_id[5] == 22 && Ap_id[6] == 22 && Ap_id[7] == 22)
      {
        signal = 1;
        all++;
        all_p++;
      }
      if (ndaughterAm == 8 && Am_id[0] == -3222 && Am_id[1] == 221 && Am_id[2] == -2212 && Am_id[3] == 111 &&
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

      WTrackParameter *wtrkp = new WTrackParameter(xmass[4], mdcKalTrk->getZHelixP(), mdcKalTrk->getZErrorP());
      HepLorentzVector p4 = mdcKalTrk->p4(xmass[4]);

      MyParticle tmp(goodTrack[i], p4,  mdcTrk->charge(), wtrkp);
      proton.push_back(tmp);
    }
  }

  if (abs(signal) == 1)
    Ncut0++;

  int np = trackProntonP.size();
  int npbar = trackProntonPbar.size();
  m_np = np;
  m_npbar = npbar;

  if (np == 0 && npbar == 0)
    return StatusCode::SUCCESS;

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
      cout << __LINE__ << " "
           << "emcGamma.size(): " << emcGamma.size() << endl;
    Ncut2++;
  }

  if (m_debug)
    cout << __LINE__ << endl;

#pragma region loop_gamma_check_eta_pi0______________________________________________________________

  // Vint igam1, igam2, igam3, igam4;
  // igam1.clear(), igam2.clear(), igam3.clear(), igam4.clear();

  // Vp4 pgam1, pgam2, pgam3, pgam4, pgam1_1C, pgam2_1C, pgam3_1C, pgam4_1C;
  // pgam1.clear(), pgam2.clear(), pgam3.clear(), pgam4.clear(), pgam1_1C.clear(), pgam2_1C.clear(), pgam3_1C.clear(),
  //     pgam4_1C.clear();

  // Vdouble chi2;
  // chi2.clear();

  std::vector<MyMotherParticleFit> eta;
  eta.clear();

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
      cout << __LINE__ << " 00000000 " << p2geta.px() << " " << p2geta.py() << " " << p2geta.pz() << " " << p2geta.e()
           << " " << p2geta.m() << endl;

      if (p2geta.m() < m_EtaMinMass || p2geta.m() > m_EtaMaxMass)
        continue;
      // if (m_test1C == 1)
      // {
      //   kmfit->init();
      //   kmfit->setChisqCut(m_chisqMax);

      //   kmfit->AddTrack(0, 0.0, emcTrki);
      //   kmfit->AddTrack(1, 0.0, emcTrkj);
      //   // 0.547862
      //   kmfit->AddResonance(0, 0.547862, 0, 1);
      //   bool oksq = kmfit->Fit();
      //   if (oksq)
      //   {
      //     kmfit->BuildVirtualParticle(0);

      MyMotherParticleFit tmp(emcGamma[i], emcGamma[j]);
      tmp.setLorentzVector(p2geta);
      eta.push_back(tmp);
      // }
      // }
    }
  }

  // Loop each gamma pair, check pi0 mass ---------------------
  std::vector<MyMotherParticleFit> pi0;
  pi0.clear();

  for (int k = 0; k < emcGamma.size() - 1; k++)
  {
    RecEmcShower *emcTrkk = emcGamma[k].getRecEmcShower();
    HepLorentzVector ptrkk = emcGamma[k].getLorentzVector();
    for (int l = k + 1; l < emcGamma.size(); l++)
    {
      RecEmcShower *emcTrkl = emcGamma[l].getRecEmcShower();
      HepLorentzVector ptrkl = emcGamma[l].getLorentzVector();

      HepLorentzVector p2gpi = ptrkk + ptrkl;

      if (p2gpi.m() < m_Pi0MinMass || p2gpi.m() > m_Pi0MaxMass)
        continue;
      // if (m_test1C == 1)
      // {
      //   kmfit->init();
      //   kmfit->setChisqCut(m_chisqMax);

      //   kmfit->AddTrack(0, 0.0, emcTrkk);
      //   kmfit->AddTrack(1, 0.0, emcTrkl);
      //   // 0.1349770
      //   kmfit->AddResonance(0, 0.1349770, 0, 1);
      //   bool oksq = kmfit->Fit();

      //   if (oksq)
      //   {
      //     kmfit->BuildVirtualParticle(0);
      // double chisq_1C = kmfit->chisq(0);
      // WTrackParameter Pi0WTrk = kmfit->wVirtualTrack(0);

      MyMotherParticleFit tmp(emcGamma[k], emcGamma[l]);
      tmp.setLorentzVector(p2gpi);
      pi0.push_back(tmp);
      // }
      // }
    }
  }
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

  // if (ngam12 == 0)
  //   return StatusCode::SUCCESS;
  // if (abs(signal) == 1)
  //   Ncut3++; // Ncut3 should equal Ncut2;

  if (pi0.size() == 0 || eta.size() == 0)
    return StatusCode::SUCCESS;

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
  // double deltaE_mina = 9999;
  // double deltaE_minb = 9999;
  // double deltaE_minc = 9999;
  // double deltaE_mind = 9999;
  // double chisq = 0;
  // int Km_index = -1;
  // int Kp_index = -1;
  // int p_index = -1;
  // int pbar_index = -1;
  // int gam1_index = -1;
  // int gam3_index = -1;
  // int gam1_indexm = -1;
  // int gam3_indexm = -1;
  // int gam2_index = -1;
  // int gam4_index = -1;
  // int gam2_indexm = -1;
  // int gam4_indexm = -1;
  // int pim_index = -1;
  // int pip_index = -1;
  // int a = 0;
  // int b = 0;
  // int c = 0;
  // int d = 0;
  // HepLorentzVector p_p4(0, 0, 0, 0), gam1a_p4(0, 0, 0, 0), gam2a_p4(0, 0, 0, 0), gam1b_p4(0, 0, 0, 0),
  //     gam2b_p4(0, 0, 0, 0), gam3a_p4(0, 0, 0, 0), gam4a_p4(0, 0, 0, 0), gam3b_p4(0, 0, 0, 0), gam4b_p4(0, 0, 0, 0),
  //     pip_p4(0, 0, 0, 0), pgam1a_1C4p(0, 0, 0, 0), pgam2a_1C4p(0, 0, 0, 0), pgam1b_1C4p(0, 0, 0, 0),
  //     pgam2b_1C4p(0, 0, 0, 0), pgam3a_1C4p(0, 0, 0, 0), pgam4a_1C4p(0, 0, 0, 0), pgam3b_1C4p(0, 0, 0, 0),
  //     pgam4b_1C4p(0, 0, 0, 0), pbar_p4(0, 0, 0, 0), pim_p4(0, 0, 0, 0);

#pragma region save_lambda_c - _and_lambda_c + ________________________________________________________________________

  // select lambda_c-, 34-> pi, 12->eta
  std::vector<MyMotherParticleFit> recoilLambdac_1c, lambdac, sigma;
  recoilLambdac_1c.clear();
  lambdac.clear();
  sigma.clear();
  for (int i = 0; i < proton.size(); i++)
  {
    for (int j = 0; j < eta.size(); j++)
    {
      for (int k = 0; k < pi0.size(); k++)
      {
        if (eta[j].getChild1().getIndex() == pi0[k].getChild1().getIndex() ||
            eta[j].getChild1().getIndex() == pi0[k].getChild2().getIndex())
          continue;
        if (eta[j].getChild2().getIndex() == pi0[k].getChild1().getIndex() ||
            eta[j].getChild2().getIndex() == pi0[k].getChild2().getIndex())
          continue;

        HepLorentzVector psigma = proton[i].getLorentzVector() + pi0[k].getChild1().getLorentzVector() +
                                  pi0[k].getChild1().getLorentzVector();
        if (psigma.m() < m_SigmaMinMass || psigma.m() > m_SigmaMaxMass)
          continue;
        MyMotherParticleFit tmp0(proton[i], pi0[k]);
        tmp0.setLorentzVector(psigma);
        sigma.push_back(tmp0);

        MyMotherParticleFit tmp(proton[i], pi0[k], eta[j]);
       
        HepLorentzVector pLambda = proton[i].getLorentzVector() + pi0[k].getLorentzVector() + eta[j].getLorentzVector();
        pLambda.boost(-m_beta);
        tmp.setLorentzVector(pLambda);
        lambdac.push_back(tmp);

        // 1C
        KalmanKinematicFit *kmfit1 = KalmanKinematicFit::instance();
        kmfit1->init();
        kmfit1->setChisqCut(1e3);
        kmfit1->setIterNumber(10);
        
        kmfit1->AddTrack(0, proton[i].getTrackParameter());
        kmfit1->AddTrack(1, pi0.getChild1().getRecEmcShower());
        kmfit1->AddTrack(2, pi0.getChild2().getRecEmcShower());
        kmfit1->AddTrack(3, eta.getChild1().getRecEmcShower());
        kmfit1->AddTrack(4, eta.getChild2().getRecEmcShower());

        kmfit1->AddMissTrack(5, 2.28646);
        kmfit1->AddFourMomentum(0, HepCMS);

        bool okvs1 = kmfit1->Fit();
        if (m_debug)
          std::cerr << "okvs_Lc_Recoil_1C_fit=" << okvs1 << std::endl;
        if (okvs1)
        {
          kmfit1->BuildVirtualParticle(0);
          // LcWTrk_1C = kmfit1->wVirtualTrack(0);
          // double Recoil_Chisq = kmfit1->chisq();

          MyMotherParticleFit tmp2(proton[i], pi0[k], eta[j], kmfit1);

          recoilLambdac_1c.push_back(tmp2);

          // proton_1C = kmfit1->pfit(0);
          // pi0_1C = kmfit1->pfit(1);
          // ks_1C = kmfit1->pfit(2);
          // Lc_1C = proton_1C + pi0_1C + ks_1C;

          // double mLc_1C = Lc_1C.m(); // Recoil 1C update, no Beam constrain
          // Lc_1C.boost(-m_beta);
          // double dE_Lc_1C = Lc_1C.t() - ebeam;
          // double mLc_1C_bc2 = ebeam * ebeam - Lc_1C.v().mag2();
          // double mLc_1C_bc = mLc_1C_bc2 > 0 ? sqrt(mLc_1C_bc2) : -10; // Recoil 1C update, Beam constrain

          // tmp_lc.chisq1 = Recoil_Chisq;
          // tmp_lc.p4_proton_1c = proton_1C;
          // tmp_lc.p4_pi0_1c = pi0_1C;
          // tmp_lc.p4_ks_1c = ks_1C;
          // tmp_lc.p4_Lc_1c = Lc_1C;
          // tmp_lc.m_1c = mLc_1C;
          // tmp_lc.m_bc_1c = mLc_1C_bc;
          // tmp_lc.m_de_1c = dE_Lc_1C;
        }

        /*/   if (fabs(deltaEb) < fabs(deltaE_minb))
        //   {
        //     b = 1;
        //     deltaE_minb = deltaEb;

        //     pbar_index = trackProntonPbar[i];
        //     gam3_indexm = igam3[k];
        //     gam4_indexm = igam4[k];
        //     pbar_p4 = ppbar[i];
        //     gam3b_p4 = pgam3[k];
        //     gam4b_p4 = pgam4[k];
        //     gam1b_p4 = pgam1[j];
        //     gam2b_p4 = pgam2[j];
        //     chisq = chi2[k];
        //     pgam3b_1C4p = pgam3_1C[k];
        //     pgam4b_1C4p = pgam4_1C[k];
        //     pgam1b_1C4p = pgam1_1C[j];
        //     pgam2b_1C4p = pgam2_1C[j];
        //   }*/
      }
    }
  }
  /*/ cout << __LINE__ << (gam3b_p4 + gam4b_p4).m() << endl;
  // cout << __LINE__ << " 00000000 " << gam3b_p4.px() << " " << gam3b_p4.py() << " " << gam3b_p4.pz() << " "
  //      << gam3b_p4.e() << " " << gam3b_p4.m() << endl;

  // cout << __LINE__ << " 00000000 " << gam4b_p4.px() << " " << gam4b_p4.py() << " " << gam4b_p4.pz() << " "
  //      << gam4b_p4.e() << " " << gam4b_p4.m() << endl;

  // // save lambda_c+ a ....
  // for (int i = 0; i < np; i++)
  // {
  //   for (int j = 0; j < ngam12; j++)
  //   {
  //     for (int k = 0; k < ngam34; k++)
  //     {
  //       HepLorentzVector psigma = pp[i] + pgam3_1C[k] + pgam4_1C[k];

  //       if (psigma.m() < m_SigmaMinMass || psigma.m() > m_SigmaMaxMass)
  //         continue;
  //       if (igam1[j] == igam3[k] || igam1[j] == igam4[k])
  //         continue;
  //       if (igam2[j] == igam3[k] || igam2[j] == igam4[k])
  //         continue;

  //       HepLorentzVector pLambda = pp[i] + pgam3_1C[k] + pgam4_1C[k] + pgam1_1C[j] + pgam2_1C[j];
  //       pLambda.boost(-m_beta);
  //       // double deltaE = fabs(pLambda.t() - ebeam);
  //       double deltaEa = pLambda.t() - ebeam;
  //       if (m_debug)
  //         cout << "fabs(deltaEa): " << fabs(deltaEa) << ", fabs(deltaE_mina): " << fabs(deltaE_mina) << endl;

  //       if (fabs(deltaEa) < fabs(deltaE_mina))
  //       {
  //         a = 1;
  //         deltaE_mina = deltaEa;

  //         p_index = trackProntonP[i];
  //         gam3_index = igam3[k];
  //         gam4_index = igam4[k];
  //         p_p4 = pp[i];
  //         gam3a_p4 = pgam3[k];
  //         gam4a_p4 = pgam4[k];
  //         gam1a_p4 = pgam1[j];
  //         gam2a_p4 = pgam2[j];
  //         chisq = chi2[k];
  //         pgam3a_1C4p = pgam3_1C[k];
  //         pgam4a_1C4p = pgam4_1C[k];
  //         pgam1a_1C4p = pgam1_1C[j];
  //         pgam2a_1C4p = pgam2_1C[j];
  //       }
  //     }
  //   }
  // }*/
#pragma endregion

#pragma region write__________________________________________________________________
  // if (a == 0 && b == 0)
  // {
  //   if (m_debug)
  //     cout << "a=0; b=0; "
  //          << "np: " << np << ", npbar: " << npbar << endl;
  //   return StatusCode::SUCCESS;
}
if (true)
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

  // Proton
  m_nProton = proton.size();
  for (int i = 0; i < proton.size(); i++)
  {
    m_Proton_ID[i] = proton[i].getIndex();
    m_Proton_Charge[i] = proton[i].getCharge();
    for (int j = 0; j < 4; j++)
    {
      m_Proton_P4[i][j] = proton[i].getLorentzVector()[j];
    }
    // for (int j = 0; j < 5; j++)
    // {
    //   m_Proton_PID[i][j] = proton[i].PidProb[j];
    // }
  }

  // Pi0 and 2 gamma
  m_nPi0 = pi0.size();
  for (int i = 0; i < pi0.size(); i++)
  {
    m_Pi0_Gam1_ID[i] = pi0[i].getChild1().getIndex();
    m_Pi0_Gam2_ID[i] = pi0[i].getChild2().getIndex();
    m_Pi0[i] = pi0.getMass();

    for (int j = 0; j < 4; j++)
    {
      m_Pi0_P4[i][j] = pi0[i].getLorentzVector()[j];
      // m_Pi0_Gam2_P4_1c[i][j] = pi0[i].pi0[i].getChild2().getLorentzVector()[j];
      // m_Pi0_P4_1c[i][j] = pi0[i].pi0_p4_1C[j];
    }
  }

  // eta and 2 gamma
  m_neta = eta.size();
  for (int i = 0; i < eta.size(); i++)
  {
    m_Eta_Gam1_ID[i] = eta[i].getChild1().getIndex();
    m_Eta_Gam2_ID[i] = eta[i].getChild2().getIndex();
    m_Eta[i] = eta.getMass();
    // m_eta_Chisq[i] = pi0[i].getFit()->chisq(0);

    for (int j = 0; j < 4; j++)
    {
      m_Eta_P4[i][j] = eta[i].getLorentzVector()[j];
    }
  }

  // Sigma+ and p
  m_nSigmap = sigma.size();5l
  for (int i = 0; i < sigma.size(); i++)
  {
    m_Sigmap_Proton_ID[i] = sigma[i].getChild1().getIndex();
    m_Sigmap_Pi0_ID[i] = sigma[i].getChild2().getIndex();
    m_Sigmap[i] = sigma[i].getMass();
  }

   m_nLc = lambdac.size();
  for (int i = 0; i < lambdac.size(); i++)
  {
    m_Lc_Charge[i] = lambdac[i].getChild1().getCharge();
    // m_Lc_Sigmap_ID[i] = lambdac[i].sigma_id;
    // m_Lc_Ks_ID[i] = lambdac[i].ks_id;

    m_Lc_Mass[i] = lambdac[i].getMass();
    HepLorentzVector pLambda = lambdac[i].getLorentzVector();
    double deltaEb = pLambda.t() - ebeam;
    double mbc2 = ebeam * ebeam - pLambda.v().mag2();
    m_bc = mbc2 > 0 ? sqrt(mbc2) : -10;
    m_Lc_MBC[i] = m_bc;
    m_Lc_De[i] = deltaEb;
  }

  m_nLc = lambdac.size();
  // cout << 
  // for (int i = 0; i < lambdac.size(); i++)
  // {
  //   m_Lc_Chisq_1c[i] = lambdac[i].getFit()->chisq();
  //   HepLorentzVector pLambda_1c = lambdac[i].getLorentzVector();
  //   double deltaEb = pLambda.t() - ebeam;
  //   double mbc2 = ebeam * ebeam - pLambda.v().mag2();
  //   m_bc = mbc2 > 0 ? sqrt(mbc2) : -10;

  //   m_Lc_Mass_1c[i] = lambdac[i].m_1c;
  //   m_Lc_MBC_1c[i] = lc_info[i].m_bc_1c;
  //   m_Lc_De_1c[i] = lc_info[i].m_de_1c;

  //   // m_Lc_Chisq_2c[i] = lc_info[i].chisq2;
  //   // m_Lc_Mass_2c[i] = lc_info[i].m_2c;
  //   // m_Lc_MBC_2c[i] = lc_info[i].m_bc_2c;
  //   // m_Lc_De_2c[i] = lc_info[i].m_de_2c;
  // }

    // }

    // m_bg = bg;
    // for (int jj = 0; jj < 4; jj++)
    //   m_pbar_p4[jj] = pbar_p4[jj];
    // for (int jj = 0; jj < 4; jj++)
    //   m_gam3b_p4[jj] = gam3b_p4[jj];
    // for (int jj = 0; jj < 4; jj++)
    //   m_gam4b_p4[jj] = gam4b_p4[jj];
    // for (int jj = 0; jj < 4; jj++)
    //   m_gam1b_p4[jj] = gam1b_p4[jj];
    // for (int jj = 0; jj < 4; jj++)
    //   m_gam2b_p4[jj] = gam2b_p4[jj];

    // m_pbarindex = pbar_index;
    // m_p4index = 4;
    // m_chi2 = chisq;
    // m_pi0m1c = (pgam3b_1C4p + pgam4b_1C4p).m();
    // m_etam1c = (pgam1b_1C4p + pgam2b_1C4p).m();

    // m_pi0m = (gam3b_p4 + gam4b_p4).m();
    // m_etam = (gam1b_p4 + gam2b_p4).m();
    // m_ebeam = ebeam;
    // m_Sigmam = (pgam3b_1C4p + pgam4b_1C4p + pbar_p4).m();
    // //	cout<<"deltaE="<<deltaE_minb<<endl;
    // m_deltaE_min = deltaE_minb;
    // //	cout<<"deltaE_min="<<m_deltaE_min<<endl;
    // HepLorentzVector pLambda = pbar_p4 + pgam3b_1C4p + pgam4b_1C4p + pgam1b_1C4p + pgam2b_1C4p;
    // pLambda.boost(-m_beta);
    // double mbc2 = ebeam * ebeam - pLambda.v().mag2();

    // m_bc = mbc2 > 0 ? sqrt(mbc2) : -10;
    // m_rightflag = 2;
    // m_np = np;
    // m_npbar = npbar;
    m_tuple1->write();
  }
/*/ if (a == 1)
// {

//   m_mode1 = mm_mode1;
//   m_mode2 = mm_mode2;
//   m_mode3 = mm_mode3;
//   m_idxmc = numParticle;
//   for (int i = 0; i < numParticle; i++)
//   {
//     m_pdgid[i] = M_pdgid[i];
//     m_motheridx[i] = M_motheridx[i];
//   }

//   m_ndaughterAp = ndaughterAp;
//   for (int aa = 0; aa < ndaughterAp; aa++)
//     m_Ap_id[aa] = Ap_id[aa];
//   for (int aa = 0; aa < ndaughterAp; aa++)
//     for (int ll = 0; ll < 4; ll++)
//       m_Ap_ptruth[aa][ll] = Ap_ptruth[aa][ll];

//   m_ndaughterAm = ndaughterAm;
//   for (int aa = 0; aa < ndaughterAm; aa++)
//     m_Am_id[aa] = Am_id[aa];
//   for (int aa = 0; aa < ndaughterAm; aa++)
//     for (int ll = 0; ll < 4; ll++)
//       m_Am_ptruth[aa][ll] = Am_ptruth[aa][ll];

//   m_mcparticle_p = numParticle_p;
//   m_mcparticle_m = numParticle_m;
//   for (int i = 0; i < numParticle_p; i++)
//   {
//     m_pdgid_p[i] = M_pdgid_p[i];
//     m_motheridx_p[i] = M_motheridx_p[i];
//   }
//   for (int i = 0; i < numParticle_m; i++)
//   {
//     m_pdgid_m[i] = M_pdgid_m[i];
//     m_motheridx_m[i] = M_motheridx_m[i];
//   }

//   m_signal = signal;
//   m_bg = bg;
//   for (int jj = 0; jj < 4; jj++)
//     m_p_p4[jj] = p_p4[jj];
//   for (int jj = 0; jj < 4; jj++)
//     m_gam3a_p4[jj] = gam3a_p4[jj];
//   for (int jj = 0; jj < 4; jj++)
//     m_gam4a_p4[jj] = gam4a_p4[jj];
//   for (int jj = 0; jj < 4; jj++)
//     m_gam1a_p4[jj] = gam1a_p4[jj];
//   for (int jj = 0; jj < 4; jj++)
//     m_gam2a_p4[jj] = gam2a_p4[jj];

//   m_pindex = p_index;
//   m_p4index = 4;
//   m_chi2 = chisq;
//   m_pi0m1c = (pgam3a_1C4p + pgam4a_1C4p).m();
//   m_etam1c = (pgam1a_1C4p + pgam2a_1C4p).m();

//   m_pi0m = (gam3a_p4 + gam4a_p4).m();
//   m_etam = (gam1a_p4 + gam2a_p4).m();
//   m_ebeam = ebeam;
//   m_Sigmam = (pgam3a_1C4p + pgam4a_1C4p + p_p4).m();
//   m_deltaE_min = deltaE_mina;
//   //	cout<<"deltaE="<<deltaE_mina<<endl;
//   HepLorentzVector pLambda = p_p4 + pgam3a_1C4p + pgam4a_1C4p + pgam1a_1C4p + pgam2a_1C4p;
//   pLambda.boost(-m_beta);
//   double mbc2 = ebeam * ebeam - pLambda.v().mag2();
//   m_bc = mbc2 > 0 ? sqrt(mbc2) : -10;
//   m_rightflag = 1;
//   m_np = np;
//   m_npbar = npbar;
//   m_tuple1->write();
// }*/
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
  cout << "all       " << all << endl;
  //	cout << "al       " <<al << endl;
  //	cout << "yes     " << yes << endl;
  //	cout << "Npp>0:  " << Npp << endl;
  //	cout << "Npm>0:  " << Npm << endl;

  return StatusCode::SUCCESS;
}

// add your code here,for other member-functions
