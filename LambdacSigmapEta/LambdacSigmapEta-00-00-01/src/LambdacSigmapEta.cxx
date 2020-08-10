//This analysis algorithm is for decay:
//Lc+ -> Sigma+ eta
//Sigma+ -> p+ pi0
//eta -> gamma gamma

#include "LambdacSigmapEta/LambdacSigmapEta.h"
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
#include "SimplePIDSvc/ISimplePIDSvc.h"
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
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include "CLHEP/Geometry/Point3D.h"
#include "SimplePIDSvc/SimplePIDSvc.h"

#include "McDecayModeSvc/McDecayModeSvc.h"

#include "GaudiKernel/StatusCode.h"

#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include <vector>
//this is the mass of                            pi                   p
const double xmass[5] = {0.000511, 0.105658, 0.139570, 0.493677, 0.9382723}; // 0.9382723 proton
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef struct
{
    int index;
    double ang_antiproton;
} Gam_info;

class group
{
public:
    //chisquare of 4C
    double chisquare_3C;
    double chisquare_4C;
    //p4 after 4C
    HepLorentzVector p4_gam1_pi0_sigma_4C, p4_gam2_pi0_sigma_4C, p4_proton_sigma_4C;
    HepLorentzVector p4_gam1_eta_4C, p4_gam2_eta_4C;
    //p4 after 3C
    HepLorentzVector p4_gam1_pi0_sigma_3C, p4_gam2_pi0_sigma_3C;
    double msigma;
    //p4 before 3C
    HepLorentzVector p4_gam1_pi0_sigma_old, p4_gam2_pi0_sigma_old, p4_proton_sigma_old;
    HepLorentzVector p4_gam1_eta_old, p4_gam2_eta_old;
    //the minimum angle of all gamma contained
    double angle_gamma_antiproton;
};

int Ntotal(0), Ncut0(0), Ncut1(0), Ncut2(0), Ncut3(0), Ncut4(0), Ncut5(0), Ncut6(0), Ncut7(0), Ncut8(0), Ncut9(0);
int Ncharge(0), Nproton(0), Nanti_proton(0), Ngamma(0), Nlambdac(0), Nanti_lambdac(0);
int i, j, k, runNo, eventNo;

HepLorentzVector LambdacSigmapEta::getP4(RecEmcShower *gTrk, Hep3Vector origin)
{
    Hep3Vector Gm_Vec(gTrk->x(), gTrk->y(), gTrk->z());
    Hep3Vector Gm_Mom = Gm_Vec - origin;
    Gm_Mom.setMag(gTrk->energy());
    HepLorentzVector pGm(Gm_Mom, gTrk->energy());
    return pGm;
}

LambdacSigmapEta::LambdacSigmapEta(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
    declareProperty("Vr0cut", m_pvr0cut = 1.0);
    declareProperty("Vz0cut", m_pvz0cut = 10.0);
    declareProperty("Costheta", m_costheta = 0.93);
    declareProperty("PhotonMinEnergy", m_minEnergy = 0.025);
    declareProperty("GammathCut", m_gammathCut = 14.0);
    declareProperty("GammatlCut", m_gammatlCut = 0.0);
    declareProperty("PhotonMaxCosThetaBarrel", m_maxCosThetaBarrel = 0.8);
    declareProperty("PhotonMinCosThetaEndcap", m_minCosThetaEndcap = 0.86);
    declareProperty("PhotonMaxCosThetaEndcap", m_maxCosThetaEndcap = 0.92);
    declareProperty("PhotonMinEndcapEnergy", m_minEndcapEnergy = 0.050);

    declareProperty("SigmaMinMassCut", m_SigmaminMass = 1.10); //tight cut 1.176 , loose cut 1.10
    declareProperty("SigmaMaxMassCut", m_SigmamaxMass = 1.27); //tight cut 1.20 , loose cut 1.27
    declareProperty("SigmaMassFit", m_SigmaMassFit = 1.18937);

    declareProperty("Pi0MinMassCut", m_Pi0minMass = 0.08);  //tight cut 0.115 , loose cut 0.100
    declareProperty("Pi0MaxMassCut", m_Pi0maxMass = 0.180); //tight cut 0.150 , loose cut 0.180
    declareProperty("Pi0MassFit", m_Pi0MassFit = 0.134977);

    declareProperty("EtaMinMassCut", m_EtaminMass = 0.46); //tight cut  , loose cut 0.46
    declareProperty("EtaMaxMassCut", m_EtamaxMass = 0.60); //tight cut  , loose cut 0.60
    declareProperty("EtaMassFit", m_EtaMassFit = 0.547862);
    declareProperty("Debug", m_debug = true);
    declareProperty("GetMCInfo", m_getmcinfo = true);
    declareProperty("BeamE", m_beamE = 2.3);
    declareProperty("UseOverallTof", m_use_Total_TOF = false);
    declareProperty("ChisqMax", m_chisqMax = 200);
}
LambdacSigmapEta::~LambdacSigmapEta()
{
    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!!m_debug begin deconstructor !!!!!!!!!!" << std::endl;
    //add your code for deconstructor
}
StatusCode LambdacSigmapEta::initialize()
{
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "LambdacSigmapEta::initialize()" << endreq;
    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!!m_debug begin initialize !!!!!!!!!!" << std::endl;
    m_beta.setX(0.011);
    m_beta.setY(0);
    m_beta.setZ(0);
    //add your code here
    StatusCode status;
    NTuplePtr nt1(ntupleSvc(), "FILE1/SigmapEta");
    if (nt1)
        m_tuple1 = nt1;
    else
    {
        m_tuple1 = ntupleSvc()->book("FILE1/SigmapEta", CLID_ColumnWiseTuple, "Lambdac to SigmapEta");
        if (m_tuple1)
        {
            status = m_tuple1->addItem("run", m_run);

            status = m_tuple1->addItem("M_pi0", m_pi0);
            status = m_tuple1->addItem("M_eta", m_eta);
            status = m_tuple1->addItem("M_sigmap", m_sigmap);

            status = m_tuple1->addItem("event", m_event);
            status = m_tuple1->addItem("rightflag", m_rightflag);

            status = m_tuple1->addItem("isig", m_isig, 0, 2000);
            status = m_tuple1->addIndexedItem("charge", m_isig, m_charge);
            status = m_tuple1->addIndexedItem("chi2_3C", m_isig, m_chisq_3C);
            status = m_tuple1->addIndexedItem("chi2_4C", m_isig, m_chisq_4C); //for choose best candidate
            status = m_tuple1->addIndexedItem("angle", m_isig, m_angle);
            status = m_tuple1->addIndexedItem("deltaE_min", m_isig, m_deltaE_min);
            status = m_tuple1->addIndexedItem("M_BC", m_isig, m_bc);
            //pi0
            // status = m_tuple1->addIndexedItem("p4_of_Gamma1_Pi0_before_4C", m_isig, 4, m_Gam1_Pi0_p4_old); //cut pi0 mass
            // status = m_tuple1->addIndexedItem("p4_of_Gamma2_Pi0_before_4C", m_isig, 4, m_Gam2_Pi0_p4_old); //cut pi0 mass
            // status = m_tuple1->addIndexedItem("p4_of_Gamma1_Pi0_after_4C", m_isig, 4, m_Gam1_Pi0_p4);      //for mbc
            // status = m_tuple1->addIndexedItem("p4_of_Gamma2_Pi0_after_4C", m_isig, 4, m_Gam2_Pi0_p4);      //for mbc
            //eta
            status = m_tuple1->addIndexedItem("p4_of_Gamma1_Eta_before_4C", m_isig, 4, m_Gam1_Eta_p4_old); //cut eta mass
            status = m_tuple1->addIndexedItem("p4_of_Gamma2_Eta_before_4C", m_isig, 4, m_Gam2_Eta_p4_old); //cut eta mass
            status = m_tuple1->addIndexedItem("p4_of_Gamma1_Eta_after_4C", m_isig, 4, m_Gam1_Eta_p4);      //for mbc
            status = m_tuple1->addIndexedItem("p4_of_Gamma2_Eta_after_4C", m_isig, 4, m_Gam2_Eta_p4);      //for mbc
            //sigma+
            status = m_tuple1->addIndexedItem("p4_of_gamma1_pi0_Sigma_after_4C", m_isig, 4, m_gam1_pi0_Sigma_p4); //for mbc
            status = m_tuple1->addIndexedItem("p4_of_gamma2_pi0_Sigma_after_4C", m_isig, 4, m_gam2_pi0_Sigma_p4); //for mbc
            status = m_tuple1->addIndexedItem("p4_of_proton_Sigma_after_4C", m_isig, 4, m_p_Sigma_p4);            //for mbc

            status = m_tuple1->addIndexedItem("p4_of_gamma1_pi0_Sigma_before_4C", m_isig, 4, m_gam1_pi0_Sigma_p4_old); //cut pi0 mass
            status = m_tuple1->addIndexedItem("p4_of_gamma2_pi0_Sigma_before_4C", m_isig, 4, m_gam2_pi0_Sigma_p4_old); //cut pi0 mass
            status = m_tuple1->addIndexedItem("p4_of_proton_Sigma_before_4C", m_isig, 4, m_p_Sigma_p4_old);            //cut sigma mass
            status = m_tuple1->addIndexedItem("p4_of_gamma1_pi0_Sigma_after_3C", m_isig, 4, m_gam1_pi0_Sigma_p4_3C);   //cut sigma mass
            status = m_tuple1->addIndexedItem("p4_of_gamma2_pi0_Sigma_after_3C", m_isig, 4, m_gam2_pi0_Sigma_p4_3C);   //cut sigma mass

            if (m_getmcinfo)
            {
                status = m_tuple1->addItem("mode1", m_mode1);
                status = m_tuple1->addItem("mode2", m_mode2);
                status = m_tuple1->addItem("mode3", m_mode3);

                status = m_tuple1->addItem("Lmdc_P", m_Lmdc_P);
                status = m_tuple1->addItem("Lmdc_M", m_Lmdc_M);

                status = m_tuple1->addItem("indexmc", m_idxmc, 0, 2000);
                status = m_tuple1->addIndexedItem("pdgid", m_idxmc, m_pdgid);
                status = m_tuple1->addIndexedItem("motheridx", m_idxmc, m_motheridx);

                status = m_tuple1->addItem("indexmc_p", m_mcparticle_p, 0, 2000);
                status = m_tuple1->addIndexedItem("pdgid_p", m_mcparticle_p, m_pdgid_p);
                status = m_tuple1->addIndexedItem("motheridx_p", m_mcparticle_p, m_motheridx_p);
                status = m_tuple1->addItem("indexmc_m", m_mcparticle_m, 0, 2000);
                status = m_tuple1->addIndexedItem("pdgid_m", m_mcparticle_m, m_pdgid_m);
                status = m_tuple1->addIndexedItem("motheridx_m", m_mcparticle_m, m_motheridx_m);

                status = m_tuple1->addItem("ndaughterAp", m_ndaughterAp, 0, 15);
                status = m_tuple1->addIndexedItem("Ap_id", m_ndaughterAp, m_Ap_id);
                status = m_tuple1->addIndexedItem("Ap_ptruth", m_ndaughterAp, 4, m_Ap_ptruth);

                status = m_tuple1->addItem("ndaughterAm", m_ndaughterAm, 0, 15);
                status = m_tuple1->addIndexedItem("Am_id", m_ndaughterAm, m_Am_id);
                status = m_tuple1->addIndexedItem("Am_ptruth", m_ndaughterAm, 4, m_Am_ptruth);
            }
        }
        else
        {
            log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
            return StatusCode::FAILURE;
        }
    }
    if (m_getmcinfo)
    {
        NTuplePtr nt2(ntupleSvc(), "FILE1/SigmapEta_truth");
        if (nt2)
            m_tuple2 = nt2;
        else
        {
            m_tuple2 = ntupleSvc()->book("FILE1/SigmapEta_truth", CLID_ColumnWiseTuple, "Lambdac to SigmapEta");
            if (m_tuple2)
            {
                status = m_tuple2->addItem("run", m_run_t);
                status = m_tuple2->addItem("event", m_event_t);
                status = m_tuple2->addItem("mode1", m_mode1_t);
                status = m_tuple2->addItem("mode2", m_mode2_t);
                status = m_tuple2->addItem("mode3", m_mode3_t);
                status = m_tuple2->addItem("Lmdc_P", m_Lmdc_P_t);
                status = m_tuple2->addItem("Lmdc_M", m_Lmdc_M_t);
                status = m_tuple2->addItem("Flag", m_Flag_t);
                status = m_tuple2->addItem("ndaughterAp", m_ndaughterAp_t, 0, 15);
                status = m_tuple2->addIndexedItem("Ap_id", m_ndaughterAp_t, m_Ap_id_t);
                status = m_tuple2->addIndexedItem("Ap_ptruth", m_ndaughterAp_t, 4, m_Ap_ptruth_t);
                status = m_tuple2->addItem("ndaughterAm", m_ndaughterAm_t, 0, 15);
                status = m_tuple2->addIndexedItem("Am_id", m_ndaughterAm_t, m_Am_id_t);
                status = m_tuple2->addIndexedItem("Am_ptruth", m_ndaughterAm_t, 4, m_Am_ptruth_t);
            }
            else
            {
                log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple2) << endmsg;
                return StatusCode::FAILURE;
            }
        }
    }
    return StatusCode::SUCCESS;
}

StatusCode LambdacSigmapEta::beginRun()
{
    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!!m_debug begin beginRun !!!!!!!!!!" << std::endl;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "LambdacSigmapEta::beginRun()" << endreq;
    //add your code here
    return StatusCode::SUCCESS;
}

bool gammaidentifier(Vint a)
{
    int x = 1;
    for (int i = 0; i < a.size() - 1; i++)
    {
        for (int j = i + 1; j < a.size(); j++)
        {
            if (a[i] == a[j])
                x = 0;
        }
    }
    if (x == 1)
        return true;
    else
        return false;
}

StatusCode LambdacSigmapEta::execute()
{
    m_rightflag = -100;
    m_run = 5000;
    m_event = -10;

    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!!m_debug begin execute !!!!!!!!!!" << std::endl;
    Ntotal++;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "LambdacSigmapEta::execute()" << endreq;
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");

    runNo = eventHeader->runNumber();
    eventNo = eventHeader->eventNumber();
    if (m_debug)
        std::cerr << __LINE__ << "runNo=" << runNo << "  ,  eventNo=" << eventNo << std::endl;

    log << MSG::DEBUG << "run, evtnum = " << runNo << " , " << eventNo << endreq;

    IMcDecayModeSvc *i_svc;
    StatusCode sc_DecayModeSvc = service("McDecayModeSvc", i_svc);
    if (sc_DecayModeSvc.isFailure())
    {
        log << MSG::FATAL << "Could not load McDecayModeSvc!" << endreq;
        return sc_DecayModeSvc;
    }
    m_svc = dynamic_cast<McDecayModeSvc *>(i_svc);

    int mode1 = ((eventHeader->flag1() / 1000000)) % 1000;
    int mode2 = (eventHeader->flag1() / 1000) % 1000;
    int mode3 = eventHeader->flag1() % 1000;

    int numParticle_p = 0;
    int numParticle_m = 0;
    int numParticle = 0;

    int M_pdgid_p[2000];
    int M_motheridx_p[2000];

    int M_pdgid_m[2000];
    int M_motheridx_m[2000];

    int M_pdgid[2000];
    int M_motheridx[2000];

    //MC info
    bool Lmdc_P = false;
    bool Lmdc_M = false;

    int ndaughterAp = 0;
    int Ap_id[15];
    double Ap_ptruth[15][4];
    // initialize the array with 0;
    for (int aa = 0; aa < 15; aa++)
        for (int ll = 0; ll < 4; ll++)
            Ap_ptruth[aa][ll] = 0;
    for (int aa = 0; aa < 15; aa++)
        Ap_id[aa] = 0;

    int ndaughterAm = 0;
    int Am_id[15];
    double Am_ptruth[15][4];
    // initialize the array with 0;
    for (int aa = 0; aa < 15; aa++)
        for (int ll = 0; ll < 4; ll++)
            Am_ptruth[aa][ll] = 0;
    for (int aa = 0; aa < 15; aa++)
        Am_id[aa] = 0;

    int Flag = 0;

    if (m_debug)
        std::cerr << __LINE__ << " mode1 = " << mode1 << " mode2 = " << mode2 << " mode3 = " << mode3 << std::endl;
    if (m_debug)
        std::cerr << __LINE__ << " begin get truth info" << std::endl;
    if (eventHeader->runNumber() < 0)
    {
        SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
        if (!mcParticleCol)
        {
            std::cerr << "Could not retrieve McParticelCol" << std::endl;
            return StatusCode::FAILURE;
        }
        else
        {
            if (m_debug)
                std::cerr << __LINE__ << " begin get pdg number" << std::endl;
            std::vector<int> pdgid;
            std::vector<int> motherindex;
            pdgid.clear();
            motherindex.clear();

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
                int mmmotherpdg = ((((*iter_mc)->mother()).mother()).mother()).particleProperty();
                { //test
                    if (pdg == 2212)
                    {
                        Flag = 1;
                    }
                    if (pdg == -2212)
                    {
                        Flag = -1;
                    }
                }
                if (m_debug)
                    std::cerr << __LINE__ << "pdg:" << pdg << "  mpgd:" << motherpdg << "  mmpgd:" << mmotherpdg << "  mmmpgd:" << mmmotherpdg << std::endl;

                if (pdg == 4122 && motherpdg == 9030443)
                { // psi -> Lambda_c+
                    Lmdc_P = true;
                    const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList(); //lambdac+ -> sigma+ pi0 pi0
                    for (unsigned int ii = 0; ii < gc.size(); ii++)
                    {
                        // 22 is gamma
                        if (gc[ii]->particleProperty() == -22)
                            continue; //why -22?
                        Ap_id[ndaughterAp] = gc[ii]->particleProperty();
                        for (int ll = 0; ll < 4; ll++)
                            Ap_ptruth[ndaughterAp][ll] = gc[ii]->initialFourMomentum()[ll];

                        ndaughterAp++;
                    } // End of "gc.size() > 0" IF
                }

                if (pdg == 3222 && motherpdg == 4122 && mmotherpdg == 9030443)
                { // psi -> Lambda_c+ -> Sigma+
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

                if (pdg == 221 && motherpdg == 4122 && mmotherpdg == 9030443)
                { // psi -> Lambda_c+ -> eta
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

                if (pdg == 111 && motherpdg == 3222 && mmotherpdg == 4122 && mmmotherpdg == 9030443)
                { // psi -> Lambda_c+ -> Sigma+ -> pi0
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

                if (pdg == -4122 && motherpdg == 9030443)
                { // psi -> Lambda_c-
                    Lmdc_M = true;
                    const SmartRefVector<Event::McParticle> &gc = (*iter_mc)->daughterList(); //Lambdac- ->sigma- pi0 pi0
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

                if (pdg == -3222 && motherpdg == -4122 && mmotherpdg == 9030443)
                { // psi -> Lambda_c- -> Sigma-
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

                if (pdg == 221 && motherpdg == -4122 && mmotherpdg == 9030443)
                { // psi -> Lambda_c- -> eta
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

                if (pdg == 111 && motherpdg == -3222 && mmotherpdg == -4122 && mmmotherpdg == 9030443)
                { // psi -> Lambda_c- -> Sigma- -> pi0
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

            for (iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++)
            {
                if ((*iter_mc)->primaryParticle())
                    continue;
                if (!(*iter_mc)->decayFromGenerator())
                    continue;

                if ((*iter_mc)->particleProperty() == -4122)
                { //4122 Lambdac+    -4122 Lambdac-
                    int mode = m_svc->extract(*iter_mc, pdgid, motherindex);
                    numParticle = pdgid.size();
                    for (int i = 0; i != pdgid.size(); i++)
                    {
                        M_pdgid[i] = pdgid[i];
                        //	if(m_debug) cout<<"M_pdgid[i]="<<M_pdgid[i]<<endl;
                        M_motheridx[i] = motherindex[i];
                        //	if(m_debug) cout<<"M_motheridx[i]="<<M_motheridx[i]<<endl;
                    }
                }
                pdgid.clear();
                motherindex.clear();

                //trace Lambda_C+
                if ((*iter_mc)->particleProperty() == 4122)
                {
                    int mode = m_svc->extract(*iter_mc, pdgid, motherindex);
                    numParticle_p = pdgid.size();
                    for (int i = 0; i != pdgid.size(); i++)
                    {
                        M_pdgid_p[i] = pdgid[i];
                        //						if(m_debug) cout<<"M_pdgid_p[i]="<<M_pdgid_p[i]<<endl;
                        M_motheridx_p[i] = motherindex[i];
                        //						if(m_debug) cout<<"M_motheridx_p[i]="<<M_motheridx_p[i]<<endl;
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
                        //						if(m_debug) cout<<"M_pdgid_m[i]="<<M_pdgid_m[i]<<endl;
                        M_motheridx_m[i] = motherindex[i];
                        //						if(m_debug) cout<<"M_motheridx_m[i]="<<M_motheridx_m[i]<<endl;
                    }
                }
            }

            if (m_debug)
                cerr << "Lmdc_P=" << Lmdc_P << " Lmdc_M=" << Lmdc_M << endl;
            m_run_t = runNo;
            m_event_t = eventNo;
            m_mode1_t = mode1;
            m_mode2_t = mode2;
            m_mode3_t = mode3;
            m_Lmdc_P_t = Lmdc_P;
            m_Lmdc_M_t = Lmdc_M;
            m_Flag_t = Flag;

            m_ndaughterAp_t = ndaughterAp;
            for (i = 0; i < ndaughterAp; i++)
                m_Ap_id_t[i] = Ap_id[i];
            for (i = 0; i < ndaughterAp; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    m_Ap_ptruth_t[i][j] = Ap_ptruth[i][j];
                }
            }

            m_ndaughterAm_t = ndaughterAm;
            for (i = 0; i < ndaughterAm; i++)
                m_Am_id_t[i] = Am_id[i];
            for (i = 0; i < ndaughterAm; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    m_Am_ptruth_t[i][j] = Am_ptruth[i][j];
                }
            }
            m_tuple2->write();
        }
    }

    //*************************Begin Good Charged Track Selection*****************************
    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!! begin good charged track selection !!!!!!!!!!" << std::endl;

    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);

    Vint iGood;
    iGood.clear();

    //	int nCharge = 0;
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

    for (i = 0; i < evtRecEvent->totalCharged(); i++)
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
        double costheta = cos(mdcTrk->theta());
        double Rvxy0 = fabs(vecipa[0]); //the nearest distance to IP in xy plane
        double Rvz0 = vecipa[3];        //the nearest distance to IP in z direction

        if (costheta >= m_costheta)
            continue;
        if (fabs(Rvz0) >= m_pvz0cut)
            continue;
        if (fabs(Rvxy0) >= m_pvr0cut)
            continue;
        iGood.push_back(i);
        //		nCharge += mdcTrk->charge();
    }
    //**************************Finish Good Charged Track Selection**************************

    int nGood = iGood.size();
    if (nGood < 1)
        return StatusCode::SUCCESS;
    Ncut0++;
    Ncharge += nGood;

    //**************************Find proton and anti_proton**************************
    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!!begin Proton and Anti-Proton Selection!!!!!!!!!!" << std::endl;

    Vint ipp;
    ipp.clear();
    Vint ipm;
    ipm.clear();

    for (i = 0; i < nGood; i++)
    {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
        RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
        ISimplePIDSvc *m_simplePIDSvc;
        Gaudi::svcLocator()->service("SimplePIDSvc", m_simplePIDSvc);
        m_simplePIDSvc->preparePID(*itTrk);
        if (m_simplePIDSvc->isproton())
        {
            if (mdcTrk->charge() == 1)
                ipp.push_back(iGood[i]);
            if (mdcTrk->charge() == -1)
                ipm.push_back(iGood[i]);
        }
    }
    //**************************Finish finding proton and anti_proton**************************

    int npp = ipp.size();
    int npm = ipm.size();
    if (npp == 0 || npm == 0) // suppress qqbar 
        return StatusCode::SUCCESS;
    Ncut1++;
    Nproton += npp;
    Nanti_proton += npm;

    //**************************obtain 4-momentum of proton and anti_proton**************************
    Vp4 p_pp, p_pm;
    p_pp.clear();
    p_pm.clear();

    for (int i = 0; i < npp; i++)
    {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ipp[i];
        RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
        mdcKalTrk->setPidType(RecMdcKalTrack::proton);
        HepLorentzVector p_pp4 = mdcKalTrk->p4(xmass[4]);
        p_pp.push_back(p_pp4);
    }
    for (int i = 0; i < npm; i++)
    {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + ipm[i];
        RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
        mdcKalTrk->setPidType(RecMdcKalTrack::proton);
        HepLorentzVector ppm4 = mdcKalTrk->p4(xmass[4]);
        p_pm.push_back(ppm4);
    }
    //**************************finish obtaining 4-momentum of proton and anti_proton**************************

    //*************************Begin to choose good shower for gamma*********************
    vector<Gam_info> iGam;
    iGam.clear();
    Gam_info tmp_gam;

    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!!begin choose good gamma !!!!!!!!!!" << std::endl;

    for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++)
    {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isEmcShowerValid())
            continue;
        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
        HepLorentzVector shP4 = getP4(emcTrk, xorigin);
        double cosThetaSh = shP4.vect().cosTheta();

        //cut angle between anti-proton and gamma , do nothing with gamma and other charged track
        double dang = 200.;
        for (int j = 0; j < evtRecEvent->totalCharged(); j++)
        {
            EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
            if (!(*jtTrk)->isExtTrackValid())
                continue;
            RecExtTrack *extTrk = (*jtTrk)->extTrack();
            if (extTrk->emcVolumeNumber() == -1)
                continue;

            // RecMdcTrack *mdcTrack = (*jtTrk)->mdcTrack();
            // int charge = mdcTrack->charge();

            // ISimplePIDSvc *m_simplePIDSvc;
            // Gaudi::svcLocator()->service("SimplePIDSvc", m_simplePIDSvc);
            // m_simplePIDSvc->preparePID(*jtTrk);

            // if (charge > 0 || !(m_simplePIDSvc->isproton()))
                // continue;

            Hep3Vector extpos = extTrk->emcPosition();
            double angd = extpos.angle(emcpos);
            if (fabs(angd) < dang)
                dang = fabs(angd);
        }

        double eraw = emcTrk->energy();
        double getTime = emcTrk->time();
        dang = dang * 180 / (CLHEP::pi);

        if (getTime > m_gammathCut || getTime < m_gammatlCut)
            continue;
        if (!((fabs(cosThetaSh) < m_maxCosThetaBarrel && eraw > m_minEnergy) || ((fabs(cosThetaSh) > m_minCosThetaEndcap) && (fabs(cosThetaSh) < m_maxCosThetaEndcap) && (eraw > m_minEndcapEnergy))))
            continue;

        tmp_gam.index = i;
        tmp_gam.ang_antiproton = dang;
        iGam.push_back(tmp_gam);
    }
    //********************************* Finish Good Photon Selection*****************************

    int nGam = iGam.size();
    if (nGam < 4 || nGam > 15)
        return StatusCode::SUCCESS;
    Ncut2++;
    Ngamma += nGam;

    //*******************  begin 3C fit for 2-pi0 & 1-eta and 4C for 2-pi0 & 1-eta & 1-sigma+ **********************
    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!! begin 3+4C fit  !!!!!!!!!!" << std::endl;
    //there are seven 'for' cycle.
    bool writeflag = false;
    group tmp_groupp;
    vector<group> lcp;
    lcp.clear();

    //j,k cycle is for pi0 in sigma
    for (int j = 0; j < nGam - 1; j++)
    {
        EvtRecTrackIterator itTrkj = evtRecTrkCol->begin() + iGam[j].index;
        RecEmcShower *emcTrkj = (*itTrkj)->emcShower();
        Hep3Vector emcposj(emcTrkj->x(), emcTrkj->y(), emcTrkj->z());
        HepLorentzVector ptrkj = getP4(emcTrkj, xorigin);

        for (int k = j + 1; k < nGam; k++)
        {
            EvtRecTrackIterator itTrkk = evtRecTrkCol->begin() + iGam[k].index;
            RecEmcShower *emcTrkk = (*itTrkk)->emcShower();
            Hep3Vector emcposk(emcTrkk->x(), emcTrkk->y(), emcTrkk->z());
            HepLorentzVector ptrkk = getP4(emcTrkk, xorigin);

            //cut pi0 mass
            HepLorentzVector p4_pi0_sigma = ptrkj + ptrkk;
            double Mass_pi0_sigma = p4_pi0_sigma.m();
            if ((Mass_pi0_sigma <= m_Pi0minMass) || (Mass_pi0_sigma >= m_Pi0maxMass))
                continue;

            m_pi0 = Mass_pi0_sigma > 0 ? Mass_pi0_sigma : -10;

            // n, o for eta
            for (int n = 0; n < nGam - 1; n++)
            {
                EvtRecTrackIterator itTrkn = evtRecTrkCol->begin() + iGam[n].index;
                RecEmcShower *emcTrkn = (*itTrkn)->emcShower();
                Hep3Vector emcposn(emcTrkn->x(), emcTrkn->y(), emcTrkn->z());
                HepLorentzVector ptrkn = getP4(emcTrkn, xorigin);

                for (int o = n + 1; o < nGam; o++)
                {
                    EvtRecTrackIterator itTrko = evtRecTrkCol->begin() + iGam[o].index;
                    RecEmcShower *emcTrko = (*itTrko)->emcShower();
                    Hep3Vector emcposo(emcTrko->x(), emcTrko->y(), emcTrko->z());
                    HepLorentzVector ptrko = getP4(emcTrko, xorigin);
                    //cut eta mass
                    HepLorentzVector p4_eta = ptrkn + ptrko;
                    double Mass_eta = p4_eta.m();
                    if ((Mass_eta <= m_EtaminMass) || (Mass_eta >= m_EtamaxMass))
                        continue;

                    //prevent one gamma being used for more than once.
                    Vint allgamma;
                    allgamma.clear();
                    allgamma.push_back(iGam[j].index);
                    allgamma.push_back(iGam[k].index);
                    allgamma.push_back(iGam[n].index);
                    allgamma.push_back(iGam[o].index);
                    if (!gammaidentifier(allgamma))
                        continue;

                    m_eta = Mass_eta > 0 ? Mass_eta : -10;
                    //fit
                    KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
                    kmfit->init();
                    kmfit->setChisqCut(m_chisqMax);
                    kmfit->setIterNumber(4);
                    kmfit->AddTrack(0, 0.0, emcTrkj);
                    kmfit->AddTrack(1, 0.0, emcTrkk);
                    kmfit->AddTrack(2, 0.0, emcTrkn);
                    kmfit->AddTrack(3, 0.0, emcTrko);

                    kmfit->AddResonance(0, m_Pi0MassFit, 0, 1);
                    kmfit->AddResonance(1, m_EtaMassFit, 2, 3);

                    bool okvs = kmfit->Fit();

                    if (!okvs)
                        continue;
                    else
                    {
                        double chisq_3C = kmfit->chisq();
                        HepLorentzVector gamma1_sigma_3C, gamma2_sigma_3C;
                        gamma1_sigma_3C = kmfit->pfit(0);
                        gamma2_sigma_3C = kmfit->pfit(1);
                        //i cycle is for sigma
                        for (int i = 0; i < npp; i++)
                        {
                            EvtRecTrackIterator itTrkpp = evtRecTrkCol->begin() + ipp[i];
                            RecMdcKalTrack *ppTrk = (*itTrkpp)->mdcKalTrack();
                            RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);
                            WTrackParameter wpptrk(xmass[4], ppTrk->getZHelixP(), ppTrk->getZErrorP());

                            //cut sigma mass
                            HepLorentzVector p4_Sigmap = p_pp[i] + gamma1_sigma_3C + gamma2_sigma_3C;
                            double Mass_sigma = p4_Sigmap.m();
                            if ((Mass_sigma <= m_SigmaminMass) || (Mass_sigma >= m_SigmamaxMass))
                                continue;

                            m_sigmap = Mass_sigma > 0 ? Mass_sigma : -10;
                            //fit
                            KalmanKinematicFit *kmfit2 = KalmanKinematicFit::instance();
                            kmfit2->init();
                            kmfit2->setChisqCut(m_chisqMax);
                            kmfit2->setIterNumber(5);
                            kmfit2->AddTrack(0, wpptrk);
                            kmfit2->AddTrack(1, 0.0, emcTrkj);
                            kmfit2->AddTrack(2, 0.0, emcTrkk);
                            kmfit2->AddTrack(3, 0.0, emcTrkn);
                            kmfit2->AddTrack(4, 0.0, emcTrko);

                            kmfit2->AddResonance(0, m_SigmaMassFit, 0, 1, 2);
                            kmfit2->AddResonance(1, m_Pi0MassFit, 1, 2);
                            kmfit2->AddResonance(2, m_EtaMassFit, 3, 4);

                            bool okvs2 = kmfit2->Fit();

                            if (!okvs2)
                                continue;
                            else
                            {

                                //save chisq
                                tmp_groupp.chisquare_3C = chisq_3C;
                                tmp_groupp.chisquare_4C = kmfit2->chisq();
                                //save p4 after 3C
                                tmp_groupp.p4_gam1_pi0_sigma_3C = gamma1_sigma_3C;
                                tmp_groupp.p4_gam2_pi0_sigma_3C = gamma2_sigma_3C;
                                //save p4 after 4C
                                tmp_groupp.p4_proton_sigma_4C = kmfit2->pfit(0);

                                tmp_groupp.p4_gam1_pi0_sigma_4C = kmfit2->pfit(1);
                                tmp_groupp.p4_gam2_pi0_sigma_4C = kmfit2->pfit(2);
                                tmp_groupp.p4_gam1_eta_4C = kmfit2->pfit(3);
                                tmp_groupp.p4_gam2_eta_4C = kmfit2->pfit(4);
                                //save p4 old
                                tmp_groupp.p4_proton_sigma_old = p_pp[i];
                                tmp_groupp.p4_gam1_pi0_sigma_old = ptrkj;
                                tmp_groupp.p4_gam2_pi0_sigma_old = ptrkk;
                                // tmp_groupp.p4_gam1_pi0_old = ptrkl;
                                // tmp_groupp.p4_gam2_pi0_old = ptrkm;
                                tmp_groupp.p4_gam1_eta_old = ptrkn;
                                tmp_groupp.p4_gam2_eta_old = ptrko;
                                //save minimum angle
                                Vdouble allangle;
                                allangle.clear();
                                allangle.push_back(iGam[j].ang_antiproton);
                                allangle.push_back(iGam[k].ang_antiproton);
                                allangle.push_back(iGam[n].ang_antiproton);
                                allangle.push_back(iGam[o].ang_antiproton);
                                double tmp_angle = 100.;

                                for (int b = 0; b < 4; b++)
                                {
                                    if (allangle[b] < tmp_angle)
                                        tmp_angle = allangle[b];
                                };
                                tmp_groupp.angle_gamma_antiproton = tmp_angle;

                                lcp.push_back(tmp_groupp);
                                writeflag = true;
                            }
                        }
                    }
                }
            }
            //     }
            // }
        }
    }
    //*******************  finish loop sigma+  eta  **********************

    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!! finish 3+4C fit for lambdac+  !!!!!!!!!!" << std::endl;

    //*******************  begin loop anti-sigma- pi0  **********************
    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!! begin 3+4C fit for lambdac-  !!!!!!!!!!" << std::endl;
    bool mwriteflag = false;
    group tmp_groupm;
    vector<group> lcm;
    lcm.clear();

    //j,k cycle is for pi0 in sigma
    for (int j = 0; j < nGam - 1; j++)
    {
        EvtRecTrackIterator itTrkj = evtRecTrkCol->begin() + iGam[j].index;
        RecEmcShower *emcTrkj = (*itTrkj)->emcShower();
        Hep3Vector emcposj(emcTrkj->x(), emcTrkj->y(), emcTrkj->z());
        HepLorentzVector ptrkj = getP4(emcTrkj, xorigin);

        for (int k = j + 1; k < nGam; k++)
        {
            EvtRecTrackIterator itTrkk = evtRecTrkCol->begin() + iGam[k].index;
            RecEmcShower *emcTrkk = (*itTrkk)->emcShower();
            Hep3Vector emcposk(emcTrkk->x(), emcTrkk->y(), emcTrkk->z());
            HepLorentzVector ptrkk = getP4(emcTrkk, xorigin);
            //cut pi0 mass
            HepLorentzVector p4_pi0_sigma = ptrkj + ptrkk;
            double Mass_pi0_sigma = p4_pi0_sigma.m();
            if ((Mass_pi0_sigma <= m_Pi0minMass) || (Mass_pi0_sigma >= m_Pi0maxMass))
                continue;

            m_pi0 = Mass_pi0_sigma > 0 ? Mass_pi0_sigma : -10;

            //n,o cycle is for eta

            for (int n = 0; n < nGam - 1; n++)
            {
                EvtRecTrackIterator itTrkn = evtRecTrkCol->begin() + iGam[n].index;
                RecEmcShower *emcTrkn = (*itTrkn)->emcShower();
                Hep3Vector emcposn(emcTrkn->x(), emcTrkn->y(), emcTrkn->z());
                HepLorentzVector ptrkn = getP4(emcTrkn, xorigin);

                for (int o = n + 1; o < nGam; o++)
                {
                    EvtRecTrackIterator itTrko = evtRecTrkCol->begin() + iGam[o].index;
                    RecEmcShower *emcTrko = (*itTrko)->emcShower();
                    Hep3Vector emcposo(emcTrko->x(), emcTrko->y(), emcTrko->z());
                    HepLorentzVector ptrko = getP4(emcTrko, xorigin);
                    //cut eta mass
                    HepLorentzVector p4_eta = ptrkn + ptrko;
                    double Mass_eta = p4_eta.m();
                    if ((Mass_eta <= m_EtaminMass) || (Mass_eta >= m_EtamaxMass))
                        continue;

                    //prevent one gamma being used for more than once.
                    Vint allgamma;
                    allgamma.clear();
                    allgamma.push_back(iGam[j].index);
                    allgamma.push_back(iGam[k].index);
                    allgamma.push_back(iGam[n].index);
                    allgamma.push_back(iGam[o].index);
                    if (!gammaidentifier(allgamma))
                        continue;

                    m_eta = Mass_eta > 0 ? Mass_eta : -10;

                    //fit
                    KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
                    kmfit->init();
                    kmfit->setChisqCut(m_chisqMax);
                    kmfit->setIterNumber(5);
                    kmfit->AddTrack(0, 0.0, emcTrkj);
                    kmfit->AddTrack(1, 0.0, emcTrkk);
                    kmfit->AddTrack(2, 0.0, emcTrkn);
                    kmfit->AddTrack(3, 0.0, emcTrko);

                    kmfit->AddResonance(0, m_Pi0MassFit, 0, 1);
                    kmfit->AddResonance(1, m_EtaMassFit, 2, 3);
                    bool okvs = kmfit->Fit();

                    if (!okvs)
                        continue;
                    else
                    {
                        double chisq_3C = kmfit->chisq();
                        HepLorentzVector gamma1_sigma_3C, gamma2_sigma_3C;
                        gamma1_sigma_3C = kmfit->pfit(0);
                        gamma2_sigma_3C = kmfit->pfit(1);
                        //i cycle is for sigma
                        for (int i = 0; i < npm; i++)
                        {
                            EvtRecTrackIterator itTrkpm = evtRecTrkCol->begin() + ipm[i];
                            RecMdcKalTrack *pmTrk = (*itTrkpm)->mdcKalTrack();
                            RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);
                            WTrackParameter wpmtrk(xmass[4], pmTrk->getZHelixP(), pmTrk->getZErrorP());

                            //cut sigma mass
                            HepLorentzVector p4_Sigmam = p_pm[i] + gamma1_sigma_3C + gamma2_sigma_3C;
                            double Mass_sigma = p4_Sigmam.m();
                            if ((Mass_sigma <= m_SigmaminMass) || (Mass_sigma >= m_SigmamaxMass))
                                continue;

                            m_sigmap = Mass_sigma > 0 ? Mass_sigma : -10;
                            //fit
                            KalmanKinematicFit *kmfit2 = KalmanKinematicFit::instance();
                            kmfit2->init();
                            kmfit2->setChisqCut(m_chisqMax);
                            kmfit2->setIterNumber(5);
                            kmfit2->AddTrack(0, wpmtrk);
                            kmfit2->AddTrack(1, 0.0, emcTrkj);
                            kmfit2->AddTrack(2, 0.0, emcTrkk);
                            kmfit2->AddTrack(3, 0.0, emcTrkn);
                            kmfit2->AddTrack(4, 0.0, emcTrko);

                            kmfit2->AddResonance(0, m_SigmaMassFit, 0, 1, 2);
                            kmfit2->AddResonance(1, m_Pi0MassFit, 1, 2);
                            kmfit2->AddResonance(2, m_EtaMassFit, 3, 4);

                            bool okvs2 = kmfit2->Fit();
                            if (!okvs2)
                                continue;
                            else
                            {
                                //save chisq
                                tmp_groupm.chisquare_3C = chisq_3C;
                                tmp_groupm.chisquare_4C = kmfit2->chisq();
                                //save p4 after 3C
                                tmp_groupm.p4_gam1_pi0_sigma_3C = gamma1_sigma_3C;
                                tmp_groupm.p4_gam2_pi0_sigma_3C = gamma2_sigma_3C;
                                //save p4 after 4C
                                tmp_groupm.p4_proton_sigma_4C = kmfit2->pfit(0);
                                ;
                                tmp_groupm.p4_gam1_pi0_sigma_4C = kmfit2->pfit(1);
                                tmp_groupm.p4_gam2_pi0_sigma_4C = kmfit2->pfit(2);
                                tmp_groupm.p4_gam1_eta_4C = kmfit2->pfit(3);
                                tmp_groupm.p4_gam2_eta_4C = kmfit2->pfit(4);
                                //save p4 old
                                tmp_groupm.p4_proton_sigma_old = p_pm[i];
                                tmp_groupm.p4_gam1_pi0_sigma_old = ptrkj;
                                tmp_groupm.p4_gam2_pi0_sigma_old = ptrkk;
                                tmp_groupm.p4_gam1_eta_old = ptrkn;
                                tmp_groupm.p4_gam2_eta_old = ptrko;
                                //save minimum angle
                                Vdouble allangle;
                                allangle.clear();
                                allangle.push_back(iGam[j].ang_antiproton);
                                allangle.push_back(iGam[k].ang_antiproton);
                                allangle.push_back(iGam[n].ang_antiproton);
                                allangle.push_back(iGam[o].ang_antiproton);
                                double tmp_angle = 100.;
                                for (int b = 0; b < 4; b++)
                                {
                                    if (allangle[b] < tmp_angle)
                                        tmp_angle = allangle[b];
                                };
                                tmp_groupm.angle_gamma_antiproton = tmp_angle;

                                lcm.push_back(tmp_groupm);
                                mwriteflag = true;
                            }
                        }
                    }
                }
            }
        }
    }
    //*******************  finish loop anti-sigma- pi0 eta  **********************

    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!! finish 3C + 4C fit for lambdac-  !!!!!!!!!!" << std::endl;

    int isig_Lcp = 0;
    isig_Lcp = lcp.size();
    if (isig_Lcp > 2000)
        writeflag = 0;

    int isig_Lcm = 0;
    isig_Lcm = lcm.size();
    if (isig_Lcm > 2000)
        mwriteflag = 0;

    if (m_debug)
        std::cerr << __LINE__ << "isig_Lcp = " << isig_Lcp << " isig_Lcm = " << isig_Lcm << std::endl;
    if (m_debug)
        std::cerr << __LINE__ << "writeflag=" << writeflag << " , mwriteflag=" << mwriteflag << std::endl;
    if ((!writeflag) && (!mwriteflag))
        return StatusCode::SUCCESS;
    Ncut3++;
    if (writeflag)
        Ncut4++;
    if (mwriteflag)
        Ncut5++;

    ///////////////////////
    //  // get beam energy and beta
    // ///////////////////////

    /*	if(m_ReadBeamEFromDB)
    {
        if(m_usecalibBeamE) m_readDb.setcalib(true);
        m_beamE=m_readDb.getbeamE(m_run,m_beamE);
        if(m_run>0) m_beta=m_readDb.getbeta();
    }
*/

    //******************* begin for Lambdac+ Decay****************************
    if (m_debug)
    std:
        cerr << __LINE__ << "!!!!!!!!!!m_debug Begin loop Lambdac+ !!!!!!!!!!" << endl;

    double ebeam = m_beamE;

    if (writeflag)
    {
        m_run = runNo;
        m_event = eventNo;
        m_rightflag = 1;

        if (m_getmcinfo)
        {
            m_mode1 = mode1;
            m_mode2 = mode2;
            m_mode3 = mode3;

            m_Lmdc_P = Lmdc_P;
            m_Lmdc_M = Lmdc_M;

            m_idxmc = numParticle;
            for (int i = 0; i < numParticle; i++)
            {
                m_pdgid[i] = M_pdgid[i];
                m_motheridx[i] = M_motheridx[i];
            }

            m_mcparticle_p = numParticle_p;
            for (int i = 0; i < numParticle_p; i++)
            {
                m_pdgid_p[i] = M_pdgid_p[i];
                m_motheridx_p[i] = M_motheridx_p[i];
            }

            m_mcparticle_m = numParticle_m;
            for (int i = 0; i < numParticle_m; i++)
            {
                m_pdgid_m[i] = M_pdgid_m[i];
                m_motheridx_m[i] = M_motheridx_m[i];
            }

            m_ndaughterAp = ndaughterAp;
            for (i = 0; i < ndaughterAp; i++)
                m_Ap_id[i] = Ap_id[i];
            for (i = 0; i < ndaughterAp; i++)
                for (j = 0; j < 4; j++)
                    m_Ap_ptruth[i][j] = Ap_ptruth[i][j];
        }

        m_isig = isig_Lcp;
        for (int isig = 0; isig < isig_Lcp; isig++)
        {
            //pi0
            for (i = 0; i < 4; i++)
                m_Gam1_Eta_p4_old[isig][i] = lcp[isig].p4_gam1_eta_old[i];
            for (i = 0; i < 4; i++)
                m_Gam2_Eta_p4_old[isig][i] = lcp[isig].p4_gam2_eta_old[i];
            for (i = 0; i < 4; i++)
                m_Gam1_Eta_p4[isig][i] = lcp[isig].p4_gam1_eta_4C[i];
            for (i = 0; i < 4; i++)
                m_Gam2_Eta_p4[isig][i] = lcp[isig].p4_gam2_eta_4C[i];
            //sigma
            for (i = 0; i < 4; i++)
                m_gam1_pi0_Sigma_p4[isig][i] = lcp[isig].p4_gam1_pi0_sigma_4C[i];
            for (i = 0; i < 4; i++)
                m_gam2_pi0_Sigma_p4[isig][i] = lcp[isig].p4_gam2_pi0_sigma_4C[i];
            for (i = 0; i < 4; i++)
                m_p_Sigma_p4[isig][i] = lcp[isig].p4_proton_sigma_4C[i];

            for (i = 0; i < 4; i++)
                m_gam1_pi0_Sigma_p4_old[isig][i] = lcp[isig].p4_gam1_pi0_sigma_old[i];
            for (i = 0; i < 4; i++)
                m_gam2_pi0_Sigma_p4_old[isig][i] = lcp[isig].p4_gam2_pi0_sigma_old[i];
            for (i = 0; i < 4; i++)
                m_p_Sigma_p4_old[isig][i] = lcp[isig].p4_proton_sigma_old[i];

            for (i = 0; i < 4; i++)
                m_gam1_pi0_Sigma_p4_3C[isig][i] = lcp[isig].p4_gam1_pi0_sigma_3C[i];
            for (i = 0; i < 4; i++)
                m_gam2_pi0_Sigma_p4_3C[isig][i] = lcp[isig].p4_gam2_pi0_sigma_3C[i];

            //mass,chisq,angle,mbc,deltaE
            m_charge[isig] = 1;
            m_chisq_3C[isig] = lcp[isig].chisquare_3C;
            m_chisq_4C[isig] = lcp[isig].chisquare_4C;
            m_angle[isig] = lcp[isig].angle_gamma_antiproton;

            HepLorentzVector p4_sigma_4C = lcp[isig].p4_gam1_pi0_sigma_4C + lcp[isig].p4_gam2_pi0_sigma_4C + lcp[isig].p4_proton_sigma_4C;
            // HepLorentzVector p4_pi0_4C = lcp[isig].p4_gam1_pi0_4C + lcp[isig].p4_gam2_pi0_4C;
            HepLorentzVector p4_eta_4C = lcp[isig].p4_gam1_eta_4C + lcp[isig].p4_gam2_eta_4C;
            // HepLorentzVector pLambdac = p4_sigma_4C + p4_pi0_4C + p4_eta_4C;
            HepLorentzVector pLambdac = p4_sigma_4C + p4_eta_4C;
            pLambdac.boost(-m_beta);
            m_deltaE_min[isig] = pLambdac.t() - ebeam;

            double mbc2 = ebeam * ebeam - pLambdac.v().mag2();
            m_bc[isig] = mbc2 > 0 ? sqrt(mbc2) : -10;

            // m_sigmap = p4_sigma_4C.m() > 0 ? p4_sigma_4C.m() : -10;
            // m_eta = p4_eta_4C.m() > 0 ? p4_eta_4C.m() : -10;
            // HepLorentzVector p_pi0 = lcp[isig].p4_gam1_pi0_sigma_4C + lcp[isig].p4_gam2_pi0_sigma_4C;
            // m_pi0 = p_pi0.m() > 0 ? p_pi0.m() : -10;
            //after give value for all branch
        }

        m_tuple1->write();
    };
    //******************* End Lambdac+ Decay****************************

    if (m_debug)
        std::cerr << __LINE__ << "******************** finish write lambdac+ ******************" << std::endl;

    //******************* begin for Anti-Lambdac- Decay****************************
    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!!m_debug Begin loop anti-Lambdac- !!!!!!!!!!" << endl;

    if (mwriteflag)
    {
        m_run = runNo;
        m_event = eventNo;
        m_rightflag = 1;

        if (m_getmcinfo)
        {
            m_mode1 = mode1;
            m_mode2 = mode2;
            m_mode3 = mode3;

            m_Lmdc_P = Lmdc_P;
            m_Lmdc_M = Lmdc_M;

            m_idxmc = numParticle;
            for (int i = 0; i < numParticle; i++)
            {
                m_pdgid[i] = M_pdgid[i];
                m_motheridx[i] = M_motheridx[i];
            }

            m_mcparticle_p = numParticle_p;
            for (int i = 0; i < numParticle_p; i++)
            {
                m_pdgid_p[i] = M_pdgid_p[i];
                m_motheridx_p[i] = M_motheridx_p[i];
            }

            m_mcparticle_m = numParticle_m;
            for (int i = 0; i < numParticle_m; i++)
            {
                m_pdgid_m[i] = M_pdgid_m[i];
                m_motheridx_m[i] = M_motheridx_m[i];
            }

            m_ndaughterAm = ndaughterAm;
            for (i = 0; i < ndaughterAm; i++)
                m_Am_id[i] = Am_id[i];
            for (i = 0; i < ndaughterAm; i++)
                for (j = 0; j < 4; j++)
                    m_Am_ptruth[i][j] = Am_ptruth[i][j];
        }

        m_isig = isig_Lcm;
        for (int isig = 0; isig < isig_Lcm; isig++)
        {
            //p4
            //pi0
            for (i = 0; i < 4; i++)
                m_Gam1_Eta_p4_old[isig][i] = lcm[isig].p4_gam1_eta_old[i];
            for (i = 0; i < 4; i++)
                m_Gam2_Eta_p4_old[isig][i] = lcm[isig].p4_gam2_eta_old[i];
            for (i = 0; i < 4; i++)
                m_Gam1_Eta_p4[isig][i] = lcm[isig].p4_gam1_eta_4C[i];
            for (i = 0; i < 4; i++)
                m_Gam2_Eta_p4[isig][i] = lcm[isig].p4_gam2_eta_4C[i];
            //anti-sigma
            for (i = 0; i < 4; i++)
                m_gam1_pi0_Sigma_p4[isig][i] = lcm[isig].p4_gam1_pi0_sigma_4C[i];
            for (i = 0; i < 4; i++)
                m_gam2_pi0_Sigma_p4[isig][i] = lcm[isig].p4_gam2_pi0_sigma_4C[i];
            for (i = 0; i < 4; i++)
                m_p_Sigma_p4[isig][i] = lcm[isig].p4_proton_sigma_4C[i];

            for (i = 0; i < 4; i++)
                m_gam1_pi0_Sigma_p4_old[isig][i] = lcm[isig].p4_gam1_pi0_sigma_old[i];
            for (i = 0; i < 4; i++)
                m_gam2_pi0_Sigma_p4_old[isig][i] = lcm[isig].p4_gam2_pi0_sigma_old[i];
            for (i = 0; i < 4; i++)
                m_p_Sigma_p4_old[isig][i] = lcm[isig].p4_proton_sigma_old[i];

            for (i = 0; i < 4; i++)
                m_gam1_pi0_Sigma_p4_3C[isig][i] = lcm[isig].p4_gam1_pi0_sigma_3C[i];
            for (i = 0; i < 4; i++)
                m_gam2_pi0_Sigma_p4_3C[isig][i] = lcm[isig].p4_gam2_pi0_sigma_3C[i];

            //mass,chisq,angle,mbc,deltaE
            m_charge[isig] = -1;
            m_chisq_3C[isig] = lcm[isig].chisquare_3C;
            m_chisq_4C[isig] = lcm[isig].chisquare_4C;
            m_angle[isig] = lcm[isig].angle_gamma_antiproton;

            HepLorentzVector p4_sigma_4C = lcm[isig].p4_gam1_pi0_sigma_4C + lcm[isig].p4_gam2_pi0_sigma_4C + lcm[isig].p4_proton_sigma_4C;
            // HepLorentzVector p4_pi0_4C = lcm[isig].p4_gam1_pi0_4C + lcm[isig].p4_gam2_pi0_4C;
            HepLorentzVector p4_eta_4C = lcm[isig].p4_gam1_eta_4C + lcm[isig].p4_gam2_eta_4C;
            HepLorentzVector pLambdac = p4_sigma_4C + p4_eta_4C;
            pLambdac.boost(-m_beta);
            m_deltaE_min[isig] = pLambdac.t() - ebeam;

            double mbc2 = ebeam * ebeam - pLambdac.v().mag2();
            m_bc[isig] = mbc2 > 0 ? sqrt(mbc2) : -10;

            // m_sigmap = p4_sigma_4C.m() > 0 ? p4_sigma_4C.m() : -10;
            // m_eta = p4_eta_4C.m() > 0 ? p4_eta_4C.m() : -10;
            // HepLorentzVector p_pi0 = lcm[isig].p4_gam1_pi0_sigma_4C + lcm[isig].p4_gam2_pi0_sigma_4C;
            // m_pi0 = p_pi0.m() > 0 ? p_pi0.m() : -10;

            //after give value for all branch
        }

        m_tuple1->write();
    }
    //******************* End Lambdac- Decay****************************

    if (m_debug)
        std::cerr << __LINE__ << "******************** finish write lambdac- ******************" << std::endl;

    return StatusCode::SUCCESS;
}

StatusCode LambdacSigmapEta::endRun()
{
    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!! m_debug endRun !!!!!!!!!!" << std::endl;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "LambdacSigmapEta::endRun()" << endreq;
    //add your code here
    return StatusCode::SUCCESS;
}

StatusCode LambdacSigmapEta::finalize()
{
    if (m_debug)
        std::cerr << __LINE__ << "!!!!!!!!!! m_debug begin finalize !!!!!!!!!!" << std::endl;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "LambdacSigmapEta::finalize()" << endreq;
    //add your code here
    std::cerr << "Totalevent   " << Ntotal << std::endl;
    std::cerr << "N_evnet_include_charged   " << Ncut0 << std::endl;
    std::cerr << "N_event_include_p   " << Ncut1 << std::endl;
    std::cerr << "N_event_include_gamma   " << Ncut2 << std::endl;
    std::cerr << "N_event_include_Sigma+/-_Eta   " << Ncut3 << std::endl;
    std::cerr << "N_event_include_Sigma+_Eta   " << Ncut4 << std::endl;
    std::cerr << "N_event_include_anti_Sigma-_Eta   " << Ncut5 << std::endl;
    std::cerr << "N_charged   " << Ncharge << std::endl;
    std::cerr << "N_proton   " << Nproton << std::endl;
    std::cerr << "N_anti_proton   " << Nanti_proton << std::endl;
    std::cerr << "N_gamma   " << Ngamma << std::endl;

    return StatusCode::SUCCESS;
}

//add your code here,for other member-functions
