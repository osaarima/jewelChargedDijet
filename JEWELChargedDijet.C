#include "TH1D.h"
#include "TProfile.h"
#include "TMath.h"
#include "TFile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <TStopwatch.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepPDT/defs.h>
#include <HepPDT/TableBuilder.hh>
#include <HepPDT/TempParticleData.hh>
#include <HepPDT/ParticleDataTable.hh>

#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
/*
#include "fastjet/config.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include <fastjet/SISConePlugin.hh>
#include <fastjet/CDFMidPointPlugin.hh>
#ifdef FASTJET_VERSION
#include <fastjet/Selector.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <fastjet/tools/Subtractor.hh>
*/
#include <iostream>

#include <TClonesArray.h>
#include "src/AliJCDijetHistos.h"
#include "src/AliJBaseTrack.h"
//#include "src/JHistos.h"
//#include "src/AliJCard.h"

#define DEBUG 0

using namespace fastjet;
using namespace std;
void CalculateJetsDijets(TClonesArray *inList,
                                         int    lDebug,
                                         int    lCBin,
                                         double lParticleEtaCut,
                                         double lParticlePtCut,
                                         double lJetCone,
                                         double lktJetCone,
                                         int    lktScheme,
                                         bool   lusePionMassInkt,
                                         bool   luseDeltaPhiBGSubtr,
                                         double lConstituentCut,
                                         double lLeadingJetCut,
                                         double lSubleadingJetCut,
                                         double lDeltaPhiCut);


class AliJCDijetHistos;

AliJCDijetHistos *fhistos;


int main(int argc, char **argv) {

	if(argc<5){
		cout<<"usage: " << argv[0] << " <input.hepmc> <output.root> <particle.tbl> dijetLeadingPt"<<endl;exit(1);
	}
	TStopwatch timer; 
	timer.Start();   

	char* hepmcFilename  = argv[1];
	TString outputs = argv[2];
	TString sParticle = argv[3];
    double dijetLeadingPt = atof(argv[4]);

	TFile *fout = new TFile(outputs.Data(),"RECREATE");
	fout->cd();//opening of the output file
    TDirectoryFile *fdir = new TDirectoryFile( "JCDijetBaseTask","JCDijetBaseTask" );
    fdir->cd();

	HepMC::IO_GenEvent ascii_in(hepmcFilename,std::ios::in);

    // Initialize HepPDF so that particle ID's can be used to determine charge and hadrons
    std::ifstream pdfile( sParticle.Data() );
    if( !pdfile ) { 
      std::cerr << "cannot open " << sParticle.Data() << std::endl;
      exit(-1);
    }
    HepPDT::ParticleDataTable datacol( "Generic Particle Table" );
    {
        // Construct table builder
        HepPDT::TableBuilder  tb(datacol);
        // bool  addParticleTable( std::istream&, TableBuilder&,
        //                         bool validate = false );
        // where:  validate=true => verify that the ParticleID is valid
        if( !addParticleTable( pdfile, tb, true ) ) { 
            std::cout << "error reading PDG pdt file " << std::endl; 
        }
    } // the tb destructor fills datacol

    //-------------------------------------------------------
    // Histograms and tools
    //-------------------------------------------------------
    fhistos = new AliJCDijetHistos();
    vector<double> centbins = {0.0, 100.0};
    fhistos->SetCentralityBinsHistos(centbins);
    fhistos->CreateEventTrackHistos();

    fhistos->fHMG->Print();

    TH1D *hCrossSectionInfo = new TH1D("hCrossSection","CrossSectionInfo",8,0,8);

    //------------------------------------------------------------------
    // Define jet reconstruction
    //------------------------------------------------------------------
    TClonesArray *inputList = new TClonesArray("AliJBaseTrack",1500);

    double partMinPtCut         = 0.15;// atlas 0.5 cms/alice 0.15
    double partMinEtaCut        = 0.8;
    double coneR                = 0.4; // atlas 0.6, cms 0.7 alice 0.4
	double ktconeR              = 0.4;
	double fusePionMassInktjets = false;
	double fuseDeltaPhiBGSubtr  = false;
    double jetConstituentCut    = 5.0;
    double dijetSubleadingPt    = 20.0;
    double dijetDeltaPhiCut     = 2.0; // Cut is pi/dijetDeltaPhiCut
    int centBin=0;
    int fktScheme               = 1;

    TString sktScheme;
    switch (fktScheme) {
        case 0:  sktScheme = "E_scheme";
                 break;
        case 1:  sktScheme = "pt_scheme";
                 break;
        case 2:  sktScheme = "pt2_scheme";
                 break;
        case 3:  sktScheme = "Et_scheme";
                 break;
        case 4:  sktScheme = "Et2_scheme";
                 break;
        case 5:  sktScheme = "BIpt_scheme";
                 break;
        case 6:  sktScheme = "BIpt2_scheme";
                 break;
        default: sktScheme = "Unknown, check macro arguments!";
                 break;
    }

    cout << endl;
    cout << "============= Settings =============" << endl;
    cout << "cent bin:                   " << centBin << endl;
    cout << "particle eta cut:           " << partMinEtaCut << endl;
    cout << "particle pt cut:            " << Form("%.2f",partMinPtCut) << endl;
    cout << "jet cone size:              " << coneR << endl;
    cout << "kt-jet cone size:           " << ktconeR << endl;
    cout << "Using pion mass in kt-jets: " << fusePionMassInktjets << endl;
    cout << "Using DeltaPhi in BG subtr: " << fuseDeltaPhiBGSubtr << endl;
    cout << "Using kt-jet scheme:        " << sktScheme.Data() << endl;
    cout << "jet costituent cut:         " << jetConstituentCut << endl;
    cout << "dijet leading pt cut:       " << dijetLeadingPt << endl;
    cout << "dijet subleading pt cut:    " << dijetSubleadingPt << endl;
    cout << "dijet DeltaPhi cut:         pi/" << dijetDeltaPhiCut << endl;
    cout << endl;

    if(fusePionMassInktjets && fktScheme!=0) {
        cout << "Warning: Using pion mass for kt-jets but not using E_scheme!" << endl;
        cout << endl;
    }

    // Save information about the settings used.
    fhistos->fh_info->Fill("Count", 1.0);
    fhistos->fh_info->Fill("MC", 1.0);
    fhistos->fh_info->Fill("Cent bin border 00", 0.0);
    fhistos->fh_info->Fill("Cent bin border 01", 5.0);
    fhistos->fh_info->Fill("Jet cone", coneR);
    fhistos->fh_info->Fill("kt-jet cone", ktconeR);
    fhistos->fh_info->Fill("Scheme", fktScheme);
    fhistos->fh_info->Fill("Use pion mass", fusePionMassInktjets);
    fhistos->fh_info->Fill("Particle eta cut", partMinEtaCut);
    fhistos->fh_info->Fill("Particle pt cut", partMinPtCut);
    fhistos->fh_info->Fill("Leading jet cut", dijetLeadingPt);
    fhistos->fh_info->Fill("Subleading jet cut", dijetSubleadingPt);
    fhistos->fh_info->Fill("Const. cut", jetConstituentCut);
    fhistos->fh_info->Fill("Delta phi cut pi/",dijetDeltaPhiCut);

    // Initialize fh_events so that the bin order is correct
    fhistos->fh_events[0]->Fill("events",0.0);
    fhistos->fh_events[0]->Fill("particles",0.0);
    fhistos->fh_events[0]->Fill("acc. particles",0.0);
    fhistos->fh_events[0]->Fill("no rho calc. events",0.0);
    fhistos->fh_events[0]->Fill("rho calc. events",0.0);
    fhistos->fh_events[0]->Fill("jets",0.0);
    fhistos->fh_events[0]->Fill("acc. jets",0.0);
    fhistos->fh_events[0]->Fill("const. cut jets",0.0);
    fhistos->fh_events[0]->Fill("bg. subtr. jets",0.0);
    fhistos->fh_events[0]->Fill("bg. subtr. const. cut jets",0.0);
    fhistos->fh_events[0]->Fill("kt-jets",0.0);
    fhistos->fh_events[0]->Fill("acc. kt-jets",0.0);
    fhistos->fh_events[0]->Fill("leading jet drop",0.0);
    fhistos->fh_events[0]->Fill("subleading jet drop",0.0);
    fhistos->fh_events[0]->Fill("raw dijets",0.0);
    fhistos->fh_events[0]->Fill("raw dijets leading cut",0.0);
    fhistos->fh_events[0]->Fill("raw acc. dijets",0.0);
    fhistos->fh_events[0]->Fill("raw deltaphi cut dijets",0.0);
    fhistos->fh_events[0]->Fill("bg. subtr. dijets",0.0);
    fhistos->fh_events[0]->Fill("bg. subtr. dijets leading cut",0.0);
    fhistos->fh_events[0]->Fill("bg. subtr. acc. dijets",0.0);
    fhistos->fh_events[0]->Fill("bg. subtr. deltaphi cut dijets",0.0);
    fhistos->fh_events[0]->Fill("bg. subtr. const. cut dijets",0.0);
    fhistos->fh_events[0]->Fill("bg. subtr. const. cut dijets leading cut",0.0);
    fhistos->fh_events[0]->Fill("bg. subtr. const. cut acc. dijets",0.0);
    fhistos->fh_events[0]->Fill("bg. subtr. const. cut deltaphi cut dijets",0.0);
    fhistos->fh_events[0]->Fill("const. cut dijets",0.0);
    fhistos->fh_events[0]->Fill("const. cut dijets leading cut",0.0);
    fhistos->fh_events[0]->Fill("const. cut acc. dijets",0.0);
    fhistos->fh_events[0]->Fill("const. cut deltaphi cut dijets",0.0);
    fhistos->fh_events[0]->Fill("kt dijets",0.0);
    fhistos->fh_events[0]->Fill("kt dijets leading cut",0.0);
    fhistos->fh_events[0]->Fill("kt acc. dijets",0.0);
    fhistos->fh_events[0]->Fill("kt deltaphi cut dijets",0.0);


    //--------------------------------------------------------
    //         B e g i n    e v e n t    l o o p.
    //--------------------------------------------------------
    cout<<"Let's start" <<endl; 
    //int ieout = nEvent/20;
    //if (ieout<1) ieout=1;
    int EventCounter = 0;
    //Float_t sigmaGen = 0.0;
    Float_t ebeweight = 1.0;

    HepMC::GenEvent* evt = ascii_in.read_next_event();

    //int iEvent=0;
    while ( evt ) {
        //iEvent++;
        inputList->Clear("C");
        ebeweight = 1.0; //no event-by-event weight at all. //sigmaGen/nTrial;
        hCrossSectionInfo->Fill(7.5,ebeweight);
        //if(iEvent % ieout == 0) cout << iEvent << "\t" << int(float(iEvent)/nEvent*100) << "%, nTried:" << nTried << ", nTrial:" << nTrial << ", sigma:" << sigmaGen << endl;

        for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p ){
            if (    (*p)->status() == 1 &&                                                 // Check if particle is final
                    datacol.particle(HepPDT::ParticleID((*p)->pdg_id()))->charge() != 0 && // Only charged particles are used.
                    datacol.particle(HepPDT::ParticleID((*p)->pdg_id()))->isHadron() ) {   // Only hadrons are used.
                TLorentzVector lParticle((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
                AliJBaseTrack track( lParticle );
                track.SetID((*p)->pdg_id());
                track.SetTrackEff(1.);
                new ((*inputList)[inputList->GetEntriesFast()]) AliJBaseTrack(track);
            }
        } // end of finalparticles

        // Here I call my function
        CalculateJetsDijets(inputList,
                            5, // Debug
                            centBin, // Cent bin
                            partMinEtaCut, // Particle eta cut
                            partMinPtCut, // Particle pt cut
                            coneR, // Jet cone size
                            ktconeR,
                            fktScheme,  
                            fusePionMassInktjets,
                            fuseDeltaPhiBGSubtr,
                            jetConstituentCut, // Jet constituent cut
                            dijetLeadingPt, // Dijet leading jet pt cut
                            dijetSubleadingPt, // Dijet subleading jet pt cut
                            dijetDeltaPhiCut);  // Dijet DeltaPhi cut is pi/(this-argument)

        delete evt;
        ascii_in >> evt;
        EventCounter++;
        //if(iEvent == nEvent-1) cout << nEvent << "\t" << "100%, nTried:" << pythia.info.nTried() << ", sigma:" << pythia.info.sigmaGen() << endl ;
    }//event loop

    fhistos->fh_events[0]->Fill("events",EventCounter);

    /*
       nTried = pythia.info.nTried();
       nAccepted = pythia.info.nAccepted();
       sigmaGen = pythia.info.sigmaGen();
       hCrossSectionInfo->Fill(0.5,nTried);
       cout << "nTried after loop:" << nTried << endl;// print also inside event loop and in the macro.
       hCrossSectionInfo->Fill(1.5,nAccepted);
    //cout << "nAccepted after loop:" << nAccepted << endl;
    hCrossSectionInfo->Fill(2.5,sigmaGen);
    cout << "sigma after loop:" << sigmaGen << endl;
    hCrossSectionInfo->Fill(3.5,EventCounter);
    hCrossSectionInfo->Fill(4.5,energy);
    hCrossSectionInfo->Fill(5.5,1); // for counting # of merged
    hCrossSectionInfo->Fill(6.5,pythia.info.weightSum()); // for counting # of merged
    */

    fout->Write();
    fout->Close();
    cout << EventCounter << " events are analyzed successfully."<< endl;
    timer.Print(); 
    return 0;
}


//______________________________________________________________________________
void CalculateJetsDijets(TClonesArray *inList,
                         int    lDebug,
                         int    lCBin,
                         double lParticleEtaCut,
                         double lParticlePtCut,
                         double lJetCone,
                         double lktJetCone,
                         int    lktScheme,
                         bool   lusePionMassInkt,
                         bool   luseDeltaPhiBGSubtr,
                         double lConstituentCut,
                         double lLeadingJetCut,
                         double lSubleadingJetCut,
                         double lDeltaPhiCut) {

	double const etaMaxCutForJet = lParticleEtaCut-lJetCone;
	double const MinJetPt = 10.0; // Min Jet Pt cut to disregard low pt jets
    double const ghost_maxrap = lParticleEtaCut;
    unsigned int const repeat = 1; // default
    double const ghost_area   = 0.005; // ALICE=0.005 // default=0.01
    double const pionmass = 0.139570; // From PDG
    enum jetClasses {iRaw, iBGSubtr, iBGSubtrConstCut, iConstCut, iktJets, jetClassesSize};

    TString sDijetTypes[jetClassesSize] = {"raw", "bg. subtr.", "bg. subtr. const. cut", "const. cut", "kt"};

    double phi, eta, pt, pt2, rho, rhom, area, mjj, ptpair, dPhi, dPhi2;
    bool leadingTrackOverThreshold = false;
	vector<fastjet::PseudoJet> chparticles;
	vector<fastjet::PseudoJet> ktchparticles;
	vector<fastjet::PseudoJet> jets[jetClassesSize];
	vector<fastjet::PseudoJet> rhoEstJets;
	vector<fastjet::PseudoJet> constituents;
    fastjet::RecombinationScheme ktScheme;
	fastjet::PseudoJet jetAreaVector;
	fastjet::PseudoJet jet_bgSubtracted;
    fastjet::PseudoJet dijet;

	//--------------------------------------------------------
	//         B e g i n    e v e n t    l o o p.
	//--------------------------------------------------------
	int noTracks = inList->GetEntries();

    chparticles.clear();
    for (int itrack = 0; itrack < noTracks; ++itrack) {//loop over all the particles in the event
        // Building input particle list for the jet reconstruction
		AliJBaseTrack *trk = (AliJBaseTrack*)inList->At(itrack);
		pt = trk->Pt();
		eta = trk->Eta();
        fhistos->fh_events[lCBin]->Fill("particles",1.0);
        if (pt>lParticlePtCut && TMath::Abs(eta) < lParticleEtaCut){
            fhistos->fh_events[lCBin]->Fill("acc. particles",1.0);
            phi = trk->Phi();
            fhistos->fh_eta[lCBin]->Fill(eta);
            fhistos->fh_phi[lCBin]->Fill(phi);
            fhistos->fh_etaPhi[lCBin]->Fill(eta,phi);
            fhistos->fh_pt[lCBin]->Fill(pt);
            chparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), trk->E()));
            if(lusePionMassInkt) ktchparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), TMath::Sqrt(trk->Px()*trk->Px() + trk->Py()*trk->Py() + trk->Pz()*trk->Pz() + pionmass*pionmass)));            
            else ktchparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), trk->E()));
        }
    }
    if(chparticles.size()==0) return; // We are not intereted in empty events.

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Run the clustering, Reconstruct jets
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    switch (lktScheme) {
        case 0:  ktScheme = fastjet::E_scheme;
                 break;
        case 1:  ktScheme = fastjet::pt_scheme;
                 break;
        case 2:  ktScheme = fastjet::pt2_scheme;
                 break;
        case 3:  ktScheme = fastjet::Et_scheme;
                 break;
        case 4:  ktScheme = fastjet::Et2_scheme;
                 break;
        case 5:  ktScheme = fastjet::BIpt_scheme;
                 break;
        case 6:  ktScheme = fastjet::BIpt2_scheme;
                 break;
        default: ktScheme = fastjet::external_scheme;
                 ::Error("AliJCDijetTask","Unknown recombination scheme!");
    }
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, lJetCone, fastjet::pt_scheme); //Other option: fastjet::E_scheme
    fastjet::JetDefinition jet_def_bge(fastjet::kt_algorithm, lktJetCone, ktScheme);

    fastjet::GhostedAreaSpec const area_spec(ghost_maxrap, repeat, ghost_area);
    fastjet::AreaDefinition const area_def(fastjet::active_area, area_spec);
    fastjet::AreaDefinition const area_def_bge(fastjet::active_area_explicit_ghosts, area_spec);

    // Selector selects first all jets inside rapidity acceptance and then all but two hardest jets.
    fastjet::Selector const selectorAllButTwo = (!fastjet::SelectorNHardest(2));
    fastjet::Selector const selectorEta = fastjet::SelectorAbsEtaMax(ghost_maxrap - lktJetCone);
    fastjet::Selector const selectorBoth = selectorAllButTwo * selectorEta; // Here right selector is applied first, then the left one.
    fastjet::JetMedianBackgroundEstimator bge(selectorEta, jet_def_bge, area_def_bge);

    fastjet::ClusterSequenceArea cs(chparticles, jet_def, area_def);
    fastjet::ClusterSequenceArea cs_bge(ktchparticles, jet_def_bge, area_def_bge);
    
    jets[iRaw]    = fastjet::sorted_by_pt(cs.inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet
    jets[iktJets] = fastjet::sorted_by_pt(cs_bge.inclusive_jets(0.0)); // APPLY Min pt cut for jet
    
    if(luseDeltaPhiBGSubtr) {
        bool removed = false;
        for (unsigned iktJet = 1; iktJet < jets[iktJets].size(); iktJet++) { // First jet is already skipped here.
            if (!removed && TMath::Abs(jets[iktJets][iktJet].delta_phi_to(jets[iktJets][0])) > TMath::Pi()/lDeltaPhiCut) {
                removed = true;
                continue;
            }
            rhoEstJets.push_back(jets[iktJets][iktJet]); 
        }
    } else {
        rhoEstJets = selectorBoth(jets[iktJets]);
    }

    if( rhoEstJets.size() < 1 ) {
        fhistos->fh_events[lCBin]->Fill("no rho calc. events",1.0);
        rho  = 0.0;
        rhom = 0.0;
    } else { 
        fhistos->fh_events[lCBin]->Fill("rho calc. events",1.0);
        bge.set_jets(rhoEstJets);
        rho  = bge.rho()<0   ? 0.0 : bge.rho();
        rhom = bge.rho_m()<0 ? 0.0 : bge.rho_m();
    }
    

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Loop over jets and fill various histos 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fhistos->fh_rho[lCBin]->Fill(rho);
    fhistos->fh_rhom[lCBin]->Fill(rhom);
    if(lDebug > 9) std::cout << "Testing: Rho_M = " << rhom << ", has_rho_m() = " << bge.has_rho_m() << std::endl;

    // anti-kt jets:
    for (unsigned ijet = 0; ijet < jets[iRaw].size(); ijet++) {
        eta = jets[iRaw][ijet].eta();
        fhistos->fh_events[lCBin]->Fill("jets",1.0);
        // anti-kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            fhistos->fh_events[lCBin]->Fill("acc. jets",1.0);
            pt = jets[iRaw][ijet].pt();
            phi = jets[iRaw][ijet].phi();
            area = jets[iRaw][ijet].area();
            jetAreaVector = jets[iRaw][ijet].area_4vector();
            fhistos->fh_jetEta[lCBin][iRaw]->Fill(eta);  
            fhistos->fh_jetPhi[lCBin][iRaw]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[lCBin][iRaw]->Fill(eta,phi - TMath::Pi());
            fhistos->fh_jetPt[lCBin][iRaw]->Fill(pt);
            fhistos->fh_jetArea[lCBin][iRaw]->Fill(area);
            fhistos->fh_jetAreaRho[lCBin][iRaw]->Fill(area*rho);

            leadingTrackOverThreshold=false;
            if(lDebug > 9) cout << "Jet i=" << ijet << ", jet pt=" << pt << endl;
            for(unsigned iconst=0;iconst<jets[iRaw][ijet].constituents().size(); iconst++) {
                if(lDebug > 9) cout << "Constituent i=" << iconst << ", constituent pt=" << jets[iRaw][ijet].constituents()[iconst].pt() << endl;
                if(jets[iRaw][ijet].constituents()[iconst].pt() > lConstituentCut) { // Jet leading constituent cut.
                    leadingTrackOverThreshold=true;
                    break;
                }
            }
            
            jet_bgSubtracted = fastjet::PseudoJet(jets[iRaw][ijet].px() -        rho * jetAreaVector.px(),
                                                  jets[iRaw][ijet].py() -        rho * jetAreaVector.py(),
                                                  jets[iRaw][ijet].pz() - (rho+rhom) * jetAreaVector.pz(),
                                                  jets[iRaw][ijet].E()  - (rho+rhom) * jetAreaVector.E());

            if(leadingTrackOverThreshold) {
                fhistos->fh_events[lCBin]->Fill("const. cut jets",1.0);
                fhistos->fh_jetEta[lCBin][iConstCut]->Fill(eta);
                fhistos->fh_jetPhi[lCBin][iConstCut]->Fill(phi - TMath::Pi());
                fhistos->fh_jetEtaPhi[lCBin][iConstCut]->Fill(eta,phi - TMath::Pi());
                fhistos->fh_jetPt[lCBin][iConstCut]->Fill(pt);
                fhistos->fh_jetArea[lCBin][iConstCut]->Fill(area);
                fhistos->fh_jetAreaRho[lCBin][iConstCut]->Fill(area*rho);

                jets[iConstCut].push_back(jet_bgSubtracted);
            }


            // Check eta acceptance also for bg subtracted jets.
            eta = jet_bgSubtracted.eta();
            if(TMath::Abs(eta) < etaMaxCutForJet) {
                fhistos->fh_events[lCBin]->Fill("bg. subtr. jets",1.0);
                pt2 = jet_bgSubtracted.pt();
                phi = jet_bgSubtracted.phi();
                if(ijet==0 && pt>lLeadingJetCut && pt2<=lLeadingJetCut)       fhistos->fh_events[lCBin]->Fill("leading jet drop",1.0);
                if(ijet==1 && pt>lSubleadingJetCut && pt2<=lSubleadingJetCut) fhistos->fh_events[lCBin]->Fill("subleading jet drop",1.0);
                fhistos->fh_jetEta[lCBin][iBGSubtr]->Fill(eta);
                fhistos->fh_jetPhi[lCBin][iBGSubtr]->Fill(phi - TMath::Pi());
                fhistos->fh_jetEtaPhi[lCBin][iBGSubtr]->Fill(eta,phi - TMath::Pi());
                fhistos->fh_jetPt[lCBin][iBGSubtr]->Fill(pt2);
                fhistos->fh_jetArea[lCBin][iBGSubtr]->Fill(area); // Assuming bg subtracted jet has the same area.
                fhistos->fh_jetAreaRho[lCBin][iBGSubtr]->Fill(area*rho);

                jets[iBGSubtr].push_back(jet_bgSubtracted);

                if(leadingTrackOverThreshold) {
                    fhistos->fh_events[lCBin]->Fill("bg. subtr. const. cut jets",1.0);
                    fhistos->fh_jetEta[lCBin][iBGSubtrConstCut]->Fill(eta);
                    fhistos->fh_jetPhi[lCBin][iBGSubtrConstCut]->Fill(phi - TMath::Pi());
                    fhistos->fh_jetEtaPhi[lCBin][iBGSubtrConstCut]->Fill(eta,phi - TMath::Pi());
                    fhistos->fh_jetPt[lCBin][iBGSubtrConstCut]->Fill(pt2);
                    fhistos->fh_jetArea[lCBin][iBGSubtrConstCut]->Fill(area);
                    fhistos->fh_jetAreaRho[lCBin][iBGSubtrConstCut]->Fill(area*rho);

                    jets[iBGSubtrConstCut].push_back(jet_bgSubtracted);
                }
            }
        }
    }//end of the anti-kt-jet loop

    for (unsigned ijet = 0; ijet < jets[iktJets].size(); ijet++) {
        eta = jets[iktJets][ijet].eta();
        fhistos->fh_events[lCBin]->Fill("kt-jets",1.0);
        // kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            fhistos->fh_events[lCBin]->Fill("acc. kt-jets",1.0);
            pt = jets[iktJets][ijet].pt();
            phi = jets[iktJets][ijet].phi();
            area = jets[iktJets][ijet].area();
            fhistos->fh_jetEta[lCBin][iktJets]->Fill(eta);  
            fhistos->fh_jetPhi[lCBin][iktJets]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[lCBin][iktJets]->Fill(eta,phi - TMath::Pi());
            fhistos->fh_jetPt[lCBin][iktJets]->Fill(pt);
            fhistos->fh_jetArea[lCBin][iktJets]->Fill(area);
            fhistos->fh_jetAreaRho[lCBin][iktJets]->Fill(area*rho);
        }
    } //end of the kt-jet loop



    // Dijet calculations 
    for(int idijet=0; idijet < jetClassesSize; idijet++) {
        if(jets[idijet].size()>1) {
            jets[idijet] = fastjet::sorted_by_pt(jets[idijet]); // Sort in case of bg subtr messed up the order.
            fhistos->fh_events[lCBin]->Fill(Form("%s dijets",sDijetTypes[idijet].Data()),1.0);
            if(jets[idijet][0].pt()>lLeadingJetCut) {
                fhistos->fh_events[lCBin]->Fill(Form("%s dijets leading cut",sDijetTypes[idijet].Data()),1.0);
                if(jets[idijet][1].pt()>lSubleadingJetCut) {
                    fhistos->fh_events[lCBin]->Fill(Form("%s acc. dijets",sDijetTypes[idijet].Data()),1.0);
                    dijet = jets[idijet][0] + jets[idijet][1];
                    mjj = dijet.m();
                    ptpair = dijet.pt();
                    fhistos->fh_dijetInvM[lCBin][idijet]->Fill(mjj);
                    fhistos->fh_dijetPtPair[lCBin][idijet]->Fill(ptpair);
                    dPhi = jets[idijet][1].delta_phi_to(jets[idijet][0]);
                    dPhi2  = dPhi<0 ? dPhi+TMath::TwoPi() : dPhi;
                    fhistos->fh_dijetDeltaPhi[lCBin][idijet]->Fill(dPhi2);

                    // If subleading jet is on the opposite hemisphere compared to leading jet.
                    if(TMath::Abs(dPhi2 - TMath::Pi()) < TMath::Pi()/lDeltaPhiCut) {
                        fhistos->fh_events[lCBin]->Fill(Form("%s deltaphi cut dijets",sDijetTypes[idijet].Data()),1.0);
                        fhistos->fh_dijetInvMDeltaPhiCut[lCBin][idijet]->Fill(mjj); 
                        fhistos->fh_dijetPtPairDeltaPhiCut[lCBin][idijet]->Fill(ptpair); 
                    }
                }
            }
        }
    }
}
