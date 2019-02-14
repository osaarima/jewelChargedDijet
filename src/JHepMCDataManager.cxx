#include  "JHepMCDataManager.h"


//______________________________________________________________________________
JHepMCDataManager::JHepMCDataManager(char* hepMCFile):
	fascii_in(NULL), 
	fdatacol(NULL), 
	fevt(NULL)
{
	// constructor
	fascii_in = new HepMC::IO_GenEvent(hepMCFile,std::ios::in);
    fdatacol = new HepPDT::ParticleDataTable( "Generic Particle Table" );

    // Initialize HepPDF so that particle ID's can be used to determine charge and hadrons
    std::ifstream pdfile( "particle.tbl" );
    if( !pdfile ) { 
      std::cerr << "cannot open particle.tbl" << std::endl;
      exit(-1);
    }
    {
        // Construct table builder
        HepPDT::TableBuilder  tb(*fdatacol);
        // bool  addParticleTable( std::istream&, TableBuilder&,
        //                         bool validate = false );
        // where:  validate=true => verify that the ParticleID is valid
        if( !addParticleTable( pdfile, tb, true ) ) { 
            std::cout << "error reading PDG pdt file " << std::endl; 
        }
    } // the tb destructor fills fdatacol

}

//______________________________________________________________________________
JHepMCDataManager::~JHepMCDataManager(){
	//if( fascii_in ) delete fascii_in;
	//if( fdatacol ) delete fdatacol;
	//if( fevt ) delete fevt;
}

//______________________________________________________________________________

bool JHepMCDataManager::next_event(){    
    bool breturn;
    if(fevt) delete fevt;
    fevt = fascii_in->read_next_event();
    if(fevt) breturn = true;
    else breturn = false;
    return breturn; 
}


//______________________________________________________________________________
void JHepMCDataManager::getList(TClonesArray* listToFill) { 

    listToFill->Clear("C");
    int counter = 0;

    for ( HepMC::GenEvent::particle_const_iterator p = fevt->particles_begin(); p != fevt->particles_end(); ++p ){
        if (    (*p)->status() == 1 &&                                                  // Check if particle is final
                fdatacol->particle(HepPDT::ParticleID((*p)->pdg_id()))->charge() != 0 && // Only charged particles are used.
                fdatacol->particle(HepPDT::ParticleID((*p)->pdg_id()))->isHadron() ) {   // Only hadrons are used.
            TLorentzVector lParticle((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
            AliJBaseTrack track( lParticle );
            track.SetID((*p)->pdg_id());
            track.SetTrackEff(1.);
            new ((*listToFill)[counter++]) AliJBaseTrack(track);
        }
    } // end of finalparticles
}
