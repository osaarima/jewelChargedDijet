#ifndef JHEPMCDATAMANAGER_H
#define JHEPMCDATAMANAGER_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/GenEvent.h>
#include <HepPDT/defs.h>
#include <HepPDT/TableBuilder.hh>
#include <HepPDT/TempParticleData.hh>
#include <HepPDT/ParticleDataTable.hh>

#include  "AliJBaseTrack.h"

class JHepMCDataManager  {

	public:
		JHepMCDataManager(char* hepMCFile);

		virtual ~JHepMCDataManager();		                    //destructor

        virtual bool next_event();
		virtual void getList(TClonesArray* listToFill);

	protected:
        HepMC::IO_GenEvent *fascii_in;
        HepPDT::ParticleDataTable *fdatacol;
        HepMC::GenEvent* fevt;
};

#endif

