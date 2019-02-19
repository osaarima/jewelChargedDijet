////////////////////////////////////////////////////
/*!
  \file AliJBaseEventHeader.cxx
  \brief
  \author O. Saarimaki, D.J.Kim (University of Jyvaskyla), Changwook Park (McGill University)
  \email: djkim@jyu.fi
  \version $Revision: 1.6 $
  \date $Date: 2018/08/22 15:36:00 $
  */
// Base class for event headers
////////////////////////////////////////////////////

#include <TNamed.h>
#include "AliJBaseEventHeader.h"

ClassImp(AliJBaseEventHeader)

//______________________________________________________________________________
AliJBaseEventHeader::AliJBaseEventHeader(): 
  TNamed("AliJBaseEventHeader", ""), 
  fEventID(-999),              
  fCentrality(-999),
  fVtxX(-999),
  fVtxY(-999),
  fVtxZ(-999),
  fVtxZErr(-999),
  fVtxMCX(9999), 
  fVtxMCY(9999), 
  fVtxMCZ(9999),
  fSigmaNN(9999),
  fSigmaGen(9999)
{
  // default constructor
}

//______________________________________________________________________________
AliJBaseEventHeader::AliJBaseEventHeader(int eventid, float cent, float vtxz): 
  TNamed("AliJBaseEventHeader", ""), 
  fEventID(eventid),              
  fCentrality(cent),
  fVtxX(-999),
  fVtxY(-999),
  fVtxZ(vtxz),
  fVtxZErr(-999),
  fVtxMCX(9999),
  fVtxMCY(9999), 
  fVtxMCZ(9999),
  fSigmaNN(9999),
  fSigmaGen(9999)
{
  //constructor
}
//______________________________________________________________________________
AliJBaseEventHeader::AliJBaseEventHeader(const AliJBaseEventHeader& a):
  TNamed(a),
  fEventID(a.fEventID),
  fCentrality(a.fCentrality),
  fVtxX(a.fVtxX),
  fVtxY(a.fVtxY),
  fVtxZ(a.fVtxZ),
  fVtxZErr(a.fVtxZErr),
  fVtxMCX(a.fVtxMCX), 
  fVtxMCY(a.fVtxMCY), 
  fVtxMCZ(a.fVtxMCZ),
  fSigmaNN(a.fSigmaNN),
  fSigmaGen(a.fSigmaGen)  
{
  //copy constructor
}

//______________________________________________________________________________
AliJBaseEventHeader&  AliJBaseEventHeader::operator=(const AliJBaseEventHeader& header){
  //operator=  
  if(this != &header){
    TNamed::operator=(header);
    fEventID    = header.fEventID;
    fCentrality = header.fCentrality;
    fVtxX       = header.fVtxX;
    fVtxY       = header.fVtxY;
    fVtxZ       = header.fVtxZ;
    fVtxZErr   = header.fVtxZErr;
    fVtxMCX     = header.fVtxMCX;
    fVtxMCY     = header.fVtxMCY;
    fVtxMCZ     = header.fVtxMCZ;
    fSigmaNN    = header.fSigmaNN;
    fSigmaGen   = header.fSigmaGen;
  }

  return *this;
}

