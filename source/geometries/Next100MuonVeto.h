// ----------------------------------------------------------------------------
// nexus | Next100MuonVeto.h
//
// Muon Veto surrounding castle placed at the LSC.
// Simplified version: fully active planes
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef NEXT100_MUONVETO_H
#define NEXT100_MUONVETO_H

#include "GeometryBase.h"

#include <G4Navigator.hh>

class G4GenericMessenger;


namespace nexus {

  //class BoxPointSampler;

  class Next100MuonVeto: public GeometryBase
  {
  public:
    /// Constructor
    Next100MuonVeto();

    /// Destructor
    ~Next100MuonVeto();

    /// Builder
    void Construct();


  private:

    // Dimensions
    G4double wall_x_, wall_y_, wall_z_;
    G4double door_x_, door_y_, door_z_;
    G4double roof_x_, roof_y_, roof_z_;
    
    G4bool visibility_;
    G4bool verbosity_;

    // Messenger for the definition of control commands
    G4GenericMessenger* msg_;

  };

 
} // end namespace nexus

#endif
