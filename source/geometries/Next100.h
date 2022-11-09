// ----------------------------------------------------------------------------
// nexus | Next100.h
//
// Main class that constructs the geometry of the NEXT-100 detector.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef NEXT100_H
#define NEXT100_H

#include "GeometryBase.h"

class G4LogicalVolume;
class G4GenericMessenger;

namespace nexus {class BoxPointSampler;}


namespace nexus {

  class BoxPointSampler;
  class Next100Shielding;
  class Next100Vessel;
  class Next100Ics;
  class Next100InnerElements;
  class LSCHallA;
  //class Next100MuonVeto;
  
  class Next100: public GeometryBase
  {
  public:
    /// Constructor
    Next100();

    /// Destructor
    ~Next100();

    /// Generate a vertex within a given region of the geometry
    G4ThreeVector GenerateVertex(const G4String& region) const;


  private:
    void BuildLab();
    void Construct();
    void ConstructMuonVeto();

  private:
    // Detector dimensions
    const G4double lab_size_;          /// Size of the air box containing the detector
    const G4double gate_tracking_plane_distance_, gate_sapphire_wdw_distance_;

    // Pointers to logical volumes
    G4LogicalVolume* lab_logic_;
    G4LogicalVolume* buffer_gas_logic_;
    G4LogicalVolume* hallA_logic_;

    // Detector parts
    LSCHallA* hallA_walls_;
    Next100Shielding* shielding_;
    Next100Vessel*    vessel_;
    Next100Ics*       ics_;
    Next100InnerElements* inner_elements_;
    //Next100MuonVeto* muon_veto_;
    
    BoxPointSampler* lab_gen_; ///< Vertex generator

    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;

    /// Specific vertex for AD_HOC region
    G4ThreeVector specific_vertex_;

    /// Position of gate in its mother volume
    G4double gate_zpos_in_vessel_;

    /// Whether or not to build LSC HallA.
    G4bool lab_walls_;

    /// Whether or not to build muon veto.
    G4bool veto_walls_;

    /// Muon veto dimensions
    G4double mv_wall_x_, mv_wall_y_, mv_wall_z_;
    G4double mv_door_x_, mv_door_y_, mv_door_z_;
    G4double mv_roof_x_, mv_roof_y_, mv_roof_z_;
    /// Muon veto conf
    G4bool mv_visibility_;
    G4bool mv_verbosity_;
    
  };

} // end namespace nexus

#endif
