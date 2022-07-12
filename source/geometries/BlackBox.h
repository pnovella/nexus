// ----------------------------------------------------------------------------
// nexus | BlackBox.cc
//
// Sapphire disk and DICE board with configurable mask thickness parallelly
// placed in a black box.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef BLACK_BOX_H
#define BLACK_BOX_H

#include "GeometryBase.h"
#include "BlackBoxSiPMBoard.h"

class G4Material;
class G4OpticalSurface;
class G4GenericMessenger;

namespace nexus {

  class BlackBoxSiPMBoard;

  class BlackBox: public GeometryBase
  {
  public:
    /// Constructor
    BlackBox();
    /// Destructor
    ~BlackBox();

    void SetMotherPhysicalVolume(G4VPhysicalVolume* mother_phys);

    /// Return vertex within region <region> of the chamber
    G4ThreeVector GenerateVertex(const G4String& region) const;


    void Construct();

  private:
    // Dimensions
    G4double world_z_;
    G4double world_xy_;
    G4double box_z_;
    G4double box_xy_;

    BlackBoxSiPMBoard* dice_;
    G4ThreeVector kdb_dimensions_;

    G4bool visibility_;

    //Messenger for configuration parameters
    G4GenericMessenger* msg_;

    G4ThreeVector specific_vertex_;
    //G4double specific_vertex_X_;
    //G4double specific_vertex_Y_;
    //G4double specific_vertex_Z_;
    G4double dice_board_x_pos_;
    G4double dice_board_y_pos_;
    G4double dice_board_z_pos_;
    G4double rotation_;
    G4double mask_thickn_;
    G4double membrane_thickn_;
    G4double coating_thickn_;
    G4bool   membrane_hole_;
    G4double membrane_hole_diam_;
    G4double membrane_hole_x_;
    G4double membrane_hole_y_;
    G4double hole_diameter_;
    G4double hole_x_;
    G4double hole_y_;
    G4VPhysicalVolume*  mother_phys_;

  };
  inline void BlackBox::SetMotherPhysicalVolume(G4VPhysicalVolume* mother_phys)
    { mother_phys_ = mother_phys; }

} // end namespace nexus

#endif
