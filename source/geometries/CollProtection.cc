#include "CollProtection.h"
#include "Visibilities.h"

#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4VisAttributes.hh>

namespace nexus {
  using namespace CLHEP;

 CollProtection::CollProtection() : _axis_centre(0.*mm)
  {
  }

  CollProtection::~CollProtection()
  {
  }

  void CollProtection::Construct()
  {
    G4double alum_thick = 9. * mm;
    G4Box* protection_solid = 
      new G4Box("SOURCE_PROTECTION", 100.*mm/2., 100.*mm/2.,  alum_thick/2.);
    G4LogicalVolume* protection_logic =
      new G4LogicalVolume(protection_solid, G4NistManager::Instance()->FindOrBuildMaterial("G4_Al"), "SOURCE_PROTECTION");
    this->SetLogicalVolume(protection_logic);

    _axis_centre = alum_thick/2.;

    G4VisAttributes red_col = nexus::Red();
    red_col.SetForceSolid(true);
    protection_logic->SetVisAttributes(red_col);
  }

   G4double CollProtection::GetAxisCentre()
  {
    return _axis_centre;
  }
}