// ----------------------------------------------------------------------------
// nexus | Next100MuonVeto.cc
//
// Muon veto surrounding Lead castle placed at the LSC.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Next100MuonVeto.h"
#include "MaterialsList.h"
#include "Visibilities.h"
#include "BoxPointSampler.h"

#include <G4GenericMessenger.hh>
#include <G4SubtractionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <G4NistManager.hh>
#include <G4Material.hh>
#include <Randomize.hh>
#include <G4TransportationManager.hh>
#include <G4RotationMatrix.hh>
#include <G4UserLimits.hh>

#include <G4SDManager.hh>
#include "IonizationSD.h"

#include <CLHEP/Units/SystemOfUnits.h>

namespace nexus {

  using namespace CLHEP;

  Next100MuonVeto::Next100MuonVeto():
    GeometryBase{},

    // dimensions of the muon veto panels
    
    wall_x_ {4.0 * cm},
    wall_y_ {207.5 * cm},
    wall_z_ {418.0 * cm},

    door_x_ {198.4 * cm},
    door_y_ {207.5 * cm},
    door_z_ {4.0 * cm},

    roof_x_ {258.8  * cm},
    roof_y_ {4.0 * cm},
    roof_z_ {418.0 * cm},
    
    visibility_ {0},
    verbosity_{false}

  {

    msg_ = new G4GenericMessenger(this, "/Geometry/Next100/", "Control commands of geometry Next100.");
    msg_->DeclareProperty("muonveto_vis", visibility_, "Muon Veto Visibility");
    msg_->DeclareProperty("muonveto_verbosity", verbosity_, "Verbosity");


  }


  void Next100MuonVeto::Construct()
  {

    //! air box as mother volume
    G4Box* veto_box_solid = new G4Box("VETO_BOX", door_x_/2., wall_y_/2., wall_z_/2.);
    G4LogicalVolume* veto_box_logic =
      new G4LogicalVolume(veto_box_solid,
			  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "lab_air");
    this->SetLogicalVolume(veto_box_logic);

    //! define muon veto planes
    G4Box* wallr_box_solid = new G4Box("WALLR_BOX", wall_x_/2., wall_y_/2., wall_z_/2.);
    G4Box* walll_box_solid = new G4Box("WALLL_BOX", wall_x_/2., wall_y_/2., wall_z_/2.);
    G4Box* doorf_box_solid = new G4Box("DOORF_BOX", door_x_/2., door_y_/2., door_z_/2.);
    G4Box* doorb_box_solid = new G4Box("DOORB_BOX", door_x_/2., door_y_/2., door_z_/2.);
    G4Box* roof_box_solid  = new G4Box("ROOF_BOX", roof_x_/2., roof_y_/2., roof_z_/2.);
    
    G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUEN");
    //G4_POLYSTYRENE  
      
    G4LogicalVolume* wallr_box_logic = new G4LogicalVolume(wallr_box_solid,mat,"WALLR_BOX");
    G4LogicalVolume* walll_box_logic = new G4LogicalVolume(walll_box_solid,mat,"WALLL_BOX");
    G4LogicalVolume* doorf_box_logic = new G4LogicalVolume(doorf_box_solid,mat,"DOORF_BOX");
    G4LogicalVolume* doorb_box_logic = new G4LogicalVolume(doorb_box_solid,mat,"DOORB_BOX");
    G4LogicalVolume* roof_box_logic  = new G4LogicalVolume(roof_box_solid,mat,"ROOF_BOX");

    // location of planes w.r.t. to mother volume

    new G4PVPlacement(0, G4ThreeVector(door_x_/2-wall_x_/2,0,0.),
                      wallr_box_logic, "MUON_VETO", veto_box_logic, false, 0, false);
    new G4PVPlacement(0, G4ThreeVector(-door_x_/2+wall_x_/2,0,0.),
                      walll_box_logic, "MUON_VETO", veto_box_logic, false, 0, false);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,wall_z_/2-door_z_/2),
                      doorf_box_logic, "MUON_VETO", veto_box_logic, false, 0, false);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,-wall_z_/2+door_z_/2),
                      doorb_box_logic, "MUON_VETO", veto_box_logic, false, 0, false);
    new G4PVPlacement(0, G4ThreeVector(0.,wall_y_/2-roof_y_/2,0.),
                      roof_box_logic, "MOUN_VETO", veto_box_logic, false, 0, false);

    //Set the volume sensitive detector. 

    IonizationSD* rwallsd = new IonizationSD("/NEXT100/MUON_VETO_RWALL");
    wallr_box_logic->SetSensitiveDetector(rwallsd);
    G4SDManager::GetSDMpointer()->AddNewDetector(rwallsd);
    IonizationSD* lwallsd = new IonizationSD("/NEXT100/MUON_VETO_LWALL");
    walll_box_logic->SetSensitiveDetector(lwallsd);
    G4SDManager::GetSDMpointer()->AddNewDetector(lwallsd);
    
    IonizationSD* fdoorsd = new IonizationSD("/NEXT100/MUON_VETO_FDOOR");
    doorf_box_logic->SetSensitiveDetector(fdoorsd);
    G4SDManager::GetSDMpointer()->AddNewDetector(fdoorsd);
    IonizationSD* bdoorsd = new IonizationSD("/NEXT100/MUON_VETO_BDOOR");
    doorb_box_logic->SetSensitiveDetector(bdoorsd);
    G4SDManager::GetSDMpointer()->AddNewDetector(bdoorsd);

    IonizationSD* roofsd = new IonizationSD("/NEXT100/MUON_VETO_ROOF");
    roof_box_logic->SetSensitiveDetector(roofsd);
    G4SDManager::GetSDMpointer()->AddNewDetector(roofsd);
    
    // visivility settings
    if (visibility_) { 
      wallr_box_logic->SetVisAttributes(nexus::LightGrey());
      walll_box_logic->SetVisAttributes(nexus::LightGrey());
      doorf_box_logic->SetVisAttributes(nexus::LightGrey());
      doorb_box_logic->SetVisAttributes(nexus::LightGrey());
      roof_box_logic->SetVisAttributes(nexus::LightGrey());
      veto_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    }
    else {
      wallr_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
      walll_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
      doorf_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
      doorb_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
      roof_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    }



    if (verbosity_){
      std::cout<<"Muon Veto Dimensions (cm)" <<std::endl;
      std::cout<<"RIGHT WALL: "<< wall_x_/cm << " x " << wall_y_/cm << " x " << wall_z_/cm <<std::endl;
      std::cout<<"LEFT WALL: "<< wall_x_/cm << " x " << wall_y_/cm << "x" << wall_z_/cm <<std::endl;
      std::cout<<"FRONT DOOR: "<< door_x_/cm << " x " << door_y_/cm << " x " << door_z_/cm <<std::endl;
      std::cout<<"BACK DOOR DOOR: "<< door_x_/cm << " x " << door_y_/cm << "x" << door_z_/cm <<std::endl;
      std::cout<<"ROOF: "<< roof_x_/cm << " x " << roof_y_/cm << " x " << roof_z_/cm <<std::endl;
    }
  }



  Next100MuonVeto::~Next100MuonVeto()
  {
    
  }

 
} //end namespace nexus
