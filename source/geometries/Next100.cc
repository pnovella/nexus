// ----------------------------------------------------------------------------
// nexus | Next100.cc
//
// Main class that constructs the geometry of the NEXT-100 detector.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "Next100.h"
#include "BoxPointSampler.h"
#include "MuonsPointSampler.h"
#include "LSCHallA.h"
#include "Next100Shielding.h"
#include "Next100Vessel.h"
#include "Next100Ics.h"
#include "Next100InnerElements.h"
//#include "Next100MuonVeto.h"
#include "FactoryBase.h"

// for muon veto 
#include "IonizationSD.h"
#include <G4SDManager.hh>
#include "Visibilities.h"
//---

#include <G4GenericMessenger.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4NistManager.hh>
#include <G4UserLimits.hh>


namespace nexus {

  REGISTER_CLASS(Next100, GeometryBase)

  using namespace CLHEP;

  Next100::Next100():
    GeometryBase(),
    // Lab dimensions
    //lab_size_ (5. * m),
    lab_size_ (10. * m), // increased to hold muon veto

    // common used variables in geomety components
    // 0.1 mm grid thickness
    // note that if grid thickness change it must be also changed in Next100FieldCage.cc
    gate_tracking_plane_distance_((26.1 + 0.1)   * mm),
    gate_sapphire_wdw_distance_  ((1458.2 - 0.1) * mm),
    
    specific_vertex_{},
    lab_walls_(false),
    veto_walls_(false),

    //--- muon veto dims and conf -----//
    mv_wall_x_ {4.0 * cm},
    mv_wall_y_ {207.5 * cm},
    mv_wall_z_ {418.0 * cm},
    mv_door_x_ {198.4 * cm},
    mv_door_y_ {207.5 * cm},
    mv_door_z_ {4.0 * cm},
    mv_roof_x_ {198.4  * cm},
    mv_roof_y_ {4.0 * cm},
    mv_roof_z_ {418.0 * cm},
    mv_visibility_(false),
    mv_verbosity_(false)
    //------------------------//

  {

    msg_ = new G4GenericMessenger(this, "/Geometry/Next100/",
				  "Control commands of geometry Next100.");

    msg_->DeclarePropertyWithUnit("specific_vertex", "mm",  specific_vertex_,
      "Set generation vertex.");

    msg_->DeclareProperty("lab_walls", lab_walls_, "Placement of Hall A walls");

    // muon veto properties
    msg_->DeclareProperty("veto_walls", veto_walls_, "Placement of Muon Veto");
    msg_->DeclareProperty("muonveto_vis", mv_visibility_, "Muon Veto Visibility");
    msg_->DeclareProperty("muonveto_verbosity", mv_verbosity_, "Muon VetoVerbosity");
    
   // The following methods must be invoked in this particular
  // order since some of them depend on the previous ones

  // Shielding
  shielding_ = new Next100Shielding();

  //Lab walls
  hallA_walls_ = new LSCHallA();

  // Vessel
  vessel_ = new Next100Vessel();

  // Internal copper shielding
  ics_ = new Next100Ics();

  // Inner Elements
  inner_elements_ = new Next100InnerElements();

  // Muon Veto
  //muon_veto_ = new Next100MuonVeto();
  
  }


  Next100::~Next100()
  {
    delete inner_elements_;
    delete ics_;
    delete vessel_;
    delete shielding_;
    delete lab_gen_;
    delete hallA_walls_;
    //delete muon_veto_;
  }


  void Next100::Construct()
  {
    // LAB /////////////////////////////////////////////////////////////
    // This is just a volume of air surrounding the detector so that
    // events (from calibration sources or cosmic rays) can be generated
    // on the outside.
    if (lab_walls_){
      // We want to simulate the walls (for muons in most cases).
      hallA_walls_->Construct();
      hallA_logic_ = hallA_walls_->GetLogicalVolume();
      G4double hallA_length = hallA_walls_->GetLSCHallALength();
      // Since the walls will be displaced need to make the
      // "lab" double sized to be sure.
      G4Box* lab_solid = new G4Box("LAB", hallA_length, hallA_length, hallA_length);
      G4Material *vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
      lab_logic_ = new G4LogicalVolume(lab_solid, vacuum, "LAB");
      this->SetSpan(2 * hallA_length);
    }
    else {
      G4Box* lab_solid = new G4Box("LAB", lab_size_/2., lab_size_/2., lab_size_/2.);
      lab_logic_ = new G4LogicalVolume(lab_solid,
        G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "LAB");
    }
    lab_logic_->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Set this volume as the wrapper for the whole geometry
    // (i.e., this is the volume that will be placed in the world)
    this->SetLogicalVolume(lab_logic_);

    // VESSEL (initialize first since it defines EL position)
    vessel_->SetELtoTPdistance(gate_tracking_plane_distance_);
    vessel_->Construct();
    G4LogicalVolume* vessel_logic = vessel_->GetLogicalVolume();
    G4LogicalVolume* vessel_internal_logic  = vessel_->GetInternalLogicalVolume();
    G4VPhysicalVolume* vessel_internal_phys = vessel_->GetInternalPhysicalVolume();
    G4ThreeVector vessel_displacement = shielding_->GetAirDisplacement(); // explained below
    gate_zpos_in_vessel_ = vessel_->GetELzCoord();

    // SHIELDING
    shielding_->Construct();
    shielding_->SetELzCoord(gate_zpos_in_vessel_);
    G4LogicalVolume* shielding_logic     = shielding_->GetLogicalVolume();
    G4LogicalVolume* shielding_air_logic = shielding_->GetAirLogicalVolume();

    // Recall that airbox is slighly displaced in Y dimension. In order to avoid
    // mistmatch with vertex generators, we place the vessel in the center of the world volume
    new G4PVPlacement(0, -vessel_displacement, vessel_logic,
                      "VESSEL", shielding_air_logic, false, 0);

    // INNER ELEMENTS
    inner_elements_->SetLogicalVolume(vessel_internal_logic);
    inner_elements_->SetPhysicalVolume(vessel_internal_phys);
    inner_elements_->SetELzCoord(gate_zpos_in_vessel_);
    inner_elements_->SetELtoSapphireWDWdistance(gate_sapphire_wdw_distance_);
    inner_elements_->SetELtoTPdistance         (gate_tracking_plane_distance_);
    inner_elements_->Construct();

    // INNER COPPER SHIELDING
    ics_->SetLogicalVolume(vessel_internal_logic);
    ics_->SetELzCoord(gate_zpos_in_vessel_);
    ics_->SetELtoSapphireWDWdistance(gate_sapphire_wdw_distance_);
    ics_->SetELtoTPdistance         (gate_tracking_plane_distance_);
    ics_->SetPortZpositions(vessel_->GetPortZpositions());
    ics_->Construct();

    G4ThreeVector gate_pos(0., 0., -gate_zpos_in_vessel_);
    if (lab_walls_){
      G4ThreeVector castle_pos(0., hallA_walls_->GetLSCHallACastleY(),
                               hallA_walls_->GetLSCHallACastleZ());

      new G4PVPlacement(0, castle_pos, shielding_logic,
                        "LEAD_BOX", hallA_logic_, false, 0);
      new G4PVPlacement(0, gate_pos - castle_pos, hallA_logic_,
                        "Hall_A", lab_logic_, false, 0, false);
    }
    else {
      new G4PVPlacement(0, gate_pos, shielding_logic, "LEAD_BOX", lab_logic_, false, 0);
    }

    //! muon veto !!!!
    if (veto_walls_){ ConstructMuonVeto(); }
   
    //// VERTEX GENERATORS
    lab_gen_ =
      new BoxPointSampler(lab_size_ - 1.*m, lab_size_ - 1.*m, lab_size_  - 1.*m, 1.*m,
                          G4ThreeVector(0., 0., 0.), 0);
  }


  G4ThreeVector Next100::GenerateVertex(const G4String& region) const
  {
    G4ThreeVector vertex(0.,0.,0.);

    // Air around shielding
    if (region == "LAB") {
      vertex = lab_gen_->GenerateVertex("INSIDE");
    }

    // Shielding regions
    else if ((region == "SHIELDING_LEAD")  ||
             (region == "SHIELDING_STEEL") ||
             (region == "INNER_AIR") ||
             (region == "EXTERNAL") ||
             (region == "SHIELDING_STRUCT") ||
             (region == "PEDESTAL") ||
             (region == "BUBBLE_SEAL") ||
             (region == "EDPM_SEAL")) {
      vertex = shielding_->GenerateVertex(region);
    }

    // Vessel regions
    else if ((region == "VESSEL")  ||
             (region == "PORT_1a") ||
             (region == "PORT_2a") ||
             (region == "PORT_1b") ||
             (region == "PORT_2b")) {
      vertex = vessel_->GenerateVertex(region);
    }

    // Inner copper shielding
    else if (region == "ICS"){
      vertex = ics_->GenerateVertex(region);
    }

    // Inner elements (photosensors' planes and field cage)
    else if ((region == "CENTER") ||
             (region == "ACTIVE") ||
             (region == "CATHODE_RING") ||
             (region == "BUFFER") ||
             (region == "XENON") ||
             (region == "LIGHT_TUBE") ||
             (region == "HDPE_TUBE") ||
             (region == "EL_GAP") ||
             (region == "EP_COPPER_PLATE") ||
             (region == "SAPPHIRE_WINDOW") ||
             (region == "OPTICAL_PAD") ||
             (region == "PMT_BODY") ||
             (region == "PMT") ||
             (region == "PMT_BASE") ||
             (region == "TP_COPPER_PLATE") ||
             (region == "SIPM_BOARD") ||
             (region == "DB_PLUG") ||
             (region == "EL_TABLE") ||
             (region == "FIELD_RING") ||
             (region == "GATE_RING") ||
             (region == "ANODE_RING") ||
             (region == "RING_HOLDER")) {
      vertex = inner_elements_->GenerateVertex(region);
    }

    else if (region == "AD_HOC") {
      // AD_HOC does not need to be shifted because it is passed by the user
      vertex = specific_vertex_;
      return vertex;
    }

    // Lab walls
    else if ((region == "HALLA_INNER") || (region == "HALLA_OUTER")){
      if (!lab_walls_)
        G4Exception("[Next100]", "GenerateVertex()", FatalException,
                    "This vertex generation region must be used with lab_walls == true!");
      vertex = hallA_walls_->GenerateVertex(region);
      while (vertex[1]<(-shielding_->GetHeight()/2.)){
        vertex = hallA_walls_->GenerateVertex(region);}
    }

    else {
      G4Exception("[Next100]", "GenerateVertex()", FatalException,
		  "Unknown vertex generation region!");
    }

    G4ThreeVector displacement = G4ThreeVector(0., 0., -gate_zpos_in_vessel_);
    vertex = vertex + displacement;

    return vertex;
  }

  void Next100::ConstructMuonVeto(){

    //muon_veto_->Construct();
    //G4LogicalVolume* muon_veto_logic = muon_veto_->GetLogicalVolume();
    //G4ThreeVector veto_dims = muon_veto_->GetDimensions();
    //new G4PVPlacement(0,veto_dims, muon_veto_logic,
    //		"MUON_VETO", lab_logic_, false, 0, true);

    //! define muon veto planes
    G4Box* wallr_box_solid = new G4Box("WALLR_BOX", mv_wall_x_/2., mv_wall_y_/2., mv_wall_z_/2.);
    G4Box* walll_box_solid = new G4Box("WALLL_BOX", mv_wall_x_/2., mv_wall_y_/2., mv_wall_z_/2.);
    G4Box* doorf_box_solid = new G4Box("DOORF_BOX", mv_door_x_/2., mv_door_y_/2., mv_door_z_/2.);
    G4Box* doorb_box_solid = new G4Box("DOORB_BOX", mv_door_x_/2., mv_door_y_/2., mv_door_z_/2.);
    G4Box* roof_box_solid  = new G4Box("ROOF_BOX", mv_roof_x_/2., mv_roof_y_/2., mv_roof_z_/2.);
    
    G4Material* mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE");
        
    G4LogicalVolume* wallr_box_logic = new G4LogicalVolume(wallr_box_solid,mat,"WALLR_BOX");
    G4LogicalVolume* walll_box_logic = new G4LogicalVolume(walll_box_solid,mat,"WALLL_BOX");
    G4LogicalVolume* doorf_box_logic = new G4LogicalVolume(doorf_box_solid,mat,"DOORF_BOX");
    G4LogicalVolume* doorb_box_logic = new G4LogicalVolume(doorb_box_solid,mat,"DOORB_BOX");
    G4LogicalVolume* roof_box_logic  = new G4LogicalVolume(roof_box_solid,mat,"ROOF_BOX");

    G4ThreeVector gate_pos(0., 0., -gate_zpos_in_vessel_);

    G4ThreeVector rel_pos;
    G4LogicalVolume* vol_logic;
    if (lab_walls_){
      vol_logic = hallA_logic_;
      rel_pos = G4ThreeVector(0., hallA_walls_->GetLSCHallACastleY(),hallA_walls_->GetLSCHallACastleZ());}
    else {
      vol_logic = lab_logic_;
      rel_pos = gate_pos;}
    
    new G4PVPlacement(0, G4ThreeVector(mv_door_x_/2,0,0.) + rel_pos,
		      wallr_box_logic, "MUON_VETO_RWALL", vol_logic, false, 0);
    new G4PVPlacement(0, G4ThreeVector(-mv_door_x_/2,0,0.) + rel_pos,
		      walll_box_logic, "MUON_VETO_LWALL", vol_logic, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,mv_wall_z_/2) + rel_pos,
		      doorf_box_logic, "MUON_VETO_FDOOR", vol_logic, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,-mv_wall_z_/2) + rel_pos,
		      doorb_box_logic, "MUON_VETO_BDOOR", vol_logic, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0.,mv_wall_y_/2,0.) + rel_pos,
		      roof_box_logic, "MOUN_VETO_ROOF", vol_logic, false, 0);

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
    if (mv_visibility_) { 
      wallr_box_logic->SetVisAttributes(nexus::LightGrey());
      walll_box_logic->SetVisAttributes(nexus::LightGrey());
      doorf_box_logic->SetVisAttributes(nexus::LightGrey());
      doorb_box_logic->SetVisAttributes(nexus::LightGrey());
      roof_box_logic->SetVisAttributes(nexus::LightGrey());
    }
    else {
      wallr_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
      walll_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
      doorf_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
      doorb_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
      roof_box_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    }
      
    if (mv_verbosity_){
      std::cout<<"Muon Veto Dimensions (cm)" <<std::endl;
      std::cout<<"RIGHT WALL: "<< mv_wall_x_/cm << " x " << mv_wall_y_/cm << " x " << mv_wall_z_/cm <<std::endl;
      std::cout<<"LEFT WALL: "<< mv_wall_x_/cm << " x " << mv_wall_y_/cm << "x" << mv_wall_z_/cm <<std::endl;
      std::cout<<"FRONT DOOR: "<< mv_door_x_/cm << " x " << mv_door_y_/cm << " x " << mv_door_z_/cm <<std::endl;
      std::cout<<"BACK DOOR DOOR: "<< mv_door_x_/cm << " x " << mv_door_y_/cm << "x" << mv_door_z_/cm <<std::endl;
      std::cout<<"ROOF: "<< mv_roof_x_/cm << " x " << mv_roof_y_/cm << " x " << mv_roof_z_/cm <<std::endl;
    }
    
  }

  
} //end namespace nexus
