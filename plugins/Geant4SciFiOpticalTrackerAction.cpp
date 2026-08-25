#include "DD4hep/Version.h"
#include "DDG4/Geant4SensDetAction.inl"
#include "G4OpticalPhoton.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {

namespace sim {

  /*
   Sensitive detector action for the SciFi SiPM readout. Like Geant4OpticalTrackerAction, but draws against the sensitive 
   material's EFFICIENCY property (SiPM_PDE) and only registers a hit if the photon passes. A photon that fails is still killed here
  */
  struct SciFiOpticalTracker {}; // tag type for the Geant4SensitiveAction template

  /// Define collections created by this sensitive action object
  template <> void Geant4SensitiveAction<SciFiOpticalTracker>::defineCollections() {
    m_collectionID = declareReadoutFilteredCollection<Geant4Tracker::Hit>();
  }

  /// Method for generating hits using the information of G4Step object.
  template <> bool
  Geant4SensitiveAction<SciFiOpticalTracker>::process(const G4Step* step,
                                                       G4TouchableHistory* /* hist */) {
    Geant4StepHandler h(step);  //wrapper around G4Step
    bool is_optical = (h.trackDef() == G4OpticalPhoton::OpticalPhotonDefinition());
    if (is_optical) {
      step->GetTrack()->SetTrackStatus(fStopAndKill);
    }

    // Missing property table/EFFICIENCY entry -> PDE=1 (always accept).
    if (is_optical) {
      G4Material* mat = step->GetPreStepPoint()->GetMaterial();
      G4MaterialPropertiesTable* mpt = mat ? mat->GetMaterialPropertiesTable() : nullptr;
      if (mpt) {
        auto* eff = mpt->GetProperty("EFFICIENCY");
        if (eff) {
          G4double energy = step->GetTrack()->GetKineticEnergy();
          G4double pde    = eff->Value(energy);
          bool     passed_pde = (G4UniformRand() <= pde);
          always("+++ PDE lookup: energy(eV)=%f  pde=%f  passed=%d",
                 energy / CLHEP::eV, pde, passed_pde);
          if (!passed_pde) return true; // absorbed but not detected, no hit
        }
      }
    }

    typedef Geant4Tracker::Hit Hit;
    auto      contrib = Hit::extractContribution(step);
    Direction hit_momentum = 0.5 * (h.preMom() + h.postMom());
    double    hit_deposit  = contrib.deposit;
    Hit* hit = new Hit(contrib, hit_momentum, hit_deposit);

    hit->cellID = cellID(step);
    if ( 0 == hit->cellID )  {
      hit->cellID = volumeID( step ) ;
      except("Invalid CELL ID for hit!");
    }
    collection(m_collectionID)->add(hit);
    mark(h.track);
    print("Hit with deposit:%f  Pos:%f %f %f ID=%016X",
          hit->energyDeposit,hit->position.X(),hit->position.Y(),
          hit->position.Z(),(void*)hit->cellID);
    Geant4TouchableHandler handler(step);
    print("    Geant4 path:%s",handler.path().c_str());
    return true;
  }

  typedef Geant4SensitiveAction<SciFiOpticalTracker> Geant4SciFiOpticalTrackerAction;

} // namespace sim
} // namespace dd4hep

#include "DDG4/Factories.h"
DECLARE_GEANT4SENSITIVE(Geant4SciFiOpticalTrackerAction)
