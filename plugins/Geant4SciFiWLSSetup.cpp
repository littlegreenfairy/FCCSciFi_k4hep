#include <DDG4/Geant4PhysicsList.h>
#include <G4OpWLS.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4ProcessManager.hh>

namespace dd4hep {
namespace sim {

  // Constructs G4OpWLS directly and adds it to the optical photon's process
  // manager by hand, instead of relying on Geant4OpticalPhotonPhysics's own
  // G4OpticalParameters-based activation -- which, in this stack, does not
  // result in G4OpWLS being constructed even when SetProcessActivation is
  // called immediately before Geant4OpticalPhotonPhysics checks it. This
  // sidesteps that mechanism entirely rather than depending on it.
  class Geant4SciFiWLSSetup : public Geant4PhysicsList {
  public:
    Geant4SciFiWLSSetup() = delete;
    Geant4SciFiWLSSetup(const Geant4SciFiWLSSetup&) = delete;
    Geant4SciFiWLSSetup(Geant4Context* ctxt, const std::string& nam)
      : Geant4PhysicsList(ctxt, nam) {}
    virtual ~Geant4SciFiWLSSetup() = default;

    virtual void constructProcesses(G4VUserPhysicsList* physics_list) {
      this->Geant4PhysicsList::constructProcesses(physics_list);

      G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
      if (0 == particle) {
        except("+++ Geant4SciFiWLSSetup: cannot resolve 'opticalphoton' particle definition!");
      }
      G4ProcessManager* pmanager = particle->GetProcessManager();

      G4OpWLS* wls = new G4OpWLS();
      wls->UseTimeProfile("exponential");
      pmanager->AddDiscreteProcess(wls);

      info("+++ Geant4SciFiWLSSetup: G4OpWLS constructed and added directly (time profile: exponential)");
    }
  };

}
}

#include <DDG4/Factories.h>
using namespace dd4hep::sim;
DECLARE_GEANT4ACTION(Geant4SciFiWLSSetup)
