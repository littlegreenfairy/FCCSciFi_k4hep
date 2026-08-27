// SciFi stepping action for DD4hep/DDG4
// Tracks optical photon boundary processes, fiber exit photons,
// and primary particle steps. At end of run, merges stepping data
// as EDM4hep SimTrackerHit collections into the main ddsim output file.

#include "DDG4/Geant4SteppingAction.h"
#include "DDG4/Geant4RunAction.h"

#include <G4EventManager.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4OpticalPhoton.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessVector.hh>
#include <G4Run.hh>
#include <G4RunManager.hh>
#include <G4Step.hh>
#include <G4SteppingManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VProcess.hh>

#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include "podio/ROOTWriter.h"
#include "edm4hep/SimTrackerHitCollection.h"

#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace dd4hep {
namespace sim {

  class SciFiSteppingAction : public Geant4SteppingAction {
  public:
    SciFiSteppingAction(Geant4Context* ctxt, const std::string& nam)
        : Geant4SteppingAction(ctxt, nam) {
      declareProperty("OutputFile", m_outputFileName);
      declareProperty("KillEscapedPhotons", m_killEscapedPhotons);
    }

    ~SciFiSteppingAction() override = default;

    void operator()(const G4Step* step, G4SteppingManager* /*mgr*/) override {
      if (!m_endRunRegistered) {
        context()->runAction().callAtEnd(this, &SciFiSteppingAction::endRun);
        m_endRunRegistered = true;
      }

      int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
      if (eventID != m_currentEventID) {
        if (m_currentEventID >= 0) printEventSummary();
        beginEvent(eventID);
      }

      if (step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
        processOpticalPhoton(step);
      } else if (step->GetTrack()->GetParentID() == 0) {
        processPrimary(step);
      }
    }

  private:

    // ── data structs ────────────────────────────────────────────────

    struct PhotonData {
      double x, y, z, px, py, pz;
      double energy_eV, time, pathLength;
      int trackId, parentId, creatorProcess;
      int reflTotalCoreClad, reflTotalCladClad;
      int reflFresnelCoreClad, reflFresnelCladClad;
      int rayleigh;
      float lengthInCore, lengthInClad1, lengthInClad2;
    };

    struct StepData {
      double x, y, z, px, py, pz;
      double edep, stepLength, time;
    };

    struct EventData {
      std::vector<PhotonData> photons;
      std::vector<StepData>   primaries;
    };

    // ── volume / process helpers ────────────────────────────────────

    G4String getPreLVName(const G4Step* step) const {
      auto* vol = step->GetPreStepPoint()->GetPhysicalVolume();
      if (!vol) return "";
      return vol->GetLogicalVolume()->GetName();
    }

    G4String getPostLVName(const G4Step* step) const {
      if (step->GetTrack()->GetTrackStatus() == fStopAndKill) return "";
      auto* vol = step->GetPostStepPoint()->GetPhysicalVolume();
      if (!vol) return "";
      return vol->GetLogicalVolume()->GetName();
    }

    int getCreatorProcessId(const G4Step* step) const {
      auto* proc = step->GetTrack()->GetCreatorProcess();
      if (!proc) return 0;
      const G4String& name = proc->GetProcessName();
      if (name == "Scintillation" || name == "ScintillationPhys") return 1;
      if (name == "Cerenkov" || name == "CerenkovPhys") return 2;
      if (name == "OpWLS") return 3;
      return 0;
    }

    static bool contains(const G4String& str, const char* sub) {
      return str.find(sub) != std::string::npos;
    }

    G4OpBoundaryProcess* findBoundaryProcess() {
      if (m_opBoundaryProcess) return m_opBoundaryProcess;
      G4ProcessManager* mgr = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
      if (!mgr) return nullptr;
      G4ProcessVector* pv = mgr->GetPostStepProcessVector(typeDoIt);
      for (std::size_t i = 0; i < pv->entries(); ++i) {
        m_opBoundaryProcess = dynamic_cast<G4OpBoundaryProcess*>((*pv)[i]);
        if (m_opBoundaryProcess) return m_opBoundaryProcess;
      }
      return nullptr;
    }

    // ── per-step processing ─────────────────────────────────────────

    void processOpticalPhoton(const G4Step* step) {
      G4int trackId    = step->GetTrack()->GetTrackID();
      G4double stepLen = step->GetStepLength();

      if (m_seenPhotonTracks.insert(trackId).second) {
        auto* proc = step->GetTrack()->GetCreatorProcess();
        G4String procName = proc ? proc->GetProcessName() : "none";
        bool isScint = (procName == "Scintillation" || procName == "ScintillationPhys");
        bool isCerenkov = (procName == "Cerenkov" || procName == "CerenkovPhys");
        if (isScint) m_nScintPhotons++;
        else if (isCerenkov) m_nCerenkovPhotons++;
        else {
          m_nOtherPhotons++;
          if (m_firstUnknownProcess.empty()) m_firstUnknownProcess = procName;
        }
      }

      G4String preLV  = getPreLVName(step);
      G4String postLV = getPostLVName(step);

      if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {

        if (postLV == "lvSciFiMatSlab2" && preLV == "lvSciFiFiberClad2") {
          if (stepLen > 0) m_reflSurf[trackId]++;
        }

        G4OpBoundaryProcessStatus boundaryStatus = Undefined;
        G4OpBoundaryProcess*      bp             = findBoundaryProcess();
        if (bp) boundaryStatus = bp->GetStatus();

        if (boundaryStatus != Undefined) {
          if ((preLV == "lvSciFiFiberClad1" && postLV == "lvSciFiFiberClad2") ||
              (preLV == "lvSciFiFiberClad2" && postLV == "lvSciFiFiberClad1")) {
            if (boundaryStatus == TotalInternalReflection) m_reflTotalCladClad[trackId]++;
            if (boundaryStatus == FresnelReflection) m_reflFresnelCladClad[trackId]++;
            if (boundaryStatus == FresnelRefraction) m_refracCladClad[trackId]++;
          }
          if ((preLV == "lvSciFiFiberCore" && postLV == "lvSciFiFiberClad1") ||
              (preLV == "lvSciFiFiberClad1" && postLV == "lvSciFiFiberCore")) {
            if (boundaryStatus == TotalInternalReflection) m_reflTotalCoreClad[trackId]++;
            if (boundaryStatus == FresnelReflection) m_reflFresnelCoreClad[trackId]++;
            if (boundaryStatus == FresnelRefraction) m_refracCoreClad[trackId]++;
          }
        }

        bool preIsSiPM  = contains(preLV, "SiPM") || contains(preLV, "Pixel");
        bool postIsSiPM = contains(postLV, "SiPM") || contains(postLV, "Pixel");
        if (postIsSiPM && !preIsSiPM) {
          storeFiberExitPhoton(step, trackId);
        }

        if (m_killEscapedPhotons &&
            (preLV == "lvSciFiFiberClad2" || preLV == "lvSciFiFiberClad1") &&
            postLV == "lvSciFiMatSlab2") {
          step->GetTrack()->SetTrackStatus(fStopAndKill);
        }

      } else {
        G4ThreeVector dMom = step->GetDeltaMomentum();
        if (dMom[0] != 0 || dMom[1] != 0 || dMom[2] != 0) {
          m_rayleigh[trackId]++;
        }
      }

      if (preLV == "lvSciFiFiberCore") m_lengthCore[trackId] += stepLen;
      if (preLV == "lvSciFiFiberClad1") m_lengthClad1[trackId] += stepLen;
      if (preLV == "lvSciFiFiberClad2") m_lengthClad2[trackId] += stepLen;
    }

    void storeFiberExitPhoton(const G4Step* step, G4int trackId) {
      PhotonData d;
      d.x  = step->GetTrack()->GetPosition().x();
      d.y  = step->GetTrack()->GetPosition().y();
      d.z  = step->GetTrack()->GetPosition().z();
      d.px = step->GetTrack()->GetMomentumDirection().x();
      d.py = step->GetTrack()->GetMomentumDirection().y();
      d.pz = step->GetTrack()->GetMomentumDirection().z();
      d.energy_eV     = step->GetTrack()->GetKineticEnergy() / CLHEP::eV;
      d.time          = step->GetTrack()->GetGlobalTime();
      d.pathLength    = step->GetTrack()->GetTrackLength();
      d.trackId       = trackId;
      d.parentId      = step->GetTrack()->GetParentID();
      d.creatorProcess = getCreatorProcessId(step);
      d.reflTotalCoreClad   = m_reflTotalCoreClad[trackId];
      d.reflTotalCladClad   = m_reflTotalCladClad[trackId];
      d.reflFresnelCoreClad = m_reflFresnelCoreClad[trackId];
      d.reflFresnelCladClad = m_reflFresnelCladClad[trackId];
      d.rayleigh            = m_rayleigh[trackId];
      d.lengthInCore  = m_lengthCore[trackId];
      d.lengthInClad1 = m_lengthClad1[trackId];
      d.lengthInClad2 = m_lengthClad2[trackId];
      m_allEvents[m_currentEventID].photons.push_back(d);
    }

    void processPrimary(const G4Step* step) {
      StepData d;
      d.x  = step->GetTrack()->GetPosition().x();
      d.y  = step->GetTrack()->GetPosition().y();
      d.z  = step->GetTrack()->GetPosition().z();
      d.px = step->GetPreStepPoint()->GetMomentum().x() / CLHEP::GeV;
      d.py = step->GetPreStepPoint()->GetMomentum().y() / CLHEP::GeV;
      d.pz = step->GetPreStepPoint()->GetMomentum().z() / CLHEP::GeV;
      d.edep       = step->GetTotalEnergyDeposit() / CLHEP::GeV;
      d.stepLength = step->GetStepLength();
      d.time       = step->GetTrack()->GetGlobalTime();
      m_allEvents[m_currentEventID].primaries.push_back(d);
    }

    // ── event / run lifecycle ───────────────────────────────────────

    void printEventSummary() {
      size_t nExit = m_allEvents.count(m_currentEventID)
                         ? m_allEvents[m_currentEventID].photons.size()
                         : 0;
      info("+++ Event %d: %d scintillation, %d Cerenkov, %d other (first='%s') photons, %zu total, %zu fiber exits",
           m_currentEventID, m_nScintPhotons, m_nCerenkovPhotons,
           m_nOtherPhotons, m_firstUnknownProcess.c_str(),
           m_seenPhotonTracks.size(), nExit);
    }

    void beginEvent(int eventID) {
      m_currentEventID = eventID;
      m_seenPhotonTracks.clear();
      m_nScintPhotons = 0;
      m_nCerenkovPhotons = 0;
      m_nOtherPhotons = 0;
      m_firstUnknownProcess.clear();
      m_reflSurf.clear();
      m_reflTotalCladClad.clear();
      m_reflTotalCoreClad.clear();
      m_reflFresnelCladClad.clear();
      m_reflFresnelCoreClad.clear();
      m_refracCladClad.clear();
      m_refracCoreClad.clear();
      m_rayleigh.clear();
      m_lengthCore.clear();
      m_lengthClad1.clear();
      m_lengthClad2.clear();
    }

    void endRun(const G4Run*) {
      if (m_currentEventID >= 0) printEventSummary();
      if (m_allEvents.empty() || m_outputFileName.empty()) return;
      mergeIntoMainOutput();
    }

    // ── merge stepping data into main ddsim output ──────────────────

    void mergeIntoMainOutput() {
      podio::ROOTReader reader;
      reader.openFile(m_outputFileName);

      std::string tmpFile = m_outputFileName + ".stepping_tmp";
      podio::ROOTWriter writer(tmpFile);

      auto categories = reader.getAvailableCategories();
      for (const auto& catView : categories) {
        std::string cat(catView);
        unsigned nEntries = reader.getEntries(cat);
        for (unsigned i = 0; i < nEntries; i++) {
          auto frame = podio::Frame(reader.readEntry(cat, i));

          if (cat == "events") {
            int evtId = static_cast<int>(i);
            if (m_allEvents.count(evtId)) {
              addSteppingCollections(frame, m_allEvents[evtId]);
            }
          }

          writer.writeFrame(frame, cat);
        }
      }

      writer.finish();

      std::remove(m_outputFileName.c_str());
      std::rename(tmpFile.c_str(), m_outputFileName.c_str());

      m_allEvents.clear();
    }

    void addSteppingCollections(podio::Frame& frame, EventData& evtData) {
      auto photonColl = edm4hep::SimTrackerHitCollection();
      for (auto& d : evtData.photons) {
        auto hit = photonColl.create();
        hit.setPosition({d.x, d.y, d.z});
        hit.setMomentum({(float)d.px, (float)d.py, (float)d.pz});
        hit.setEDep(d.energy_eV * 1e-9f);
        hit.setTime(d.time);
        hit.setPathLength(d.pathLength);
        hit.setQuality(d.creatorProcess);
        hit.setCellID(
            (static_cast<uint64_t>(d.trackId) & 0xFFFFFFFF) |
            (static_cast<uint64_t>(d.parentId) << 32));
      }
      frame.put(std::move(photonColl), "FiberExitPhotons");

      auto primaryColl = edm4hep::SimTrackerHitCollection();
      for (auto& d : evtData.primaries) {
        auto hit = primaryColl.create();
        hit.setPosition({d.x, d.y, d.z});
        hit.setMomentum({(float)d.px, (float)d.py, (float)d.pz});
        hit.setEDep(d.edep);
        hit.setTime(d.time);
        hit.setPathLength(d.stepLength);
      }
      frame.put(std::move(primaryColl), "PrimarySteps");
    }

    // ── member data ─────────────────────────────────────────────────

    std::string m_outputFileName;
    bool m_killEscapedPhotons{false};

    bool m_endRunRegistered{false};
    int  m_currentEventID{-1};
    G4OpBoundaryProcess* m_opBoundaryProcess{nullptr};

    std::map<int, EventData> m_allEvents;

    std::set<int> m_seenPhotonTracks;
    int m_nScintPhotons{0};
    int m_nCerenkovPhotons{0};
    int m_nOtherPhotons{0};
    std::string m_firstUnknownProcess;

    std::map<int, int> m_reflSurf;
    std::map<int, int> m_reflTotalCladClad;
    std::map<int, int> m_reflTotalCoreClad;
    std::map<int, int> m_reflFresnelCladClad;
    std::map<int, int> m_reflFresnelCoreClad;
    std::map<int, int> m_refracCladClad;
    std::map<int, int> m_refracCoreClad;
    std::map<int, int> m_rayleigh;
    std::map<int, float> m_lengthCore;
    std::map<int, float> m_lengthClad1;
    std::map<int, float> m_lengthClad2;
  };

} // namespace sim
} // namespace dd4hep

using dd4hep::sim::SciFiSteppingAction;
#include "DDG4/Factories.h"
DECLARE_GEANT4ACTION(SciFiSteppingAction)
