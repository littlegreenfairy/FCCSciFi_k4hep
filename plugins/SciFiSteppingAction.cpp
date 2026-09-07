// SciFi stepping action for DD4hep/DDG4
// Tracks optical photon boundary processes, fiber exit photons,
// primary particle steps, and photon loss classification.
// At end of run, merges stepping data as EDM4hep SimTrackerHit
// collections into the main ddsim output file.

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

#include <G4Material.hh>
#include <G4MaterialPropertiesTable.hh>

#include "TFile.h"
#include "TH1F.h"

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
        dumpWLSDiagnostics(step);
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

    struct PhotonCounters {
      int nTotal{0};
      int nScint{0};
      int nCerenkov{0};
      int nWLSCreated{0};
      int nFiberExits{0};
      int nDetected{0};
      int nAbsWLS{0};
      int nAbsAtten{0};
      int nAbsClad{0};
      int nAbsSiPM{0};
      int nNonCaptured{0};
      int nAbsUnknown{0};
      int nLostSurface{0};
      int nLostEscaped{0};
      int nLostWorld{0};
      int nLostOther{0};

      void clear() { *this = PhotonCounters{}; }

      // Accumulates one event's counts into a run-level total. Field list
      // has to be kept in sync with the struct above -- there's no reflection
      // to do this generically here.
      PhotonCounters& operator+=(const PhotonCounters& o) {
        nTotal += o.nTotal; nScint += o.nScint; nCerenkov += o.nCerenkov;
        nWLSCreated += o.nWLSCreated; nFiberExits += o.nFiberExits; nDetected += o.nDetected;
        nAbsWLS += o.nAbsWLS; nAbsAtten += o.nAbsAtten; nAbsClad += o.nAbsClad;
        nAbsSiPM += o.nAbsSiPM; nNonCaptured += o.nNonCaptured; nAbsUnknown += o.nAbsUnknown;
        nLostSurface += o.nLostSurface; nLostEscaped += o.nLostEscaped;
        nLostWorld += o.nLostWorld; nLostOther += o.nLostOther;
        return *this;
      }
    };

    // creator process category for a track
    enum PhotonType { kPrimary = 0, kWLS = 1 };

    // ── diagnostics ─────────────────────────────────────────────────

    void dumpWLSDiagnostics(const G4Step* step) {
      // Check if OpWLS process exists on optical photons
      G4ProcessManager* pmgr = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
      bool foundWLS = false;
      if (pmgr) {
        G4ProcessVector* pv = pmgr->GetProcessList();
        info("=== Optical photon processes ===");
        for (G4int i = 0; i < pv->entries(); ++i) {
          info("  process[%d]: %s", i, (*pv)[i]->GetProcessName().c_str());
          if ((*pv)[i]->GetProcessName() == "OpWLS") foundWLS = true;
        }
      }
      info("  OpWLS registered: %s", foundWLS ? "YES" : "NO");

      // Check core material for WLS properties
      auto* preVol = step->GetPreStepPoint()->GetPhysicalVolume();
      if (preVol) {
        G4Material* mat = preVol->GetLogicalVolume()->GetMaterial();
        info("=== Material: %s ===", mat->GetName().c_str());
        G4MaterialPropertiesTable* mpt = mat->GetMaterialPropertiesTable();
        if (mpt) {
          auto* wlsComp = mpt->GetProperty("WLSCOMPONENT");
          auto* wlsAbs  = mpt->GetProperty("WLSABSLENGTH");
          bool  hasTime = mpt->ConstPropertyExists("WLSTIMECONSTANT");
          info("  WLSCOMPONENT:    %s (entries: %zu)",
               wlsComp ? "YES" : "NO", wlsComp ? wlsComp->GetVectorLength() : 0);
          info("  WLSABSLENGTH:    %s (entries: %zu)",
               wlsAbs ? "YES" : "NO", wlsAbs ? wlsAbs->GetVectorLength() : 0);
          info("  WLSTIMECONSTANT: %s", hasTime ? "YES" : "NO");
          if (hasTime)
            info("  WLSTIMECONSTANT value: %g ns", mpt->GetConstProperty("WLSTIMECONSTANT") / CLHEP::ns);
          if (wlsAbs && wlsAbs->GetVectorLength() > 0) {
            info("  WLSABSLENGTH at 3.0 eV: %g mm", wlsAbs->Value(3.0 * CLHEP::eV) / CLHEP::mm);
            info("  WLSABSLENGTH at 3.5 eV: %g mm", wlsAbs->Value(3.5 * CLHEP::eV) / CLHEP::mm);
          }
        } else {
          info("  No MaterialPropertiesTable!");
        }
      }
    }

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

    bool isWLSPhoton(const G4Step* step) const {
      auto* proc = step->GetTrack()->GetCreatorProcess();
      if (!proc) return false;
      return proc->GetProcessName() == "OpWLS";
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

      bool wls = isWLSPhoton(step);

      if (m_seenPhotonTracks.insert(trackId).second) {
        auto* proc = step->GetTrack()->GetCreatorProcess();
        G4String procName = proc ? proc->GetProcessName() : "none";
        bool isScint = (procName == "Scintillation" || procName == "ScintillationPhys");
        bool isCerenkov = (procName == "Cerenkov" || procName == "CerenkovPhys");

        if (wls) {
          m_cWLS.nTotal++;
          m_cWLS.nWLSCreated++;
          m_cAll.nWLSCreated++;
        } else {
          m_cAll.nTotal++;
          m_cPrim.nTotal++;
          if (isScint) { m_cAll.nScint++; m_cPrim.nScint++; }
          else if (isCerenkov) { m_cAll.nCerenkov++; m_cPrim.nCerenkov++; }
        }
      }

      G4String preLV  = getPreLVName(step);

      bool killedByUs = false;

      if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {

        G4String postLV = getPostLVName(step);

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

        // Track photons escaping fiber: cladding → material
        if ((preLV == "lvSciFiFiberClad2" || preLV == "lvSciFiFiberClad1") &&
            postLV == "lvSciFiMatSlab2") {
          m_escapedFiber.insert(trackId);
        }

        bool preIsSiPM  = contains(preLV, "SiPM") || contains(preLV, "Pixel");
        bool postIsSiPM = contains(postLV, "SiPM") || contains(postLV, "Pixel");
        if (postIsSiPM && !preIsSiPM) {
          storeFiberExitPhoton(step, trackId);
          m_fiberExitTracks.insert(trackId);
          m_cAll.nFiberExits++;
          if (wls) m_cWLS.nFiberExits++;
          else m_cPrim.nFiberExits++;
        }

        // SiPM fill factor: a photon landing on the inactive frame around a pixel is
        // absorbed there. Without this, the dead shell (lvSciFiPixel) is optically
        // identical to the sensitive core (both SciFiEpoxy, so there is no refractive
        // index step between them) and photons that miss the active area simply cross
        // into it from the side, so pixelSensitiveRatio never takes effect. Exact name
        // match, so lvSciFiPixelSensitive is not caught. Mirrors the DeadDetector kill
        // in SciFiMatG4's SteppingAction.cc. Placed after the fiber-exit registration
        // above so the photon still counts as having reached the SiPM.
        if (postLV == "lvSciFiPixel") {
          step->GetTrack()->SetTrackStatus(fStopAndKill);
          killedByUs = true;
        }

        if (m_killEscapedPhotons &&
            (preLV == "lvSciFiFiberClad2" || preLV == "lvSciFiFiberClad1") &&
            postLV == "lvSciFiMatSlab2") {
          step->GetTrack()->SetTrackStatus(fStopAndKill);
          killedByUs = true;
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

      if (step->GetTrack()->GetTrackStatus() == fStopAndKill) {
        trackPhotonDeath(step, trackId, killedByUs, wls);
      }
    }

    // ── photon loss classification ──────────────────────────────────

    void classifyLoss(const G4Step* step, G4int trackId, bool killedByUs,
                      const G4String& procName, const G4String& preLV,
                      PhotonCounters& c) {
      if (m_fiberExitTracks.count(trackId)) {
        c.nDetected++;
        return;
      }
      if (killedByUs) { c.nLostEscaped++; return; }

      auto* postVol = step->GetPostStepPoint()->GetPhysicalVolume();
      if (!postVol) { c.nLostWorld++; return; }

      if (procName == "OpAbsorption" || procName == "OpWLS") {
        if (preLV == "lvSciFiFiberCore") {
          if (procName == "OpWLS") c.nAbsWLS++;
          else c.nAbsAtten++;
        }
        else if (preLV == "lvSciFiFiberClad1" || preLV == "lvSciFiFiberClad2") c.nAbsClad++;
        else if (contains(preLV, "SiPM") || contains(preLV, "Pixel")) c.nAbsSiPM++;
        else if (m_escapedFiber.count(trackId)) c.nNonCaptured++;
        else c.nAbsUnknown++;
      } else if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
        c.nLostSurface++;
      } else {
        c.nLostOther++;
      }
    }

    void trackPhotonDeath(const G4Step* step, G4int trackId, bool killedByUs, bool wls) {
      auto* proc = step->GetPostStepPoint()->GetProcessDefinedStep();
      G4String procName = proc ? proc->GetProcessName() : "unknown";
      G4String preLV = getPreLVName(step);

      classifyLoss(step, trackId, killedByUs, procName, preLV, m_cAll);
      if (wls) classifyLoss(step, trackId, killedByUs, procName, preLV, m_cWLS);
      else     classifyLoss(step, trackId, killedByUs, procName, preLV, m_cPrim);

      if (m_cAll.nLostOther > 0 && m_cAll.nLostOther <= 3 &&
          !m_fiberExitTracks.count(trackId) && !killedByUs) {
        info("  unknown loss: proc='%s' vol='%s'", procName.c_str(), preLV.c_str());
      }
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

    void printCounters(const char* label, const PhotonCounters& c) {
      info("  [%s] %d photons (%d Cerenkov, %d scintillation, %d WLS-created) | %d fiber exits, %d detected",
           label, c.nTotal, c.nCerenkov, c.nScint, c.nWLSCreated, c.nFiberExits, c.nDetected);
      info("  [%s] losses: %d WLS-abs, %d atten-abs, %d abs-clad, %d abs-SiPM, %d non-captured, %d surface, %d escaped, %d left-world, %d unknown, %d other",
           label, c.nAbsWLS, c.nAbsAtten, c.nAbsClad, c.nAbsSiPM, c.nNonCaptured,
           c.nLostSurface, c.nLostEscaped, c.nLostWorld, c.nAbsUnknown, c.nLostOther);
    }

    void printEventSummary() {
      info("+++ Event %d summary:", m_currentEventID);
      printCounters("ALL", m_cAll);
      printCounters("WLS", m_cWLS);
      PhotonCounters cSurv = m_cPrim;
      cSurv.nTotal -= m_cAll.nWLSCreated;
      printCounters("Primary", cSurv);
    }

    // Folds the just-finished event's counters into the run-level totals
    // and records its WLS counters individually for the per-event
    // histograms. Shared by beginEvent() (every event but the first) and
    // endRun() (the last event, which no later beginEvent() catches).
    void foldCurrentEventIntoTotals() {
      m_runAll  += m_cAll;
      m_runWLS  += m_cWLS;
      m_runPrim += m_cPrim;
      m_perEventWLS.push_back(m_cWLS);
    }

    void beginEvent(int eventID) {
      // Fold the event just finished into the run-level totals before its
      // per-event counters get wiped below. Skipped on the very first call
      // (m_currentEventID still at its -1 default), matching the same
      // "nothing to summarize yet" guard endRun() already uses.
      if (m_currentEventID >= 0) foldCurrentEventIntoTotals();
      m_currentEventID = eventID;
      m_seenPhotonTracks.clear();
      m_fiberExitTracks.clear();
      m_escapedFiber.clear();
      m_cAll.clear();
      m_cWLS.clear();
      m_cPrim.clear();
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
      // Fold in the last event: beginEvent() won't run again to do it.
      if (m_currentEventID >= 0) {
        foldCurrentEventIntoTotals();
        printEventSummary();
      }
      if (m_outputFileName.empty()) return;
      // Kept separate from the histogram write below: the podio merge is
      // only worth its file-rewrite cost if there's actual per-photon data
      // to add, but the counter histograms are worth writing either way.
      if (!m_allEvents.empty()) mergeIntoMainOutput();
      writeSummaryHistograms();
    }

    // ── run-level summary histograms ────────────────────────────────

    // One entry per PhotonCounters field: a pointer-to-member plus its
    // display label. Single source of truth for both the aggregate
    // per-run histograms below and the per-event WLS histograms, so the
    // 16 destinies only need to be listed once.
    struct DestinyField { int PhotonCounters::*member; const char* label; };

    static const std::vector<DestinyField>& destinyFields() {
      static const std::vector<DestinyField> fields = {
        {&PhotonCounters::nTotal, "Total"}, {&PhotonCounters::nScint, "Scintillation"},
        {&PhotonCounters::nCerenkov, "Cerenkov"}, {&PhotonCounters::nWLSCreated, "WLSCreated"},
        {&PhotonCounters::nFiberExits, "FiberExits"}, {&PhotonCounters::nDetected, "Detected"},
        {&PhotonCounters::nAbsWLS, "AbsWLS"}, {&PhotonCounters::nAbsAtten, "AbsAtten"},
        {&PhotonCounters::nAbsClad, "AbsClad"}, {&PhotonCounters::nAbsSiPM, "AbsSiPM"},
        {&PhotonCounters::nNonCaptured, "NonCaptured"}, {&PhotonCounters::nAbsUnknown, "AbsUnknown"},
        {&PhotonCounters::nLostSurface, "LostSurface"}, {&PhotonCounters::nLostEscaped, "LostEscaped"},
        {&PhotonCounters::nLostWorld, "LostWorld"}, {&PhotonCounters::nLostOther, "LostOther"},
      };
      return fields;
    }

    // Run-total histogram: one bin per destiny, bin content = total count
    // across the whole run (i.e. what's already in a PhotonCounters).
    TH1F* makeCounterHistogram(const std::string& name, const std::string& title,
                                const PhotonCounters& c) const {
      const auto& fields = destinyFields();
      auto* h = new TH1F(name.c_str(), title.c_str(), (int)fields.size(), 0, (int)fields.size());
      for (size_t i = 0; i < fields.size(); ++i) {
        h->SetBinContent(i + 1, c.*(fields[i].member));
        h->GetXaxis()->SetBinLabel(i + 1, fields[i].label);
      }
      return h;
    }

    // Per-event distributions, WLS photons only: for each destiny, one
    // histogram whose entries are events -- each event contributes one
    // fill, at the number of WLS photons that met that destiny in that
    // event. Range is sized from the actual per-event maximum since it
    // isn't known ahead of time.
    void writePerEventWLSHistograms() const {
      for (const auto& field : destinyFields()) {
        int maxVal = 0;
        for (const auto& c : m_perEventWLS) {
          int v = c.*(field.member);
          if (v > maxVal) maxVal = v;
        }
        std::string name  = std::string("PerEventWLS_") + field.label;
        std::string title = std::string("WLS photons '") + field.label + "' per event;" +
                             field.label + " (WLS photons in event);number of events";
        auto* h = new TH1F(name.c_str(), title.c_str(), maxVal + 2, 0, maxVal + 2);
        for (const auto& c : m_perEventWLS) h->Fill(c.*(field.member));
        h->Write();
      }
    }

    void writeSummaryHistograms() {
      // mergeIntoMainOutput() above (when it ran) used podio's own
      // reader/writer, which owns and closes the file handle itself; by the
      // time we get here that's already finished, so reopening directly
      // with a plain TFile in UPDATE mode is safe.
      TFile* f = TFile::Open(m_outputFileName.c_str(), "UPDATE");
      if (!f || f->IsZombie()) {
        error("+++ Could not reopen %s to write summary histograms", m_outputFileName.c_str());
        return;
      }
      makeCounterHistogram("PhotonDestinies_All", "Photon destinies (all photons)", m_runAll)->Write();
      makeCounterHistogram("PhotonDestinies_WLS", "Photon destinies (WLS-created photons only)", m_runWLS)->Write();
      makeCounterHistogram("PhotonDestinies_Primary", "Photon destinies (primary-created photons only)", m_runPrim)->Write();
      writePerEventWLSHistograms();
      f->Close();
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

    // Photon tracking sets
    std::set<int> m_seenPhotonTracks;
    std::set<int> m_fiberExitTracks;
    std::set<int> m_escapedFiber;

    // Per-event counters split by photon origin (reset every event)
    PhotonCounters m_cAll;
    PhotonCounters m_cWLS;
    PhotonCounters m_cPrim;

    // Run-level totals: same breakdown, accumulated across every event,
    // never cleared. What the summary histograms are built from.
    PhotonCounters m_runAll;
    PhotonCounters m_runWLS;
    PhotonCounters m_runPrim;

    // One entry per event: that event's WLS-only counters, kept individually
    // (not just summed) so the per-event destiny histograms can be filled.
    std::vector<PhotonCounters> m_perEventWLS;

    // Per-event, per-track working maps
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
