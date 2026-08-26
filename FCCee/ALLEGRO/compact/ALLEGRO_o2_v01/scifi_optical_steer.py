from DDSim.DD4hepSimulation import DD4hepSimulation

SIM = DD4hepSimulation()

SIM.physics.list = "QGSP_BERT"

SIM.physics.rangecut = 1.0 #like in SciFiMatG4's PhysicsList.cc 

# Registers Cerenkov + Scintillation + optical-photon transport physics, together with the standard physics lists
# Adapted from k4geo/example/arcfullsim.py's setupCerenkov()

def setupOpticalPhysics(kernel, _sim=SIM):
    from DDG4 import PhysicsList

    seq = kernel.physicsList()

    cerenkov = PhysicsList(kernel, "Geant4CerenkovPhysics/CerenkovPhys")
    cerenkov.MaxNumPhotonsPerStep = 100
    cerenkov.MaxBetaChangePerStep = 10.0
    cerenkov.VerboseLevel = 1
    cerenkov.enableUI()
    seq.adopt(cerenkov)

    scintillation = PhysicsList(kernel, "Geant4ScintillationPhysics/ScintillationPhys")
    # Property is named "ByParticleType" here, not "ScintByParticleType" --
    # confirmed via hasProperty() on the live object in this DD4hep release.
    scintillation.ByParticleType = False
    scintillation.VerboseLevel = 1
    scintillation.enableUI()
    seq.adopt(scintillation)

    ph = PhysicsList(kernel, "Geant4OpticalPhotonPhysics/OpticalGammaPhys")
    ph.addParticleConstructor("G4OpticalPhoton")
    ph.VerboseLevel = 1
    ph.BoundaryInvokeSD = False
    ph.enableUI()
    seq.adopt(ph)

    from DDG4 import SteppingAction
    scifi_step = SteppingAction(kernel, "SciFiSteppingAction/SciFiStep")
    scifi_step.OutputFile = _sim.outputFile
    scifi_step.KillEscapedPhotons = False
    scifi_step.enableUI()
    kernel.steppingAction().adopt(scifi_step)

    return None


SIM.physics.setupUserPhysics(setupOpticalPhysics)

# Use DD4hep optical tracker like for ARC (SciFiTracker is the detector name in the xml file)
SIM.action.mapActions["SciFiTracker"] = "Geant4SciFiOpticalTrackerAction"
SIM.filter.mapDetFilter["SciFiTracker"] = None

# Optical photons carry only a few eV, so DDSim's default particle-saving energy cut
# would silently drop them from the output MCParticles collection even when they are
# genuinely produced. saveProcesses alone did not surface them (process-name mismatch?),
# so force-save everything for now -- revisit with a narrower cut once confirmed working.
SIM.part.keepAllParticles = True
SIM.part.minimalKineticEnergy = "0*eV"

