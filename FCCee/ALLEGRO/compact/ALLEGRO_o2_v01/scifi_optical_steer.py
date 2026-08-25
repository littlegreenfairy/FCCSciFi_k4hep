from DDSim.DD4hepSimulation import DD4hepSimulation

SIM = DD4hepSimulation()

SIM.physics.list = "QGSP_BERT"

SIM.physics.rangecut = 1.0 #like in SciFiMatG4's PhysicsList.cc 

# Registers Cerenkov + Scintillation + optical-photon transport physics, together with the standard physics lists
# Adapted from k4geo/example/arcfullsim.py's setupCerenkov()

def setupOpticalPhysics(kernel):
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
    # WLSTimeProfile not implemented, we'll have to figure out how to set the time profile for wavelength-shifting photons later
    ph.VerboseLevel = 1
    ph.BoundaryInvokeSD = False 
    #Disables triggering the sensitive-detector hit at the geometry-boundary step
    #(fix inherited from ARC simulation)
    ph.enableUI()
    seq.adopt(ph)

    # Geometric step-by-step trace
    from DDG4 import SteppingAction
    step_trace = SteppingAction(kernel, "Geant4TestStepAction/StepTrace")
    step_trace.OutputLevel = 1
    step_trace.enableUI()
    kernel.steppingAction().adopt(step_trace)

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

