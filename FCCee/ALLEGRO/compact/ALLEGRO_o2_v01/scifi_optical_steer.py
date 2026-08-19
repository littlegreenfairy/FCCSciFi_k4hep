from DDSim.DD4hepSimulation import DD4hepSimulation

SIM = DD4hepSimulation()

# Building a modular physics list
SIM.physics.list = None

SIM.physics.rangecut = 1.0 #like in SciFiMatG4's PhysicsList.cc 

# Registers Cerenkov + Scintillation + optical-photon transport physics, together with the standard physics lists
# Adapted from k4geo/example/arcfullsim.py's setupCerenkov()

def setupOpticalPhysics(kernel):
    from DDG4 import PhysicsList

    seq = kernel.physicsList()

    em = PhysicsList(kernel, "G4EmStandardPhysics/EmStandardPhys")
    em.enableUI()
    seq.adopt(em)

    cerenkov = PhysicsList(kernel, "Geant4CerenkovPhysics/CerenkovPhys")
    cerenkov.MaxNumPhotonsPerStep = 100
    cerenkov.MaxBetaChangePerStep = 10.0
    cerenkov.VerboseLevel = 0
    cerenkov.enableUI()
    seq.adopt(cerenkov)

    scintillation = PhysicsList(kernel, "Geant4ScintillationPhysics/ScintillationPhys")
    scintillation.ScintByParticleType = False
    scintillation.VerboseLevel = 0
    scintillation.enableUI()
    seq.adopt(scintillation)

    ph = PhysicsList(kernel, "Geant4OpticalPhotonPhysics/OpticalGammaPhys")
    ph.addParticleConstructor("G4OpticalPhoton")
    ph.WLSTimeProfile = "exponential"
    ph.VerboseLevel = 0
    ph.BoundaryInvokeSD = False 
    #Disables triggering the sensitive-detector hit at the geometry-boundary step
    #(fix inherited from ARC simulation)
    ph.enableUI()
    seq.adopt(ph)

    hadronic = PhysicsList(kernel, "G4HadronPhysicsQGSP_BERT/HadronPhys")
    hadronic.enableUI()
    seq.adopt(hadronic)

    return None


SIM.physics.setupUserPhysics(setupOpticalPhysics)

# Use DD4hep optical tracker like for ARC (SciFiTracker is the detector name in the xml file)
SIM.action.mapActions["SciFiTracker"] = "Geant4OpticalTrackerAction"
SIM.filter.mapDetFilter["SciFiTracker"] = None
