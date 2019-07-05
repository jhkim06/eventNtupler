import FWCore.ParameterSet.Config as cms
process = cms.Process("ntupler")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("./ISR_GENtuple_aMCNLO.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",

    fileNames = cms.untracked.vstring(
    "file:./HIG-RunIISummer15wmLHEGS-00152_7.root"    
    )
)

process.generatorSmeared = cms.EDProducer('GeneratorSmearedProducer',
  currentTag = cms.untracked.InputTag('VtxSmeared'),
  previousTag = cms.untracked.InputTag('generator')
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi") # need Pythia particle table in the EventSeup in GenParticleProducer
process.genParticlesNew = cms.EDProducer("GenParticleProducer",
    saveBarCodes = cms.untracked.bool(True),
    src = cms.InputTag("generatorSmeared"),
    abortOnUnknownPDGCode = cms.untracked.bool(False)
)
       
process.GEN = cms.EDAnalyzer('eventNtupler',
    genSrc = cms.InputTag("genParticlesNew")
)

#
process.p = cms.Path(
    process.generatorSmeared+
    process.genParticlesNew+
    process.GEN
)
