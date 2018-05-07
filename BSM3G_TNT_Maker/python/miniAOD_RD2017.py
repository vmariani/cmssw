import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import copy
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

options = VarParsing.VarParsing('analysis')
# Variables one can control from the multicrab configuration file.
# When connecting a variable you need to tell the module certain information
# about the object.
#                   - Object name.
#                   - Default value.
#                   - Is object a single number or a list.
#                   - Object type.
#                   - Details of object.
#

# ===== Register new variables =====
options.register('optionJECAK4PFchsDATA1',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L1FastJet_AK4PFchs.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFchsDATA1 JEC file")

options.register('optionJECAK4PFchsDATA2',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2Relative_AK4PFchs.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFchsDATA2 JEC file"
)
options.register('optionJECAK4PFchsDATA3',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L3Absolute_AK4PFchs.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFchsDATA3 JEC file"
)
options.register('optionJECAK4PFchsDATA4',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2L3Residual_AK4PFchs.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFchsDATA4 JEC file"
)
options.register('optionJECAK4PFchsDATAUnc',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_Uncertainty_AK4PFchs.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFchsDATAUnc JEC file"
)


options.register('optionJECAK4PFPuppiDATA1',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L1FastJet_AK4PFPuppi.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFPuppiDATA1 JEC file"
)
options.register('optionJECAK4PFPuppiDATA2',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2Relative_AK4PFPuppi.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFPuppiDATA2 JEC file"
)
options.register('optionJECAK4PFPuppiDATA3',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L3Absolute_AK4PFPuppi.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFPuppiDATA3 JEC file"
)
options.register('optionJECAK4PFPuppiDATA4',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2L3Residual_AK4PFPuppi.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFPuppiDATA4 JEC file"
)
options.register('optionJECAK4PFPuppiDATAUnc',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_Uncertainty_AK4PFPuppi.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK4PFPuppiDATAUnc JEC file"
)



options.register('optionJECAK8PFchsDATA1',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L1FastJet_AK8PFchs.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK8PFchsDATA1 JEC file"
)
options.register('optionJECAK8PFchsDATA2',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2Relative_AK8PFchs.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK8PFchsDATA2 JEC file"
)
options.register('optionJECAK8PFchsDATA3',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L3Absolute_AK8PFchs.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK8PFchsDATA3 JEC file"
)
options.register('optionJECAK8PFchsDATA4',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2L3Residual_AK8PFchs.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK8PFchsDATA4 JEC file"
)
options.register('optionJECAK8PFchsDATAUnc',
'BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_Uncertainty_AK8PFchs.txt',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for AK8PFchsDATAUnc JEC file"
)

options.register('ofName',
'sentinel_output_name',
VarParsing.VarParsing.multiplicity.singleton,
VarParsing.VarParsing.varType.string,
"Name for output file."
)


# ===== Get & parse any command line arguments =====
options.parseArguments()


#####
##   Initial standard configs
#####
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag.globaltag = '94X_dataRun2_ReReco17_forValidation'
process.GlobalTag.globaltag = '94X_dataRun2_v6'
process.prefer("GlobalTag")
process.load('Configuration.StandardSequences.Services_cff')

#####
##   Input files
#####
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    '/store/data/Run2017B/SingleMuon/MINIAOD/17Nov2017-v1/40000/0021369B-9BD8-E711-BFE9-FA163EAA42CB.root'
  ),
  skipEvents = cms.untracked.uint32(0)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

##### JEC
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = cms.InputTag('slimmedJets'),
  labelName = 'UpdatedJEC',
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet','L2Relative','L3Absolute']), 'None')
)
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
#####
##   ELECTRON ID SECTION
#####
process.load("RecoEgamma/ElectronIdentification/ElectronIDValueMapProducer_cfi")
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
switchOnVIDPhotonIdProducer(process,DataFormat.MiniAOD)
my_id_modules = [ 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                  'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
                  'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff', 
                  'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
                  'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff'
                  ]
pho_id_modules =[ 'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V1_TrueVtx_cff',
                  'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff',  
                  ]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
for idmod in pho_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

######
### Electron smear and regression
######
process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
process.load('EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   slimmedElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                                 engineName = cms.untracked.string('TRandom3'),
                                                                                 ),
                                                   slimmedPhotons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                                               engineName = cms.untracked.string('TRandom3'),
                                                                               ),
                                                   )
process.slimmedElectrons = process.calibratedPatElectrons.clone(electrons=cms.InputTag("slimmedElectrons",processName=cms.InputTag.skipCurrentProcess()))
process.slimmedPhotons = process.calibratedPatPhotons.clone(photons=cms.InputTag("slimmedPhotons",processName=cms.InputTag.skipCurrentProcess())) 
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
process.photonMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedPhotons')
process.photonIDValueMapProducer.srcMiniAOD = cms.InputTag('slimmedPhotons')
process.egammaScaleSmearTask = cms.Task(process.slimmedElectrons,process.slimmedPhotons)
process.egammaScaleSmearSeq = cms.Sequence( process.egammaScaleSmearTask)
process.egammaScaleSmearAndVIDSeq = cms.Sequence(process.egammaScaleSmearSeq*process.egmGsfElectronIDSequence*process.egmPhotonIDSequence)

#####
##   For tt+X
#####
# Setting input particle collections to be used by the tools
genParticleCollection = "prunedGenParticles"
genJetCollection      = "slimmedGenJets"
jetFlavourInfos       = "genJetFlavourInfos"
jetAlgo               = "AntiKt"
rParam                = 0.4
genJetPtMin           = 20.
genJetAbsEtaMax       = 2.4
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenCHadron
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(particles = genParticleCollection)
process.genJetFlavourInfos = ak4JetFlavourInfos.clone(jets=genJetCollection,rParam=cms.double(rParam),jetAlgorithm = cms.string(jetAlgo))
process.matchGenBHadron = matchGenBHadron.clone(genParticles = genParticleCollection,jetFlavourInfos = jetFlavourInfos)
process.matchGenCHadron = matchGenCHadron.clone(genParticles = genParticleCollection,jetFlavourInfos = jetFlavourInfos)


#####Tau#####

from BSMFramework.BSM3G_TNT_Maker.runTauIdMVA import *
na = TauIDEmbedder(process, cms,
        debug=True,
        toKeep = ["dR0p32017v2"] # pick the one you need: ["2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1"]
        )
na.runTauID()


############### MET Re-correct ##################
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD( process, True )
#PFMet
runMetCorAndUncFromMiniAOD(process,
                           isData=True
                           )
#Puppi MET
runMetCorAndUncFromMiniAOD(process,
                           isData=True,
                           metType="Puppi",
                           postfix="Puppi",
                           jetFlavor="AK4PFPuppi",
                           )

process.puppiNoLep.useExistingWeights = False
process.puppi.useExistingWeights = False



#####
##   Output file
#####
options.ofName += ".root"
process.TFileService = cms.Service("TFileService",
  fileName = cms.string("OutTree.root")
)

#####
##   Analysis parameters
#####
process.TNT = cms.EDAnalyzer("BSM3G_TNT_Maker",
  #### Running options
  # Choose which trigger you want (do NOT need to put * as it will consider all the versions by default)
  ifevtriggers      = cms.bool(True), # True means you want to require the triggers
  maxtriggerversion = cms.double(20), # please leave it as a double
  evtriggers        = cms.vstring(
  #############################
     #'HLT_Ele115_CaloIdVT_GsfTrkIdT_v',
#     'HLT_DoubleEle33_CaloIdL_MW_v',
     #'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v',
     #'HLT_IsoMu24_v',
     #'HLT_IsoTkMu24_v',
#     'HLT_Mu50_v',
#     'HLT_TkMu50_v',
#     'HLT_Mu30_TkMu11_v',
    'HLT_Ele32_WPTight_Gsf_v',
    'HLT_Ele35_WPTight_Gsf_v',
    'HLT_IsoMu24_v',
    'HLT_IsoMu27_v',
    'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v',
    'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
    'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
    'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
    'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
    'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
    'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v',
    'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v',
    'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v',
    'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v',
    'HLT_TripleMu_12_10_5_v',
  ),
  # Choose which information you want to use
  fillgeninfo           = cms.bool(False),
  fillgenHFCategoryinfo = cms.bool(False),
  filleventinfo         = cms.bool(True),
  filltriggerinfo       = cms.bool(True),
  fillPVinfo            = cms.bool(True),
  fillmuoninfo          = cms.bool(True),
  fillelectronpatinfo   = cms.bool(True),
  filltauinfo           = cms.bool(True),
  filljetinfo           = cms.bool(True),
  filltthjetinfo        = cms.bool(False), #F
  fillBoostedJetinfo    = cms.bool(True),
  fillTopSubJetinfo     = cms.bool(False), #F
  fillTauJetnessinfo    = cms.bool(False),
  fillBJetnessinfo      = cms.bool(False),
  fillBJetnessFVinfo    = cms.bool(False),
  fillBTagReweight      = cms.bool(False),
  fillPileupReweight    = cms.bool(True),
  fillMETinfo           = cms.bool(True),
  fillphotoninfo        = cms.bool(False), #F   
  # Choose format 
  MiniAODv2 = cms.bool(True),
  is_data   = cms.bool(True),
  reHLT     = cms.bool(True),
  debug_    = cms.bool(False),
  super_TNT = cms.bool(False),
  AJVar     = cms.bool(False),
  tthlepVar = cms.bool(True),
  bjetnessselfilter = cms.bool(False),
  PuppiVar  = cms.bool(False),
  qglVar    = cms.bool(True),
  MC2016    = cms.bool(False),
  # Input tags 
  bits                = cms.InputTag("TriggerResults","","HLT"),
  prescales           = cms.InputTag("patTrigger"),
  objects             = cms.InputTag("selectedPatTrigger"),  
  vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
  beamSpot            = cms.InputTag("offlineBeamSpot"),
  muons               = cms.InputTag("slimmedMuons"),
  patElectrons        = cms.InputTag("slimmedElectrons"),
  electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-veto"),
  electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose"),
  electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium"),
  electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight"),
  eleMVATrigIdMap        = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp80"),
  eleMVAnonTrigIdMap     = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
  eleMVATrigwp90IdMap    = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp90"),
  eleMVAnonTrigwp90IdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
  eleMVATrigwpLooseIdMap    = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wpLoose"),
  eleMVAnonTrigwpLooseIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose"),
  eleHEEPIdMap                 = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
  elemvaValuesMap_Trig      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
  elemvaCategoriesMap_Trig  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
  elemvaValuesMap_nonTrig         = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
  elemvaCategoriesMap_nonTrig     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
  eleMVAHZZwpLooseIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-HZZ-V1-wpLoose"),
  elemvaValuesMap_HZZ          = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
  elemvaCategoriesMap_HZZ      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Categories"),
  ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),
  #taus                = cms.InputTag("slimmedTaus"),
  taus                = cms.InputTag("NewTauIDsEmbedded"),
  jets                = cms.InputTag("slimmedJets"),
  lepjets             = cms.InputTag("updatedPatJetsUpdatedJEC"),
  jetsPUPPI           = cms.InputTag("slimmedJetsPuppi"),
  fatjets             = cms.InputTag("slimmedJetsAK8"),
  topsubjets          = cms.InputTag("slimmedJetsCMSTopTagCHSPacked", "SubJets"),
  mets                = cms.InputTag("slimmedMETs"),
  metsPUPPI           = cms.InputTag("slimmedMETsPuppi"),
  metFilterBits       = cms.InputTag("TriggerResults", "", "RECO"),
  photons             = cms.InputTag("slimmedPhotons"),
  packedPFCandidates  = cms.InputTag("packedPFCandidates"), 
  pruned              = cms.InputTag("prunedGenParticles"),
  # JER
  jerAK4PFchs     =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt"),
  jerAK4PFchsSF   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK4PFchs.txt"),
  jerAK4PFPuppi   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt"),
  jerAK4PFPuppiSF =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt"),
  jerAK8PFchs     =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt"),
  jerAK8PFchsSF   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK8PFchs.txt"),
  jerAK8PFPuppi   =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt"),
  jerAK8PFPuppiSF =  cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JER/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt"),
  # JEC - CORRECTIONS ON FLY
  jecPayloadNamesAK4PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFPuppiMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_Uncertainty_AK8PFchs.txt"),
   #=== DATA ===
   jecPayloadNamesAK4PFchsDATA1   = cms.FileInPath(options.optionJECAK4PFchsDATA1),
   jecPayloadNamesAK4PFchsDATA2   = cms.FileInPath(options.optionJECAK4PFchsDATA2),
   jecPayloadNamesAK4PFchsDATA3   = cms.FileInPath(options.optionJECAK4PFchsDATA3),
   jecPayloadNamesAK4PFchsDATA4   = cms.FileInPath(options.optionJECAK4PFchsDATA4),
   jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath(options.optionJECAK4PFchsDATAUnc),
   jecPayloadNamesAK4PFPuppiDATA1   = cms.FileInPath(options.optionJECAK4PFPuppiDATA1),
   jecPayloadNamesAK4PFPuppiDATA2   = cms.FileInPath(options.optionJECAK4PFPuppiDATA2),
   jecPayloadNamesAK4PFPuppiDATA3   = cms.FileInPath(options.optionJECAK4PFPuppiDATA3),
   jecPayloadNamesAK4PFPuppiDATA4   = cms.FileInPath(options.optionJECAK4PFPuppiDATA4),
   jecPayloadNamesAK4PFPuppiDATAUnc = cms.FileInPath(options.optionJECAK4PFPuppiDATAUnc),
   jecPayloadNamesAK8PFchsDATA1   = cms.FileInPath(options.optionJECAK8PFchsDATA1),
   jecPayloadNamesAK8PFchsDATA2   = cms.FileInPath(options.optionJECAK8PFchsDATA2),
   jecPayloadNamesAK8PFchsDATA3   = cms.FileInPath(options.optionJECAK8PFchsDATA3),
   jecPayloadNamesAK8PFchsDATA4   = cms.FileInPath(options.optionJECAK8PFchsDATA4),
   jecPayloadNamesAK8PFchsDATAUnc = cms.FileInPath(options.optionJECAK8PFchsDATAUnc),
  # PILEUP REWEIGHTING
  PUReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/PileUpReweighting2017.root"),
  MinBiasUpReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/PileUpUpReweighting2017.root"),
  MinBiasDownReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/PileUpDownReweighting2017.root"),
  # PUPPI WEIGHT
  PuppiWeightFilePath = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUPPI/puppiCorr.root"),
  # BTAG REWEIGHTING
  BTAGReweightfile1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/csv_rwt_fit_hf_v2_final_2016_06_30test.root"),
  BTAGReweightfile2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/BTAGReweight/csv_rwt_fit_lf_v2_final_2016_06_30test.root"),
  # Object selection
  # Primary vertex cuts
  Pvtx_ndof_min   = cms.double(4.),
  Pvtx_vtx_max    = cms.double(24.),
  Pvtx_vtxdxy_max = cms.double(24.),
  # Obj primary vertex cuts
  vtx_ndof_min        = cms.int32(4),
  vtx_rho_max         = cms.int32(2),
  vtx_position_z_max  = cms.double(24.),
  # Muon cuts
  Muon_pt_min         = cms.double(5.),
  Muon_eta_max        = cms.double(50),
  # Electron cuts
  patElectron_pt_min  = cms.double(7.),
  patElectron_eta_max = cms.double(50),
  # Tau cuts
  Tau_pt_min          = cms.double(15.),
  Tau_eta_max         = cms.double(50.),
  # Jet cuts
  Jet_pt_min = cms.double(10.),
  # Photon cuts 
  Photon_pt_min   = cms.double(5.0),
  Photon_eta_max  = cms.double(5.0),
  # ttHFCategorization
  genJetPtMin               = cms.double(genJetPtMin),
  genJetAbsEtaMax           = cms.double(genJetAbsEtaMax),
  genJets                   = cms.InputTag(genJetCollection),
  genBHadJetIndex           = cms.InputTag("matchGenBHadron", "genBHadJetIndex"),
  genBHadFlavour            = cms.InputTag("matchGenBHadron", "genBHadFlavour"),
  genBHadFromTopWeakDecay   = cms.InputTag("matchGenBHadron", "genBHadFromTopWeakDecay"),
  genBHadPlusMothers        = cms.InputTag("matchGenBHadron", "genBHadPlusMothers"),
  genBHadPlusMothersIndices = cms.InputTag("matchGenBHadron", "genBHadPlusMothersIndices"),
  genBHadIndex              = cms.InputTag("matchGenBHadron", "genBHadIndex"),
  genBHadLeptonHadronIndex  = cms.InputTag("matchGenBHadron", "genBHadLeptonHadronIndex"),
  genBHadLeptonViaTau       = cms.InputTag("matchGenBHadron", "genBHadLeptonViaTau"),
  genCHadJetIndex           = cms.InputTag("matchGenCHadron", "genCHadJetIndex"),
  genCHadFlavour            = cms.InputTag("matchGenCHadron", "genCHadFlavour"),
  genCHadFromTopWeakDecay   = cms.InputTag("matchGenCHadron", "genCHadFromTopWeakDecay"),
  genCHadBHadronId          = cms.InputTag("matchGenCHadron", "genCHadBHadronId"),
)

#####
##   Dump gen particle list 
#####
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(-1),
  printVertex = cms.untracked.bool(True),
  src = cms.InputTag("prunedGenParticles")
)

#QG likelihood
process.load('BSMFramework.BSM3G_TNT_Maker.QGTagger_cfi')
process.QGTagger.srcJets       = cms.InputTag('slimmedJets')
process.QGTagger.jetsLabel     = cms.string('QGL_AK4PFchs')
#####
##   PROCESS
#####
process.p = cms.Path(
process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC *
#process.regressionApplication *
#process.calibratedPatElectrons  *
#process.egmGsfElectronIDSequence *
#process.electronIDValueMapProducer *
#process.egmPhotonIDSequence *
process.egammaScaleSmearAndVIDSeq *
process.fullPatMetSequence *
process.QGTagger *
process.rerunMvaIsolationSequence *
process.NewTauIDsEmbedded* # *getattr(process, "NewTauIDsEmbedded")
#process.selectedHadronsAndPartons*process.genJetFlavourInfos*process.matchGenCHadron*process.matchGenBHadron*
#process.primaryVertexFilter* 
#process.CSCTightHaloFilter*process.eeBadScFilter*process.HBHENoiseFilterResultProducer*process.ApplyBaselineHBHENoiseFilter*
process.TNT
)
