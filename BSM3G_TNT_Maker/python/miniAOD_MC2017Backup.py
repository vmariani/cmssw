import FWCore.ParameterSet.Config as cms
#####
##   Initial standard configs
#####
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
process.prefer("GlobalTag")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

#####
##   Input files
#####
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    '/store/mc/RunIIFall17MiniAOD/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/0CF65340-0200-E811-ABB7-0025905C53F0.root'
  ),
  skipEvents = cms.untracked.uint32(0)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


##### JEC
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
  process,
  jetSource = cms.InputTag('slimmedJets'),
  labelName = 'UpdatedJEC',
  jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet','L2Relative','L3Absolute']), 'None')
)
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

from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
process.load('RecoTauTag.Configuration.loadRecoTauTagMVAsFromPrepDB_cfi')
from RecoTauTag.RecoTau.PATTauDiscriminationByMVAIsolationRun2_cff import *
def loadMVA_WPs_run2_2017(process):
               
      for training, gbrForestName in tauIdDiscrMVA_trainings_run2_2017.items():

         process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
            cms.PSet(
               record = cms.string('GBRWrapperRcd'),
               tag = cms.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_version)),
               label = cms.untracked.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_version))
            )
         )

         for WP in tauIdDiscrMVA_WPs_run2_2017[training].keys():
            process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
               cms.PSet(
                  record = cms.string('PhysicsTGraphPayloadRcd'),
                  tag = cms.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_version, WP)),
                  label = cms.untracked.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_version, WP))
               )
            )

         process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
            cms.PSet(
               record = cms.string('PhysicsTFormulaPayloadRcd'),
               tag = cms.string("RecoTauTag_%s%s_mvaOutput_normalization" % (gbrForestName, tauIdDiscrMVA_2017_version)),
               label = cms.untracked.string("RecoTauTag_%s%s_mvaOutput_normalization" % (gbrForestName, tauIdDiscrMVA_2017_version))
            )
        )


# 2017v2 dR=0.3
tauIdDiscrMVA_2017_version = "v2"
tauIdDiscrMVA_trainings_run2_2017 = {
       'tauIdMVAIsoDBoldDMdR0p3wLT2017' : "tauIdMVAIsoDBoldDMdR0p3wLT2017",
}

tauIdDiscrMVA_WPs_run2_2017 = {
       'tauIdMVAIsoDBoldDMdR0p3wLT2017' : {
      'Eff95' : "DBoldDMdR0p3wLTEff95",
      'Eff90' : "DBoldDMdR0p3wLTEff90",
      'Eff80' : "DBoldDMdR0p3wLTEff80",
      'Eff70' : "DBoldDMdR0p3wLTEff70",
      'Eff60' : "DBoldDMdR0p3wLTEff60",
      'Eff50' : "DBoldDMdR0p3wLTEff50",
      'Eff40' : "DBoldDMdR0p3wLTEff40"
    }
}

loadMVA_WPs_run2_2017(process)

process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2raw = patDiscriminationByIsolationMVArun2v1raw.clone(
   PATTauProducer = cms.InputTag('slimmedTaus'),
   Prediscriminants = noPrediscriminants,
   loadMVAfromDB = cms.bool(True),
   mvaName = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2"),
   mvaOpt = cms.string("DBoldDMwLTwGJ"),
   requireDecayMode = cms.bool(True),
   verbosity = cms.int32(0)
)

process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VLoose = patDiscriminationByIsolationMVArun2v1VLoose.clone(
   PATTauProducer = cms.InputTag('slimmedTaus'),
   Prediscriminants = noPrediscriminants,
   toMultiplex = cms.InputTag('rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2raw'),
   key = cms.InputTag('rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2raw:category'),
   loadMVAfromDB = cms.bool(True),
   mvaOutput_normalization = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_mvaOutput_normalization"),
   mapping = cms.VPSet(
      cms.PSet(
      category = cms.uint32(0),
      cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff90"),
      variable = cms.string("pt"),
      )
   ),
   verbosity = cms.int32(0)
)

process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VVLoose = process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VLoose.clone()
process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VVLoose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff95")
process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Loose = process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VLoose.clone()
process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Loose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff80")
process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Medium = process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VLoose.clone()
process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Medium.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff70")
process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Tight = process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VLoose.clone()
process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Tight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff60")
process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VTight = process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VLoose.clone()
process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff50")
process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VVTight = process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VLoose.clone()
process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VVTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff40")


process.rerunMvaIsolationSequence = cms.Sequence(
   process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2raw
   *process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VLoose
   *process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VVLoose
   *process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Loose
   *process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Medium
   *process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Tight
   *process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VTight
   *process.rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VVTight
)

# embed new id's into new tau collection
embedID = cms.EDProducer("PATTauIDEmbedder",
   src = cms.InputTag('slimmedTaus'),
      tauIDSources = cms.PSet(
            byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017 = cms.InputTag('rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2raw'),
            byVVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = cms.InputTag('rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VVLoose'),
            byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = cms.InputTag('rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VLoose'),
            byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = cms.InputTag('rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Loose'),
            byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = cms.InputTag('rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Medium'),
            byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = cms.InputTag('rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2Tight'),
            byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = cms.InputTag('rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VTight'),
            byVVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017 = cms.InputTag('rerunDiscriminationByIsolationOldDMdR0p3MVArun2017v2VVTight')
   ),
 )

setattr(process, "NewTauIDsEmbedded", embedID)

# inclusion in the process
#process.p += process.rerunMvaIsolationSequence
#process.p += getattr(process, "NewTauIDsEmbedded")

####
## Output
#####

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("OutTree.root")
)

#####
##   Analysis parameters
#####
process.TNT = cms.EDAnalyzer("BSM3G_TNT_Maker",
  #### Running options
  # Choose which trigger you want (do NOT need to put * as it will consider all the versions by default)
  ifevtriggers      = cms.bool(False), # True means you want to require the triggers
  maxtriggerversion = cms.double(20), # please leave it as a double
  evtriggers        = cms.vstring(
     #'HLT_Ele115_CaloIdVT_GsfTrkIdT_v',
     'HLT_DoubleEle33_CaloIdL_MW_v',
     #'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v',
     #'HLT_IsoMu24_v',
     #'HLT_IsoTkMu24_v',
     'HLT_Mu50_v',
     'HLT_TkMu50_v',
     'HLT_Mu30_TkMu11_v',
  ),
  # Choose which information you want to use
  fillgeninfo           = cms.bool(True),
  fillgenHFCategoryinfo = cms.bool(False),
  filleventinfo         = cms.bool(True),
  filltriggerinfo       = cms.bool(True),
  fillPVinfo            = cms.bool(True),
  fillmuoninfo          = cms.bool(True),
  fillelectronpatinfo   = cms.bool(True),
  filltauinfo           = cms.bool(True),#FF
  filljetinfo           = cms.bool(True), #T
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
  is_data   = cms.bool(False),
  reHLT     = cms.bool(True),
  debug_    = cms.bool(False),
  super_TNT = cms.bool(False),
  AJVar     = cms.bool(False),
  tthlepVar = cms.bool(True),#FF
  bjetnessselfilter = cms.bool(False),
  PuppiVar  = cms.bool(False),
  qglVar    = cms.bool(True),#FF
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
  metFilterBits       = cms.InputTag("TriggerResults", "", "PAT"),
  photons             = cms.InputTag("slimmedPhotons"),
  packedPFCandidates  = cms.InputTag("packedPFCandidates"), 
  pruned              = cms.InputTag("prunedGenParticles"),
  # Tau
  #byIsolationMVArun2017v2DBoldDMwLTraw2017 = cms.string('byIsolationMVArun2017v2DBoldDMwLTraw2017'),
  #byVVLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byVVLooseIsolationMVArun2017v2DBoldDMwLT2017'),
  #byVLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byVLooseIsolationMVArun2017v2DBoldDMwLT2017'),
  #byLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byLooseIsolationMVArun2017v2DBoldDMwLT2017'),
  #byMediumIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byMediumIsolationMVArun2017v2DBoldDMwLT2017'),
  #byTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byTightIsolationMVArun2017v2DBoldDMwLT2017'),
  #byVTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byVTightIsolationMVArun2017v2DBoldDMwLT2017'),
  #byVVTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byVVTightIsolationMVArun2017v2DBoldDMwLT2017'),
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
  jecPayloadNamesAK4PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L1FastJet_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2Relative_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L3Absolute_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2L3Residual_AK4PFchs.txt"),
  jecPayloadNamesAK4PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_Uncertainty_AK4PFchs.txt"),
  jecPayloadNamesAK4PFPuppiMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L1FastJet_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2Relative_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L3Absolute_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2L3Residual_AK4PFPuppi.txt"),
  jecPayloadNamesAK4PFPuppiDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_Uncertainty_AK4PFPuppi.txt"),
  jecPayloadNamesAK8PFchsMC1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMC3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsMCUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_Uncertainty_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA1   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L1FastJet_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA2   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2Relative_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA3   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L3Absolute_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATA4   = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2L3Residual_AK8PFchs.txt"),
  jecPayloadNamesAK8PFchsDATAUnc = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/JEC/DATA/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_Uncertainty_AK8PFchs.txt"),
  # PILEUP REWEIGHTING
  PUReweightfile      = cms.FileInPath("BSMFramework/BSM3G_TNT_Maker/data/PUReweight/PileUpReweighting2017.root"),
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
  patElectron_pt_min  = cms.double(5.),
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
process.egmGsfElectronIDSequence *
process.electronIDValueMapProducer *
#process.egammaScaleSmearAndVIDSeq *
process.QGTagger *
process.rerunMvaIsolationSequence *
getattr(process, "NewTauIDsEmbedded")*
#process.selectedHadronsAndPartons*process.genJetFlavourInfos*process.matchGenCHadron*process.matchGenBHadron*
#process.primaryVertexFilter* 
#process.CSCTightHaloFilter*process.eeBadScFilter*process.HBHENoiseFilterResultProducer*process.ApplyBaselineHBHENoiseFilter*
process.TNT
)
