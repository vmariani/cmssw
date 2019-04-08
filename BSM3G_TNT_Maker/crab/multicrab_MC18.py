if __name__ == '__main__':
 #####
 ##   Multicrab configuration
 #####
 import sys
 from CRABClient.UserUtilities import config, getUsernameFromSiteDB
 config = config()
 from CRABAPI.RawCommand import crabCommand
 from CRABClient.ClientExceptions import ClientException
 from httplib import HTTPException
 config.General.workArea = 'Crab_projects'

 def submit(config):
  try:
   crabCommand('submit', config = config)
  except HTTPException as hte:
   print "Failed submitting task: %s" % (hte.headers)
  except ClientException as cle:
   print "Failed submitting task: %s" % (cle)
 #####
 ##   Crab configuration
 #####
 datasetnames  = [
# signal
#'Fall17V2_TTHnobb',
# prompt background
# TT+ST
'Autumn18_TT',
'Autumn18_ST',
'Autumn18_SaT',
# WW+WZ+ZZ
'Autumn18_WW',
'Autumn18_WZ',
'Autumn18_ZZ'
# Wjets
'Autumn18_WJetsHT70To100',
'Autumn18_WJetsHT100To200',
'Autumn18_WJetsHT200To400',
'Autumn18_WJetsHT400To600',
'Autumn18_WJetsHT600To800',
'Autumn18_WJetsHT800To1200',
'Autumn18_WJetsHT1200To2500',
'Autumn18_WJetsHT2500ToInf',
'Autumn18_WJets_',
#DY
'Autumn18_DY',                 
]
 datasetinputs = [
# signal
#'/ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM',
# prompt background
# TT
'/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
'/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM',
'/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM',
#WW
'/WW_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM',
'/WZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM',
'/ZZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM',
#WJETS
'/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
'/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
'/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
'/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
'/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
'/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
'/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
'/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
'/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM',
#DY
'/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM',
                ]

for d in range(36,len(datasetnames)):
#for d in range(0,38):
    print 'multicrab.py: Running datasetname: ', datasetnames[d]

    config.section_('General')
    config.General.requestName = datasetnames[d]
    config.General.workArea    = datasetnames[d]
    config.General.transferLogs = True

    config.section_('JobType')
    config.JobType.pluginName  = 'Analysis'
    # List of parameters to pass to CMSSW parameter-set configuration file:
    config.JobType.psetName    = '/afs/cern.ch/user/v/vmariani/CMSSW_10_2_5/src/BSMFramework/BSM3G_TNT_Maker/python/miniAOD_MC2018.py'
    #config.JobType.inputFiles = ['/afs/cern.ch/work/b/binghuan/private/BSMFWTest/CMSSW_9_4_10/src/BSMFramework/BSM3G_TNT_Maker/data/L1Prefire/L1PrefiringMaps_new.root']
    config.JobType.sendExternalFolder = True
    config.JobType.maxMemoryMB = 2000 # Default == 2Gb : maximum guaranteed to run on all sites
    #config.JobType.allowUndistributedCMSSW = True
    ofParam = 'ofName=' + datasetnames[d]
    config.JobType.pyCfgParams = [ofParam]
    config.section_('Data')
    config.Data.allowNonValidInputDataset = True
    config.Data.inputDataset   = datasetinputs[d]
    config.Data.inputDBS       = 'global'
    config.Data.splitting      = 'FileBased'
    config.Data.totalUnits     = 40000 #With 'FileBased' splitting tells how many files to analyse
    config.Data.unitsPerJob    = 1
    config.Data.outLFNDirBase = '/store/user/vmariani/'# First part of LFN for output files (must be /store/user/<username>/ or /store/group/<username>/  )
    config.Data.outputDatasetTag = datasetnames[d]

    print 'multicrab.py: outLFNDirBase = /store/user/vmariani/'
    #config.Data.publication = True

    config.section_('Site')
    config.Site.storageSite    = 'T2_IT_Bari'#'T2_CH_CERN' # Site to which output is permenantly copied by crab3
    print 'multicrab.py: Submitting Jobs'
    submit(config)
