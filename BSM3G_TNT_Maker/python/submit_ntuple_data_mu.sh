arr=(/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD /SingleMuon/Run2017C-31Mar2018-v1/MINIAOD /SingleMuon/Run2017D-31Mar2018-v1/MINIAOD /SingleMuon/Run2017E-31Mar2018-v1/MINIAOD /SingleMuon/Run2017F-31Mar2018-v1/MINIAOD)
sam=(SingleMu_Run2017B SingleMu_Run2017C SingleMu_Run2017D SingleMu_Run2017E SingleMu_Run2017F)

work="crab_projects_skim"

mkdir -p $work  

i=0
rm log.out
touch log.out

for item in ${arr[*]}

do
        DATASET=$item
#arr[*]


        echo "
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = '${sam[$i]}'
config.General.workArea = '$work'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniAOD_RD2017.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.priority = 10;

config.Data.inputDataset = '$DATASET'
config.Data.inputDBS = 'global'
config.Data.splitting='Automatic'
config.Data.publication = False
config.Data.lumiMask = \"/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_MuonPhys.txt\"

config.Site.storageSite = \"T2_IT_Pisa\"
" > crab_cfg.py


crab submit -c crab_cfg.py 

i=$((i+1))
done

