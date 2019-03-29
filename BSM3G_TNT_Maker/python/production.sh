first=1
last=20
inType="TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/TT/171020_102139/"
outType="TT"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/v/vmariani/CMSSW_9_4_9_cand2/src/BSMFramework/BSM3G_TNT_Maker/python 
#/home/roberto.leonardi/CMSSW_7_4_12/src/rootple
eval `scramv1 runtime -sh`

#chmod 600 /home/roberto.leonardi/CMSSW_7_4_12/src/rootple/tmp_proxy
export X509_USER_PROXY=/afs/cern.ch/user/v/vmariani/proxy
voms-proxy-info

pathin=root://xrootd-redic.pi.infn.it//store/user/vmariani/$inType
pathout=/eos/user/v/vmariani/HN/Rootuple/
rootplizer=/afs/cern.ch/user/v/vmariani/CMSSW_9_4_9_cand2/src/BSMFramework/BSM3G_TNT_Maker/python/Rootplizer_HeavyNeutrino_SigTopDY.cc
#rootplizer=/home/roberto.leonardi/CMSSW_7_4_12/src/rootple/Rootplizer_HeavyNeutrino_QCD.cc
#FileName=eleB1

for ((f = $first ; f <= $last ; f++));
do
  if [ "$f" -lt "1000" ];
  then
      block="0000"
  else
      block="0001"
  fi
  string="(\"$pathin"$block"/OutTree_"$f".root\",\"$pathout""/"$outType"_"$f".root\")"
  echo $string
  #root -b -q -l $rootplizer$string
done
