################# PISA ###################

for i in `gfal-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/vmariani/ZZ_TuneCP5_13TeV-pythia8/crab_ZZ/181007_082456/0000`;

do gfal-copy srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/vmariani/ZZ_TuneCP5_13TeV-pythia8crab_ZZ/181007_082456/0000/${i//\/cms/} /eos/user/v/vmariani/NTuples/ZZ/; done


