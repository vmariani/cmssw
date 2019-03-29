rm listWJetsToLNu.txt;

################# PISA ###################

for i in `gfal-ls srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/vmariani/mu_L20_M8000/private_mu_L20M8000/170731_150839/0000`;

do echo "root://xrootd.ba.infn.it//store/user/roleonar/mu_L20_M8000/private_mu_L20M8000/170731_150839/0000/${i//\/cms/}" >> listSign_mumu_L20_M8000.txt; done

sed -i -e "s/\"/\\\\\"/g" listSign_mumu_L20_M8000.txt
sed -i -e "s/\//\\\\\//g" listSign_mumu_L20_M8000.txt

