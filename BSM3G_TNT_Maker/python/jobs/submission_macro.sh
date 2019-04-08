
TASK="Sign_mumu_L13_M5000_" #name of the output directory and the files
INPUT_DATASET=(`cat listSign_mumu_L13_M5000.txt`) 

i=0
JOB_FOLDER=jobs_test #folder for jobs

OUTPUT_FOLDER=/eos/user/v/vmariani/NTuples/HN/$TASK
echo "$OUTPUT_FOLDER"

RESUBMIT=$1 

if [ ! -d "$OUTPUT_FOLDER" ] || [[ $1 = "resubmit" ]]; then
mkdir -p $JOB_FOLDER
mkdir -p $OUTPUT_FOLDER

echo "SUBMITTING: $TASK"

#this is for the resubmission
if [[ $1 = "resubmit" ]]; then
JOB_FOLDER=resubmit_$JOB_FOLDER
mkdir -p $JOB_FOLDER
rm -fr $JOB_FOLDER/*
fi

for file in ${INPUT_DATASET[*]}; do

i=$((i+1))

if [ ! -e $OUTPUT_FOLDER/${TASK}${i}.root ]; then

# Preparing job code
#name of your macro
cp Rootplizer_HeavyNeutrino_SigTopDY.cc $JOB_FOLDER/${TASK}${i}.cc

#this substitute inputfile and outpute file in your macro with the proper path
#echo "\"s/inputFile/${file}/\" $JOB_FOLDER/${TASK}${i}.cc" 
sed -i -e "s/inputFile/${file}/" $JOB_FOLDER/${TASK}${i}.cc
sed -i -e "s/outputFile/${TASK}${i}.root/" $JOB_FOLDER/${TASK}${i}.cc
sed -i -e "s/filename_/${TASK}${i}/" $JOB_FOLDER/${TASK}${i}.cc

echo "#!/bin/bash 
cd $PWD 
export SCRAM_ARCH=slc6_amd64_gcc481
eval \`scram runtime -sh\`

export X509_USER_PROXY=/afs/cern.ch/user/v/vmariani/proxy

cd -

root -b ${TASK}${i}.cc

" > $JOB_FOLDER/${TASK}${i}.sh 

chmod +x $JOB_FOLDER/${TASK}${i}.sh

echo "#!/bin/bash
mv ${TASK}${i}.root $OUTPUT_FOLDER/${TASK}${i}.root
" > $JOB_FOLDER/copy_${TASK}${i}.sh
 
chmod +x $JOB_FOLDER/copy_${TASK}${i}.sh

if [[ $((i%30)) -eq 0 ]]; then
echo "prepared $JOB_FOLDER/${TASK}${i}.sh"
fi

if [[ $1 = "resubmit" ]]; then
echo "prepared $JOB_FOLDER/${TASK}${i}.sh"
fi

## Uncomment to limit the number of jobs
#if [[ $i = 30 ]]; then
#break
#fi

else

echo "${i} CHECKED, skipping"

fi

done


## submit job
echo "Submitting task $TASK"
echo "executable      = $JOB_FOLDER/\$Fn(filename).sh
arguments             = \$(ClusterID) \$(ProcId)
output                = /dev/null
error                 = condor${TASK}.err 
log                   = condor${TASK}.log
transfer_input_files    = $JOB_FOLDER/copy_\$Fn(filename).sh,$JOB_FOLDER/\$Fn(filename).cc
+PostCmd  = \"copy_\$Fn(filename).sh\"
+JobFlavour             = \"longlunch\"
queue filename matching ($JOB_FOLDER/${TASK}*.sh)
" > $JOB_FOLDER/condor${TASK}.sub

condor_submit $JOB_FOLDER/condor${TASK}.sub

else
echo "DIR $OUTPUT_FOLDER ALREADY EXISTS, SKIPPING"

fi
