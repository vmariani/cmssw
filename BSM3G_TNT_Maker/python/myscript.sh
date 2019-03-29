# Actually, let's do something more useful than a simple hello world... this will print the input arguments passed to the script
echo "Here there are all the input arguments"
echo $@

# If you are curious, you can have a look at the tweaked PSet. This however won't give you any information...
echo "================= PSet.py file =================="
cat PSet.py

# This is what you need if you want to look at the tweaked parameter set!!
echo "================= Dumping PSet ===================="
python -c "import PSet; PSet.process.source.fileNames=[x.replace('/store', 'root://ds-202-06-19.cr.cnaf.infn.it:31094//store') for x in PSet.process.source.fileNames]; import pickle; pickle.dump(PSet.process, open('PSet.pkl','w'))"
#python -c "import PSet; PSet.process.source.fileNames=[x.replace('/store', 'root://cloud-vm90.cloud.cnaf.infn.it:31194//store') for x in PSet.process.source.fileNames]; import pickle; pickle.dump(PSet.process, open('PSet.pkl','w'))"

# Ok, let's stop fooling around and execute the job:
cmsRun -j FrameworkJobReport.xml -p PSet.py
