Per qualche ragione non riesco a pushare in gitub i JEC 2018, quindi li ho messi qui
/afs/cern.ch/user/v/vmariani/public/JEC/

Una volta che il framework funziona correttamente le Ntuple possono essere prodotte lanciato i job su crab da

/afs/cern.ch/user/v/vmariani/CMSSW_10_2_5/src/BSMFramework/BSM3G_TNT_Maker/crab/

con multicrab_MC18.py per i sample MC e multicrab_DATA18.py per i dati.

N.B.1 I sample single ele 2018 sono al momento mancanti (da cercare su DAS)
N.B.2 Tuttu i samples/json file/correzioni/ sono state updatate al 2018
N.B.3 sostituire i path di vmariani con i propri locali (in config.JobType.psetName)
      Nel cfg file miniAOD_MC2018.txt sostituire la riga 
       "fileName = cms.string("OutTree.root")" 
      con 
       fileName = cms.string(options.ofName+".root") in modo da attribuire il corretto nome al file di output.

!! Il cfg dei dati miniAOD_RD2018.txt va ancora aggiornato al 2018 !!

Una volta che i job per la produzione ntuple saranno finiti si possono produrre le rootuple con 
Rootplizer_HeavyNeutrino_SigTopDY.cc
Modificare dalla riga 1343 gli scale factor (muoni ed elettroni) e i pesi relativi ai vari campioni.
Lo script makelist.sh permette di creare un list*.txt che contiene l'elenco dei file su cui girare.
Sottomettere i job con lo script submission_macro.sh 

Unire i file di output e produrre i plot con la macro Analisi_Mu_HM.C che applicherà i tagli finali per i plot di analisi. 
