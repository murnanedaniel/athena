#---------------
# install dumper
#---------------
download the installer in athena/Tracking/TrkDump/scripts/installer.sh
execute the installer from the locatation where you want to install athena
sh installer.sh

#---------------------
#to recompile the code
#---------------------
execute command line
athena/Tracking/TrkDumpAlgs/scripts/dumper-rebuild.sh

#-----------------
#to run the dumper
#-----------------
first parameter the files you want to generate in the file athena/Tracking/TrkDumpAlgs/scipts/run_reco.sh
set parameter to True or False to generate or not the files
'RAWtoESD:dumpObjects.csvFile=True/False'
'RAWtoESD:dumpObjects.rootFile=True/False'
athena/Tracking/TrkDumpAlgs/scripts/dumper-run.sh

#------------------------
#to compile the converter
#------------------------


#to execute the converter do

#to compare txt et dat files extrated from RDO do

