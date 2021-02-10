#!/usr/bin/python
import json
import subprocess
#import psutil
import shutil
import os

def createCase(caseName,dataDir):
    print('Creating '+caseName)
    shutil.copytree(dataDir,caseName)
    return

def checkIfProcessRunning(processName):
    '''
    Check if there is any running process that contains the given name processName.
    '''
    #Iterate over the all the running process
    for proc in psutil.process_iter():
        try:
            # Check if process name contains the given name string.
            if processName.lower() in proc.name().lower():
                return True
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    return False;

data = {}
file = 'Hammarkullen2020.json'
with open(file) as f:
  data = json.load(f)

print(data)
variable='ElevationModelResolution'
if not variable in data:
    print(variable+" doesn't exist in "+file+", exiting!")
    quit()
ids = [0.5, 1, 2, 4, 8, 16]
print("Variable found, running...")
cwd=os.getcwd();
for id in ids:
    data[variable]=id;
    data['DataDirectory']=cwd+'/myCase'+str(id)
    with open(variable+str(id)+'.json', 'w') as outfile:
        json.dump(data, outfile)
    createCase('myCase'+str(id), '/home/dtcc/core/data/Hammarkullen2020/')

    #process = subprocess.Popen(['echo', 'More output'],
    #                 stdout=subprocess.PIPE, 
    #                 stderr=subprocess.PIPE)
    #stdout, stderr = process.communicate()
    #stdout, stderr
    string=variable+str(id)+'.json'
    print("Spawning demo for "+variable+": "+str(id))
    with open(variable+str(id)+'.log', 'w') as f:
#        process = subprocess.Popen(['../bin/dtcc-generate-elevation-model',string], stdout=f, stderr=f)
        process = subprocess.Popen(['./dtcc-demo-pipeline-1arg',string], stdout=f, stderr=f)

