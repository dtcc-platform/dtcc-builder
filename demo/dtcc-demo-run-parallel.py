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
file = 'Majorna2021.json'
with open(file) as f:
  data = json.load(f)

print(data)
variable='MeshResolution'
if not variable in data:
    print(variable+" doesn't exist in "+file+", exiting!")
    quit()
#ids = []
#ids = [5.66, 11.31, 22.63, 45.25, 64.00, 90.51, 128.00]; 
ids= [4.00, 5.66, 8.00, 11.31, 16.00, 22.63, 32.00, 45.25, 64.00, 90.51, 128.00]
#ids = [ 181.02, 256.00, 362.04, 512.00, 724.08, 1024.00, 1448.15, 2048.00, 2896.31]
#ids=[28, 29];
print("Variable found, running...")
cwd=os.getcwd();
for id in ids:
    data[variable]=id;
    data['DataDirectory']=cwd+'/myCase'+str(id)
    with open(variable+str(id)+'.json', 'w') as outfile:
        json.dump(data, outfile)
    createCase('myCase'+str(id), '/home/dtcc/core/data/Majorna2021/')

    #process = subprocess.Popen(['echo', 'More output'],
    #                 stdout=subprocess.PIPE, 
    #                 stderr=subprocess.PIPE)
    #stdout, stderr = process.communicate()
    #stdout, stderr
    string=variable+str(id)+'.json'
    print("Spawning demo for "+variable+": "+str(id))
    with open(variable+str(id)+'.log', 'w') as f:
#        process = subprocess.Popen(['../bin/dtcc-generate-elevation-model',string], stdout=f, stderr=f)
        process = subprocess.Popen(['./demo-generate-mesh-arg',string], stdout=f, stderr=f)
        process.wait()

