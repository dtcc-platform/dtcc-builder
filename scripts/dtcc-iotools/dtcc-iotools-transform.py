import subprocess
import sys
if len(sys.argv)<5 or len(sys.argv)>5:
	print("Usage vc-transform filein fileout coordin coordout")	
	exit()
filein=sys.argv[1];
fileout=sys.argv[2];
coordin=sys.argv[3];
coordout=sys.argv[4];
string2='gdaltransform -s_srs EPSG:'+coordin+' -t_srs EPSG:'+coordout+' < '+filein+' > '+fileout+' -output_xy';
print(string2);
p = subprocess.Popen(string2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
text = p.communicate()[0]
print(text)
retval=p.wait()
