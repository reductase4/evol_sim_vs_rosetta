import os, re
import subprocess
from subprocess import Popen
from numpy import *
#from pylab import *
from scipy.stats import pearsonr as pearson
#from scipy.stats import polyfit

#list of files that contain L vs RSA data
data = []

#search string to use
searchStr = "^align_data_array_ordered" + "[a-zA-Z0-9_\.\-]*" + "_evolved.dat"

x = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

fpO = open("raw_mle_lines_ordered_evolved.csv","w")

#find all csv files that match the search string
for path, names, filename in os.walk('.',False):
    
    for file in filename:
        #print "Testing file: %s" %file
        if(re.search(searchStr, file)!=None):

        	print "Found file: %s" % file
        	#print Popen(["which Rscript"],shell=True, stdout=subprocess.PIPE).communicate()
        	#output = Popen(["/share/apps/R-2.11.1/bin/Rscript MLE_calc.R " + file],shell=True, stdout=subprocess.PIPE).communicate()
	    	#output = Popen(["/home/elj299/RSA_from_ceres/data/rsa_analysis/ MLE_calc.R " + file], shell=True, stdout=subprocess.PIPE).communicate()
        	#print output
	    	output = subprocess.Popen(["Rscript /Users/qian/Desktop/evol_sim_vs_rosetta/r_scripts/MLE_calc.R " + file], shell=True, stdout=subprocess.PIPE).communicate()
	    	print output

        	result = re.findall("\-?[0-9]+\.?[0-9]*", re.split("\n",output[0])[-2])
        	print result
        	slop = result[0]
        	int = result[1]
            
        	print file + " y(x) = " + str(int) + " + x" + str(slop)
        	fpO.write(file+","+str(int)+","+str(slop)+"\n")
        	data.append([int, slop])

fpO.close()

#print transpose(matrix(data))

#plot(transpose(matrix(data)),'ob')
#xlim(-2,2)
#ylim(-2,2)
#savefig("mle_output.svg")
#raw_input()
