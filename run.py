#!/usr/bin/python
import os
sources = ['Ge68','Co60','Cs137','Mn54','K40']
#sources = ['K40']
#file = open('log0','r') 
#lines = file.readlines()
lines = ['0_0']
for line in lines:
	for source in sources:
		R = line.split('_')[0]
		Z = line.split('_')[1].strip('\n')
		cmd = './fit_cb '+source+' '+R+' '+Z
		os.system(cmd)
		cmd = './fit_emc '+source+' '+R+' '+Z
		os.system(cmd)
		cmd = './fit_nocom '+source+' '+R+' '+Z
		os.system(cmd)
