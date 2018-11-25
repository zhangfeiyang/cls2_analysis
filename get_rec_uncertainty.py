#!/usr/bin/python

data_file0 = open('data0','r') # data of the ideal case of calibration
data_file2 = open('data2','r') # data of the more realistic case of calibration

data0s = data_file0.readlines()
data2s = data_file2.readlines()

for data0,data2 in zip(data0s,data2s):
	vares = [data0.split(),data2.split()]
	#[[energy0,mean0,emean0,sigma0,esigma0],[energy1,mean1,emean1,sigma1,esigma1]] = map(float,vares)
	[energy0,mean0,emean0,sigma0,esigma0] = map(float,vares[0])
	[energy1,mean1,emean1,sigma1,esigma1] = map(float,vares[1])
	#uncertainty1 = (sigma1/mean1/(sigma0/mean0) - 1)*100
	#uncertainty2 = esigma1/mean1/(sigma0/mean0)*100
	uncertainty1 = (sigma1/(sigma0)-1)*100
	uncertainty2 = esigma1/(sigma0)*100
	uncertainty = (uncertainty1**2+uncertainty2**2)**0.5
	print uncertainty

