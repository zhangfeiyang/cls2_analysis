#!/usr/bin/env python
from ROOT import *

file0 = open('../Cs137/0_0/result_nocom','r')
file1 = open('../Cs137/0_0/result_emc','r')

line0s = file0.readlines()
line1s = file1.readlines()

h = TH1F('h','',50,0,0)
for i,line0 in enumerate(line0s):
	data0 = line0.split()[4]
	data1 = line1s[i].split()[6]
	print(data0,data1)
	#h.Fill(float(data0)-float(data1))
	h.Fill(float(data0))

h.Draw()

a = input()


