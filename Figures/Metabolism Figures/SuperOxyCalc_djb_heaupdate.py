#! /usr/bin/env python

Usage='''
Something
Usage:
	Group2_djb.py HeaderLine.txt Infiles
'''
#####
#This update to the SuperOxyCalc_djb.py script takes each input .txt file from the Oxy
#and spits out a 'cleaned' version that includes the raw oxygen values corrected for the 
#calculations that are done by the original PreSens_cracked.xlsx spreadsheet
#This script no longer spits out the big compiled file


Debug=False

import sys, os, time
from math import *

StartTime=time.time()
def TimeConvert(TimeString):
	if "PM" in TimeString:
		TimeList=TimeString.split(" ")[0].split(":")
		if TimeList[0]=='12':
			MilTime='%s' %(":".join(TimeList))
		else:
			MilTime='%s:%s:%s' %(12+int(TimeList[0]),TimeList[1],TimeList[2])				
	else:
		if "AM" in TimeString:
			TimeList=TimeString.split(" ")[0].split(":")
			if TimeList[0]=='12':
				MilTime='0:%s' %(":".join(TimeList[1:]))
			else:
				if TimeList[0] in range(1,9):
					MilTime='0%s' %(":".join(TimeList))
				else:
					MilTime=TimeString.split(" ")[0]
	return MilTime

def OxyCalcPt1(Phase, Temp):
	Salin = float('35') #change this for different salinity
	Phase = float(Phase) #phase from datarow column 5 (python 4) E19 for first row
	Cal0 = float(60.62) #$B$6 from presense xcel sheet CHANGE FOR NEW SPOTS
	Cal100 = float(27.68) #$B$7 from presense xcel sheet CHANGE FOR NEW SPOTS
	AirPres = float(1013) #$B$8 from presense xcel sheet CHANGE FOR NEW SPOTS
	F1const = float(0.808)  #$B$11 from presense xcel sheet
	DelPhiK = float(-0.08255) #$B$12 from presense xcel sheet
	DelKsvK = float(0.000492) #$B$13 from presense xcel sheet
	Mconst = float(29.87) #$B$14 from presense xcel sheet
	Temp = float(Temp) #temp from datarow column 7 (python 6) G19 for first row
	T0Temp = float(20) #$E$6 from presense xcel sheet CHANGE FOR NEW SPOTS
	T100Temp = float(20) #$E$7 from presense xcel sheet CHANGE FOR NEW SPOTS
	KsvT100 = float(0.0514271) #$H$11 from presense xcel sheet


	PercAirSat = (-((tan(Phase*pi/180))/(tan((Cal0+(DelPhiK*(Temp-T0Temp)))*pi/180))*(KsvT100+(DelKsvK*(Temp-T100Temp)))+(tan(Phase*pi/180))/(tan((Cal0+(DelPhiK*(Temp-T0Temp)))*pi/180))*1/Mconst*(KsvT100+(DelKsvK*(Temp-T100Temp)))-F1const*1/Mconst*(KsvT100+(DelKsvK*(Temp-T100Temp)))-(KsvT100+(DelKsvK*(Temp-T100Temp)))+F1const*(KsvT100+(DelKsvK*(Temp-T100Temp))))+(sqrt((pow(((tan(Phase*pi/180))/(tan((Cal0+(DelPhiK*(Temp-T0Temp)))*pi/180))*(KsvT100+(DelKsvK*(Temp-T100Temp)))+(tan(Phase*pi/180))/(tan((Cal0+(DelPhiK*(Temp-T0Temp)))*pi/180))*1/Mconst*(KsvT100+(DelKsvK*(Temp-T100Temp)))-F1const*1/Mconst*(KsvT100+(DelKsvK*(Temp-T100Temp)))-(KsvT100+(DelKsvK*(Temp-T100Temp)))+F1const*(KsvT100+(DelKsvK*(Temp-T100Temp)))),2))-4*((tan(Phase*pi/180))/(tan((Cal0+(DelPhiK*(Temp-T0Temp)))*pi/180))*1/Mconst*pow((KsvT100+(DelKsvK*(Temp-T100Temp))),2))*((tan(Phase*pi/180))/(tan((Cal0+(DelPhiK*(Temp-T0Temp)))*pi/180))-1))))/(2*((tan(Phase*pi/180))/(tan((Cal0+(DelPhiK*(Temp-T0Temp)))*pi/180))*1/Mconst*pow((KsvT100+(DelKsvK*(Temp-T100Temp))),2)))
	PercO2 = PercAirSat*20.95/100
	P02HpA = (AirPres-exp(52.57-6690.9/(273.15+Temp)-4.681*log(273.15+Temp)))*PercAirSat/100*0.2095
	P02Torr = P02HpA/1.33322
	CO2MgL = ((AirPres-exp(52.57-6690.9/(273.15+Temp)-4.681*log(273.15+Temp)))/1013)*PercAirSat/100*0.2095*(48.998-1.335*Temp+0.02755*pow(Temp,2)-0.000322*pow(Temp,3)+0.000001598*pow(Temp,4))*32/22.414
	CO2umolL = CO2MgL*31.25
	Sheet1Data = ['%.6f'%PercAirSat, '%.7f'%PercO2, '%.6f'%P02HpA, '%.6f'%P02Torr, '%.8f'%CO2MgL, '%.6f'%CO2umolL]
	
	Chlorin = (Salin-0.03)/1.805
	FinalCO2MgL = ((AirPres-exp(52.57-6690.9/(273.15+Temp)-4.681*log(273.15+Temp)))/1013)*PercAirSat/100*0.2095*((49-1.335*Temp+0.02759*pow(Temp,2)-0.0003235*pow(Temp,3)+0.000001614*pow(Temp,4))-(Chlorin*(5.516*pow(10,-1)-1.759*pow(10,-2)*Temp+2.253*pow(10,-4)*pow(Temp,2)-2.654*pow(10,-7)*pow(Temp,3)+5.363*pow(10,-8)*pow(Temp,4))))*32/22.414
	FinalCO2umolL = FinalCO2MgL*31.25
	Sheet2FinalData = ['%.0f'%Salin,'%.7f'%Chlorin, '%.8f'%FinalCO2MgL, '%.6f'%FinalCO2umolL]
	return Sheet1Data, Sheet2FinalData

if len(sys.argv)<2:
	print Usage
else:
	InFileList=sys.argv[2:]
	InDelim=";"
	FileCount=0
	OutList=[]
	HeaderFile=open(sys.argv[1],'rU')
	HeaderLineCount=0
	#BigOut=open('CompiledOutfile.txt', 'w')
	for HeaderLine in HeaderFile:
		HeaderLineCount+=1
		if HeaderLineCount==1:
			Header=HeaderLine.rstrip()
		if HeaderLineCount==2:
			OxyCalcHeader=HeaderLine.rstrip()
	if Debug:
		print Header
	for InFile in InFileList:
		FileCount+=1
		LineCount=0
		FileToClean=open(InFile, 'rU')
		CleanedOut=open('%s_cleaned.txt' %(os.path.splitext(InFile)[0]),'w')
		StartParsing=0
		for Line in FileToClean:
			LineCount+=1
			if LineCount==2:
				ChamberID=Line.rstrip().split(" ")[1].split('-')[1]
				if Debug:
					print ChamberID
			if Line.startswith(Header[0:10]):
				if Debug:
					print Line.rstrip()
					print '%s_cleaned.txt' %(os.path.splitext(sys.argv[1])[0])
				CleanedOut.write('Chamber\t%s' %("\t".join(Line.rstrip().split(InDelim)[1:]))) # This prints the header line for the bigoutfile based on the first header line of the first input file
				CleanedOut.write('\t%s\n' %(OxyCalcHeader))
				#CleanedOut.write("%s\n" %("\t".join(Line.rstrip().split(InDelim)[1:])))
				StartParsing=LineCount
				#if FileCount==1:
					#CleanedOut.write('Chamber\t%s' %("\t".join(Line.rstrip().split(InDelim)[1:]))) # This prints the header line for the bigoutfile based on the first header line of the first input file
					#CleanedOut.write('\t%s\n' %(OxyCalcHeader))
			if StartParsing!=0 and LineCount>StartParsing:
				try:
					Cols=Line.rstrip().split(";")
					int(Cols[0].split(".")[0])
					#CleanedOut.write("%s\n" %("\t".join(Cols[1:])))
					CleanedOut.write('%s\t%s\t%s' %(ChamberID, TimeConvert(Cols[1]),"\t".join(Cols[2:])))
					Calcd, FinalCalcd=OxyCalcPt1(Cols[4], Cols[6])
					CleanedOut.write('\t%s\t%s\n' %("\t".join(Calcd), "\t".join(FinalCalcd)))
				except ValueError:
					StartParsing=LineCount
		FileToClean.close()
		CleanedOut.close()
	#BigOut.close()
	print "Compiled %d files into %s" %(FileCount, 'CompiledOutfile.txt')
	print "Time Elapsed: %f" %(time.time()-StartTime)
