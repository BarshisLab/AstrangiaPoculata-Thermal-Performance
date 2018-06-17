#! /usr/bin/env python

###Takes a list of hobo .csv files and cleans them, outputs as _clean.txt files, and writes the first few lines of the r script for merging in a new file called yourfilename.r
###Usage hobocleaner.py yourRfilename.R anynumberoffiles

import os, sys

filelist=sys.argv[2:]
rscript=open(sys.argv[1],'w')
filecount=0
for file in filelist:
	filecount+=1
	infile=open(file, 'r')
	outfile=open('%s_clean.txt' %(file[:-4]), 'w')
	linecount=0
	for line in infile:
		linecount+=1
		if linecount==2:
			outfile.write('%s\t%s\t%s\n'%('RecNum','DateTime',file[:-4]))
		if linecount>2:
			cols=line.rstrip().split('	')
			if cols[2]=='':
				continue
			else:
				outfile.write('%s\n'%('\t'.join(cols[0:3])))
	rscript.write('boing%d<-read.delim("%s_clean.txt")\nboing%d$DateTime<-strptime(boing%d$DateTime, format="%s")\n' %(filecount, file[:-4], filecount, filecount, '%m/%d/%y %I:%M:%S %p'))
	infile.close()
	outfile.close()
	if filecount==2:
		rscript.write('spoing<-merge(boing1[2:3], boing%d[2:3], by="DateTime", all=T)\n'%(filecount))
	if filecount>2:
		rscript.write('spoing<-merge(spoing, boing%d[2:3], by="DateTime", all=T)\n'%(filecount))