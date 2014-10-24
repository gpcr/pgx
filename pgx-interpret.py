#!/usr/bin/python

import sys
import argparse
import vcf
import vcf.utils
import pgxClasses
import re

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-v","--verbose")
  parser.add_argument("--ref_vcf", type=argparse.FileType('r'), default='pgx_genotyper_input.vcf')
  parser.add_argument("--ref_site_filename", type=str, default='pgxRefSites.txt')
  parser.add_argument("--gene_filename", type=str, default='pgxGenes.txt')
  parser.add_argument("--out_filename", type=str, default='out.txt')
  parser.add_argument("-d","--depth", type=int, default=30)
  parser.add_argument('vcf_file', nargs='+', type=argparse.FileType('r'))
  args = parser.parse_args()

  refSites = readInRefSites(args.ref_site_filename)
  delimRE = r'[\s,;\/\\]+'
  genes = readInGenes(args.gene_filename, delimRE)
  outFile = open(args.out_filename,'w')
  outFile.write('Gene\tVariant ID (rs#)\tVariant (genomic position)\tVariant (transcript' 
  + ' position)\tZygosity\tRead Depth by reference base, variant\tAssociated *Alleles\n')
  
  # Set up reference tables
  tagsites = list(vcf.Reader(args.ref_vcf))
  for samplefile in args.vcf_file:
    # Process a sample
    sample = vcf.Reader(samplefile)
    if len(sample.samples) > 1:
      print 'Interpretation of multi-sample VCF files not implemented'
      sys.exit(1)

    ID = sample.samples[0]
#    print ID
    for (refrecord, record) in vcf.utils.walk_together(iter(tagsites),sample, kwargs=lambda r: (r.CHROM, r.POS, r.REF, R.ALT)):
      if refrecord is None:
         print 'ERROR: novel variant?', record
         sys.exit(1)
      if record is None:
#         if args.verbose or 'AscAlleles' in refrecord.INFO: 
#          print 'WARNING: tagsite missing', refrecord
         continue 

      depth = 0
      if 'DP' in record.INFO:
        depth = record.INFO['DP']
      if depth < args.depth:
         if args.verbose or 'AscAlleles' in refrecord.INFO:
           print 'WARNING: Low reads on tag site!', record, record.genotype(ID)

      if 'AscAlleles' in refrecord.INFO:
        call = record.genotype(ID)
        if call.gt_type is 0:
          spacer = '\t'
        else:
          spacer = ''
        #print spacer, ID, refrecord.INFO['Gene'], refrecord.INFO['AscAlleles'], call
    	#print refrecord.ALT, record.ALT
    
    
    	#find reference with the same set of *alleles as the record.
    	if (call.gt_type is not 0):
    		recordReadDepthList = call['AD'] 		
    		#find associated alleles for just the variant allele with the most reads
    		assocAllelesRecordSet = set(processStrToList(
    		refrecord.INFO['AscAlleles'][recordReadDepthList[1:].index(max(recordReadDepthList[1:]))], delimRE))
    		for refSite in refSites:
				assocAllelesRefSet = set(processStrToList(refSite.associatedAlleles, delimRE))
				if (assocAllelesRecordSet == assocAllelesRefSet
				and refSite.pgxGene == refrecord.INFO['Gene'] and call['AD']):							
					#output results to file	
					writeResults(outFile, refSite, call['AD'])	
						
  outFile.close()



#reads in tab-separated reference site file. in each row: associatedAlleles, cypVariantID, 
#gene, rsVariantID, variantGenomicPos, variantTranscriptPos.
#returns list of pgxClasses.pgxRefSite objects.
def readInRefSites(filename):
	fields = pgxClasses.pgxRefSite().__dict__.keys()
	refSites=[]
	f = open(filename)
	for rowNum, row in enumerate(f, start=1):
		if row==['']: #ignore blank rows
			return
		row = row.strip().split('\t')
		checkRowLen(filename, row, rowNum, fields, '\\t')
		refSites += [pgxClasses.pgxRefSite(*row)]
	return refSites
	f.close()


#reads in tab-separated genes file.
def readInGenes(filename, delimRE):
	fields = pgxClasses.pgxGene().__dict__.keys()
	genes=[]
	geneStartSymbol = '>'
	f = open(filename)
	for rowNum, row in enumerate(f, start=1):
		row = row.strip().split('\t')
		if row==['']: #ignore blank rows
			continue
		if row[0][0] == geneStartSymbol:
			row[0] = row[0][1:]
			genes += [pgxClasses.pgxGene()]
 		if row[0] not in fields:
			throwError(filename, rowNum, 'First value in row must be a pgxGene field.\nAcceptable values: '	+ ', '.join(fields))	
		#dict of pgxGene field : expected row format
		rowFormatDict = dict([('name', ['name', '[geneName]']), ('defaultAllele', ['defaultAllele', '[default*Allele]']), 
		('alleleFunc', ['alleleFunc', '[funcLevelID]', '[label]', '[assoc*Alleles]']), ('phen', ['phen', '[phenotypeLabel]', '[alleleFuncPairs]'])])	
		checkRowLen(filename, row, rowNum, rowFormatDict[row[0]], '\\t')
		
		#populate gene object
		if not genes:
			throwError(filename, rowNum, 'New gene specifications must start with symbol ' + geneStartSymbol)
		if ((row[0] == 'name' and genes[-1].name != 'undefined') or (row[0] == 'defaultAllele' and genes[-1].defaultAllele != 'undefined')):
			throwError(filename, rowNum, 'Contradictory / overlapping gene specifications.\n' + 
			'New gene specifications must start with symbol ' + geneStartSymbol)
		geneObj = genes[-1]
		if row[0] == 'name':
			if (row[1] in [z.name for z in genes[:-1]]):
				throwError(filename, rowNum, 'Gene name already exists. Unique identifier required.')			
			geneObj.name = row[1]
		elif row[0] == 'defaultAllele':
			geneObj.defaultAllele = row[1]
		elif row[0] == 'alleleFunc':
			if geneObj.alleleFunc == 'undefined':
				geneObj.alleleFunc = dict([(row[1], [row[2], processStrToList(row[3], delimRE)])])
			elif row[1] not in geneObj.alleleFunc:
				geneObj.alleleFunc[row[1]] = [row[2], processStrToList(row[3], delimRE)]		
			else:
				throwError(filename, rowNum, 'Specified allele function already exists for gene.')
		elif row[0] == 'phen':
			#process and validate allele function pairs
			funcPairs = processStrToList(row[2], delimRE)
			for x in funcPairs:
				if x not in geneObj.alleleFunc:
					throwError(filename, rowNum, 'Unrecognized allele function ID: ' + x + 
					'\nAllele functions must be specified before phenotypes.')
			if len(funcPairs)%2 != 0:
				throwError(filename, rowNum, 'Couldn\'t parse allele function pairs.\n' + 
				'Expected format: [funcID x] [funcID y], [funcID z] [funcID w], ...')		
			funcPairs = [[funcPairs[i], funcPairs[i+1]] for i in range(0, len(funcPairs)-1, 2)]
			if geneObj.phen == 'undefined':
				geneObj.phen = dict([(row[1], funcPairs)])
			elif row[1] not in geneObj.phen:
				geneObj.phen[row[1]] = funcPairs
			else:
				throwError(filename, rowNum, 'Specified phenotype already exists for gene.')
	return genes		
			
			
#raises exception, displaying filename and row number
def throwError(filename, rowNum, message):
	print (filename + ', Row ' + str(rowNum) + ': ' + message)
	sys.exit(1)

#throws error if row is the wrong length.
#inputs: filename (string), row data (list), row number in file (int), expected row format (list), row delimiter (string).
def checkRowLen(filename, row, rowNum, format, delim):
	if len(row) != len(format):
		throwError(filename, rowNum, 'Number of values per row should be ' 
		+ str(len(format)) + ', is ' + str(len(row)) + '.\nExpected format: ' + delim.join(format))
		
#converts string with elements separated by delimRE, into list 
def processStrToList(str, delimRE):
	str = re.sub(delimRE, '$', str)
	if str[0] == '$':
		str = str[1:]
	if str[-1] == '$':
		str = str[:-1]
	return str.split('$')

#returns top n values in given list, in descending order.
def calcMaxNValuesInList(list1, n):
	list2 = list(list1)	
	maxVals=[]
	for i in range(min(len(list2), n)):				
		highestVal = max(list2)
		maxVals += [highestVal]
		list2.remove(highestVal)
	return maxVals
	
#calculates zygosity given a list of reference and variant allele read depths.
#returns string.
def calcZygosity(readDepthList):
	highestCount, secondHighestCount = calcMaxNValuesInList(readDepthList, 2)
	if highestCount != 0:
		if float(secondHighestCount) / highestCount >= 0.5:
			return 'HET'
		elif readDepthList[0] == highestCount:
			return 'HOM-REF'
		else:
			return 'HOM-ALT'
	return 'undefined'

#writes info to output file
def writeResults(outFile, refSite, readDepths):
	zygosity = calcZygosity(readDepths)
	outFile.write(refSite.pgxGene + '\t' + refSite.rsVariantID + '\t' +
	refSite.variantGenomicPos + '\t' + refSite.variantTranscriptPos +
	'\t' + zygosity + '\t')
	isFirstReadDepth = 1
	for readDepth in readDepths:				
		if readDepth in calcMaxNValuesInList(readDepths, 2):
			outFile.write(str(readDepth))
			if isFirstReadDepth:
				outFile.write(',')
			isFirstReadDepth = 0			
	outFile.write('\t' + str(refSite.associatedAlleles) + '\n')						

				
if __name__ == '__main__':
  main()
