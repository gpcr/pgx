#!/usr/bin/python

import os
import sys
import argparse
import vcf
import vcf.utils
import pgxClasses
import re

def main():
	scriptPath = os.path.realpath(__file__)[:os.path.realpath(__file__).index(os.path.basename(__file__))]
	parser = argparse.ArgumentParser()
	parser.add_argument("-v","--verbose")
	parser.add_argument("--ref_vcf", type=argparse.FileType('r'), default=scriptPath+'/pgx_genotyper_input.vcf')
	parser.add_argument("--ref_site_filename", type=str, default = scriptPath+'/pgxRefSites.txt')
	parser.add_argument("--gene_filename", type=str, default=scriptPath+'/pgxGenes.txt')
	parser.add_argument("--out_filename", type=str, default='out.txt')
	parser.add_argument("--discovery_gt_dir", type=str, default='.')
	parser.add_argument("-d","--depth", type=int, default=30)
	parser.add_argument('vcf_file', nargs='+', type=argparse.FileType('r'))

	args = parser.parse_args()

	delimREGenesFile = r'[\s,;\/\\]+'
	delimREAssocAlleles = r'[^0-9A-Z*]+'
	genes = readInGenes(args.gene_filename, delimREGenesFile)
	allRefSites = readInRefSites(args.ref_site_filename, [z.name for z in genes], delimREAssocAlleles)
	finalRefSites = [] 	#ref sites to output
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
#		print ID
		for (refrecord, record) in vcf.utils.walk_together(iter(tagsites),sample, kwargs=lambda r: (r.CHROM, r.POS, r.REF, R.ALT)):
			if refrecord is None:
				print 'ERROR: novel variant?', record
				sys.exit(1)
			if record is None:
#				if args.verbose or 'AscAlleles' in refrecord.INFO: 
#					print 'WARNING: tagsite missing', refrecord
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
#				print spacer, ID, refrecord.INFO['Gene'], refrecord.INFO['AscAlleles'], call
				#print refrecord.ALT, record.ALT	

				#find associated alleles for just the variant allele with the most reads
				assocAllelesRecordSet = set([])
				GT = [int(z) for z in call['GT'].split('/')]	
				assocAllelesRecordSet = set.union(*[set(processStrToList(x, delimREAssocAlleles)) for x in refrecord.INFO['AscAlleles']])
				#find reference with the same set of *alleles as the record. 
				for i in range(len(allRefSites)):
					if (assocAllelesRecordSet == allRefSites[i].assocAllelesParsed) and (allRefSites[i].geneName == refrecord.INFO['Gene']):
						allRefSites[i].zygosity = calcZygosity(GT)
						if call.gt_type is 0:
							try:
								allRefSites[i].readDepths = call['AD']	
							except AttributeError:
								pass
						else:
							try:
								allRefSites[i].readDepths = [call['AD'][0], call['AD'][max(GT)]]	
							
							#look in discovery file for missing AD
							except AttributeError:				
								try:
									discoveryF = open(args.discovery_gt_dir + '/' + samplefile.name[:samplefile.name.index('knownsitesGT.vcf')] + 'discoveryGT.vcf')
									for discoveryRecord in vcf.Reader(discoveryF):
										if discoveryRecord.CHROM == record.CHROM and discoveryRecord.POS == record.POS:
											try:
												allRefSites[i].readDepths = [discoveryRecord.genotype(ID)['AD'][0], discoveryRecord.genotype(ID)['AD'][max(GT)]]	
											except AttributeError:
												pass
									discoveryF.close()		
								except IOError:
									pass
									
							finalRefSites += [allRefSites[i]]
							printRefSite(allRefSites[i], outFile)

						
	#output genotypes, phenotypes						 
	for gene in genes:
		gene.genotype = calcGenotype(gene, [z for z in finalRefSites if z.geneName == gene.name], allRefSites, outFile)
		gene.phenotype = calcPhenotype(gene)
	outFile.write('\nGene\t*Alleles\tAllele functional status\tExpected phenotype\n')
	for gene in genes:
		outFile.write(gene.name + '\t' + ('/'.join(gene.genotype) if isinstance(gene.genotype, list) else gene.genotype)
			 + '\t' + (' / '.join(gene.funcStatus) if isinstance(gene.funcStatus, list) else gene.funcStatus)
			 +'\t' + gene.phenotype + '\n')
	outFile.close()
	
#reads in tab-separated reference site file. in each row: associatedAlleles, cypVariantID, 
#gene, rsVariantID, variantGenomicPos, variantTranscriptPos.
#returns list of pgxClasses.pgxRefSite objects.
def readInRefSites(filename, geneNames, delimRE):
	fields = pgxClasses.pgxRefSite().__dict__.keys()
	refSites=[]
	f = open(filename)
	for rowNum, row in enumerate(f, start=1):
		row = row.strip().split('\t')
		if row==['']: 
			continue
		checkRowLen(filename, row, rowNum, fields[:-3], '\\t')
		refSite = pgxClasses.pgxRefSite(*row)
		if refSite.geneName not in geneNames:
			throwError(filename, rowNum, 'Unrecognized gene name: ' + refSite.geneName +
			'\nAcceptable values: '	+ ', '.join(geneNames) + '. Specify genes in gene file.')	
		refSite.assocAllelesParsed = set(processStrToList(refSite.associatedAlleles, delimRE))
		refSites += [refSite]
	f.close()
	return refSites


#reads in tab-separated genes file.
def readInGenes(filename, delimRE):
	fields = pgxClasses.pgxGene().__dict__.keys()
	genes=[]
	geneStartSymbol = '>'
	f = open(filename)
	for rowNum, row in enumerate(f, start=1):
		row = row.strip().split('\t')
		if row==['']:
			continue
		if row[0][0] == geneStartSymbol:
			row[0] = row[0][1:]
			genes += [pgxClasses.pgxGene()]
 		if row[0] not in fields:
			throwError(filename, rowNum, 'First value in row must be a geneName field.\nAcceptable values: '	+ ', '.join(fields))	
		#dict of geneName field : expected row format
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
	f.close()
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
	if str == '':
		return []
	str = re.sub(delimRE, '$', str)
	if len(str)>0 and str[0] == '$':
		str = str[1:]
	if len(str)>0 and str[-1] == '$':
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
	
#returns string.
def calcZygosity(GT):		
	if all(x == 0 for x in GT):
		return 'HOM-REF'
	elif all(x!= 0 for x in GT):
		return 'HOM-ALT'
	else:
		return 'HET'	

#returns pair [*a, *b] indicating genotype	
#operates on refSites pertaining to single gene
def calcGenotype(gene, refSites, allRefSites, outFile):
	for refSite in refSites:
		if refSite.geneName != gene.name:
			raise Exception('Ref site does not correspond to gene')
		if len(refSite.assocAllelesParsed) == 0:
			raise Exception('Ref site has no associated alleles')		
			
	if len(refSites) == 0:
		return [gene.defaultAllele, gene.defaultAllele]
	assocAllelesRemaining = findAllelesRemaining(refSites, allRefSites, outFile) #list of sets
	assocAllelesRemaining = [set([re.sub('[A-Z]', '', z) for z in w]) for w in assocAllelesRemaining]
	if all(len(assocAllelesRemaining[i]) == 1 for i in range(len(refSites))):
		uniqueAlleles = set.union(*assocAllelesRemaining)
		if all(x.zygosity == 'HET' for x in refSites):
			if len(uniqueAlleles) == 1:
				return sortAllelePair([list(uniqueAlleles)[0], gene.defaultAllele])
			elif len(uniqueAlleles) == 2:
				return sortAllelePair([list(uniqueAlleles)[0], list(uniqueAlleles)[1]])
		elif all(x.zygosity == 'HOM-ALT' for x in refSites):
			if len(uniqueAlleles) == 1:
				return sortAllelePair([list(uniqueAlleles)[0], list(uniqueAlleles)[0]])
	return 'Unknown genotype'
	
#sorts allele pair [*x, *y] by numeric value of x and y
def sortAllelePair(allelePair):
	if int(allelePair[0][1:]) < int(allelePair[1][1:]):	
		return [allelePair[0], allelePair[1]]
	else:
		return [allelePair[1], allelePair[0]]

#inputs: gene object, 2-element list of genotype *alleles. enters func status into gene object.
#returns phenotype.
def calcPhenotype(gene): 	
	funcIDs = [[y for y in gene.alleleFunc if x in gene.alleleFunc[y][1]] for x in gene.genotype]
	if not (len(funcIDs) == 2 and len(funcIDs[0]) == 1 and len(funcIDs[1]) == 1):
		return 'Unknown phenotype'
	gene.funcStatus = [gene.alleleFunc[funcIDs[0][0]][0], gene.alleleFunc[funcIDs[1][0]][0]]
	phens = [y for y in gene.phen if set([x[0] for x in funcIDs]) in [set(z) for z in gene.phen[y]]]
	if len(phens) != 1:
		return 'Unknown phenotype'
	return phens[0]
	
#returns HOM-REF reference sites that have an associated allele that parameter reference site has.
#(which can be eliminated when deciding between multiple associated *alleles)
def findHomRef(refSite, allRefSites, outFile):
	homRefSites = []
	assocAlleles1 = set([re.sub('[A-Z]', '', z) for z in refSite.assocAllelesParsed])
	for site in allRefSites:
		assocAlleles2 = set([re.sub('[A-Z]', '', z) for z in site.assocAllelesParsed])
		if assocAlleles1.intersection(assocAlleles2) and refSite.geneName == site.geneName and site.zygosity == 'HOM-REF':
			homRefSites += [site]
	for x in homRefSites:
		printRefSite(x, outFile)
	return homRefSites

#find alleles that are not HOM-REF at other sites
#(alleles remaining when deciding between multiple associated alleles)
def findAllelesRemaining(refSites, allRefSites, outFile):
	assocAllelesRemaining = [z.assocAllelesParsed for z in refSites]
	for i in range(len(refSites)):
		assocAllelesParsed = set([re.sub('[A-Z]', '', z) for z in refSites[i].assocAllelesParsed])
		if len(assocAllelesParsed)>1:
			homRef = findHomRef(refSites[i], allRefSites, outFile)
			if len(homRef)>0:
				assocAllelesRemaining[i] = set([z for z in refSites[i].assocAllelesParsed if z not in set.union(*[w.assocAllelesParsed for w in homRef])]) 
	return assocAllelesRemaining


def printRefSite(refSite, outFile):
	outFile.write(refSite.geneName + '\t' +  refSite.rsVariantID + '\t' + refSite.variantGenomicPos + 
	'\t' + refSite.variantTranscriptPos + '\t' + refSite.zygosity + '\t' 
	+ (','.join(str(x) for x in refSite.readDepths) if refSite.readDepths != 'undefined' else 'undefined')
	+ '\t' + refSite.associatedAlleles + '\n')		


			
if __name__ == '__main__':
  main()
