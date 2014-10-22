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
  parser.add_argument("--out_filename", type=str, default='out.txt')
  parser.add_argument("-d","--depth", type=int, default=30)
  parser.add_argument('vcf_file', nargs='+', type=argparse.FileType('r'))
  args = parser.parse_args()

  refSites = readInRefSites(args.ref_site_filename)
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
    
    
    	#find zygosity and reference matching to record
    	if (call.gt_type is not 0):
    		assocAllelesRecordSet = processAssocAlleles(str(refrecord.INFO['AscAlleles']))
    		for refSite in refSites:
				assocAllelesRefSet = processAssocAlleles(refSite.associatedAlleles)
				if (assocAllelesRecordSet.intersection(assocAllelesRefSet) 
				and refSite.gene == refrecord.INFO['Gene'] and call['AD']):
					readDepthList = list(call['AD'])					
					highestCount = max(readDepthList)
					readDepthList.remove(highestCount)
					secondHighestCount = max(readDepthList)
					if highestCount != 0:
						if float(secondHighestCount) / highestCount >= 0.5:
							zygosity = 'HET'
						elif call['AD'][0] == highestCount:
							zygosity = 'HOM-REF'
						else:
							zygosity = 'HOM-ALT'
							
					#output results to file		
					outFile.write(refSite.gene + '\t' + refSite.rsVariantID + '\t' +
					refSite.variantGenomicPos + '\t' + refSite.variantTranscriptPos +
					'\t' + zygosity + '\t')
					isFirstReadDepth = 1
					for readDepth in call['AD']:
						if readDepth in [highestCount, secondHighestCount]:
							outFile.write(str(readDepth))
							if isFirstReadDepth:
								outFile.write(',')
							isFirstReadDepth = 0
					outFile.write('\t' + str(refSite.associatedAlleles) + '\n')
							
  outFile.close()

#reads in tab-separated reference site file. in each row: associatedAlleles, cypVariantID, 
#gene, rsVariantID, variantGenomicPos, variantTranscriptPos.
#returns list of pgxClasses.pgxRefSite objects.
def readInRefSites(refSiteFilename):
	numValuesPerRow = 6
	refSites=[]
	refSiteFile = open(refSiteFilename)
	for refSiteRowNum, row in enumerate(refSiteFile, start=1):
		row = row.strip().split('\t')
		if len(row) != numValuesPerRow:
			print ('Reference sites file, Row ' + str(refSiteRowNum) + ': ' +
			'Number of values per row should be ' + str(numValuesPerRow) + ', is ' +
			 str (len(row)) + '.\nExpected format: associatedAlleles\\tcypVariantID\\t' +
				'gene\\trsVariantID\\tvariantGenomicPos\\tvariantTranscriptPos')
			sys.exit(1)
		refSites += [pgxClasses.pgxRefSite(row[0], row[1], row[2], row[3], row[4], row[5])]
	return refSites
	refSiteFile.close()

#processes a string of associated alleles -- extracts *n values from string.
#returns set of *n values
def processAssocAlleles(assocAllelesStr):
	assocAlleles = re.sub(r'[^0-9A-Z*]+', '$', assocAllelesStr)
	if assocAlleles[0] == '$':
		assocAlleles = assocAlleles[1:]
	if assocAlleles[-1] == '$':
		assocAlleles = assocAlleles[:-1]
	return set(assocAlleles.split('$'))


if __name__ == '__main__':
  main()
