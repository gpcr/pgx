#!/usr/bin/python

import sys
import argparse
import vcf
import vcf.utils

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-v","--verbose")
  parser.add_argument("--ref_vcf", type=argparse.FileType('r'), default='pgx_genotyper_input.vcf')
  parser.add_argument("-d","--depth", type=int, default=30)
  parser.add_argument('vcf_file', nargs='+', type=argparse.FileType('r'))
  args = parser.parse_args()

  # Set up reference tables
  tagsites = list(vcf.Reader(args.ref_vcf))

  for samplefile in args.vcf_file:
    # Process a sample
    sample = vcf.Reader(samplefile)
    if len(sample.samples) > 1:
      print 'Interpretation of multi-sample VCF files not implemented'
      sys.exit(1)

    ID = sample.samples[0]
    print ID
    for (refrecord, record) in vcf.utils.walk_together(iter(tagsites),sample, kwargs=lambda r: (r.CHROM, r.POS, r.REF, R.ALT)):
      if refrecord is None:
         print 'ERROR: novel variant?', record
         sys.exit(1)
      if record is None:
         if args.verbose or 'AscAlleles' in refrecord.INFO: 
           print 'WARNING: tagsite missing', refrecord
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
        print spacer, ID, refrecord.INFO['Gene'], refrecord.INFO['AscAlleles'], call
        #print refrecord.ALT, record.ALT

if __name__ == '__main__':
  main()
