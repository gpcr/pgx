class pgxRefSite(object):
	def __init__(self, associatedAlleles = 'undefined', cypVariantID = 'undefined', 
	geneName = 'undefined', rsVariantID = 'undefined', variantGenomicPos = 'undefined', 
	variantTranscriptPos = 'undefined', assocAllelesParsed = 'undefined',
	zygosity = 'undefined', readDepths = 'undefined'):
		self.associatedAlleles = associatedAlleles  #string
		self.cypVariantID = cypVariantID
		self.geneName = geneName	#string
		self.rsVariantID = rsVariantID
		self.variantGenomicPos = variantGenomicPos
		self.variantTranscriptPos = variantTranscriptPos
		self.assocAllelesParsed = assocAllelesParsed #set
		self.zygosity = zygosity	#string
		self.readDepths = readDepths  #list of read depths for each record corresponding to that set of assoc. alleles
		
class pgxGene(object):						
	def __init__ (self, name= 'undefined', defaultAllele = 'undefined', alleleFunc = 'undefined', 
	phen = 'undefined', genotype = 'undefined', funcStatus = 'undefined', phenotype = 'undefined'):
		self.name = name
		self.defaultAllele = defaultAllele
		self.alleleFunc = alleleFunc	#dict of allele function ID : [allele function label, set of associated *alleles]
		self.phen = phen	#dict of phenotype : list of corresponding allele function pairs
		self.genotype = genotype  #2-element list[*allele1, *allele2]
		self.funcStatus = funcStatus  #2-element list [*allele1 function, *allele2 function]
		self.phenotype = phenotype	#string
