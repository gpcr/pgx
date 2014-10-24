class pgxRefSite(object):
	def __init__(self, associatedAlleles = 'undefined', cypVariantID = 'undefined', 
	pgxGene = 'undefined', rsVariantID = 'undefined', variantGenomicPos = 'undefined', 
	variantTranscriptPos = 'undefined'):
		self.associatedAlleles = associatedAlleles
		self.cypVariantID = cypVariantID
		self.pgxGene = pgxGene
		self.rsVariantID = rsVariantID
		self.variantGenomicPos = variantGenomicPos
		self.variantTranscriptPos = variantTranscriptPos
	
		
class pgxGene(object):						
	def __init__ (self, name= 'undefined', defaultAllele = 'undefined', 
	alleleFunc = 'undefined', phen = 'undefined'):
		self.name = name
		self.defaultAllele = defaultAllele
		self.alleleFunc = alleleFunc	#dict of allele function ID : [allele function label, set of associated *alleles]
		self.phen = phen	#dict of phenotype : list of corresponding allele function pairs
