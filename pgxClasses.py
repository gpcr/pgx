class pgxRefSite(object):
	def __init__(self, associatedAlleles, cypVariantID, gene, rsVariantID, 
	variantGenomicPos, variantTranscriptPos):
		self.associatedAlleles = associatedAlleles
		self.cypVariantID = cypVariantID
		self.gene = gene
		self.rsVariantID = rsVariantID
		self.variantGenomicPos = variantGenomicPos
		self.variantTranscriptPos = variantTranscriptPos
