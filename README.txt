pgx-interpret.py  
10/31/14

Purpose
--------
This script processes sequencing runs intended to find which allelic variants are present
at known sites in the genome. 

For each variant found, the gene, rs number, genomic position, transcript position, 
zygosity, number of reads of each allele, and associated star alleles are determined. 
For each gene, the star allele and associated functional status for each chromosome, and 
the overall expected phenotype, are determined. 

Intermediate results are output to standard output and final results are output as a 
tab-separated text file.


Requirements
-------------
PyVCF (https://github.com/jamescasbon/PyVCF)
 

Input files
------------
-Reference vcf file: Specifies reference and alternate alleles, and associated star 
 alleles and genes, at desired locations in genome.
 
-Query vcf file: Contains alleles found at given genomic locations during a sequencing 
 run.
 
-Genes file: Specifies the name, default star allele, star allele functions, and 
 phenotypes for each gene. Formatted as a tab-separated text file.
 Example gene specification:
	>name	CYP2C19	
	defaultAllele	*1
	alleleFunc	nonFunc	Non-functional	*2, *3, *4, *5, *6, *7, *8
	alleleFunc	normFunc	Normal function	*1, *9, *10, *29, *41, *49, *50, *54, *55, *59, *69, *72,
	alleleFunc	gainFunc	Gain of function	*17
	phen	Ultrarapid metabolizer	gainFunc gainFunc, normFunc gainFunc
	phen	Extensive metabolizer	normFunc normFunc
	phen	Intermediate metabolizer	normFunc nonFunc
	phen	Poor metabolizer	nonFunc nonFunc
	
-Reference sites file: Tab-separated text file. In each row: star alleles, variant ID, 
 gene, rs number, genomic position, and transcript position.
 
 
Usage
-----
Argument: Path of query GENOTYPE_GIVEN_ALLELES mode (GATK Unified Genotyper) VCF file.
Options:
-v, --verbose: Verbose output.
--ref_vcf: Reference VCF filepath. Default: [pgx-interpret.py path] + "/pgx_genotyper_input.vcf".
--ref_site_filename: Reference site specification filepath. Default: [pgx-interpret.py path] + "/pgxRefSites.txt".
--gene_filename: Gene specification filepath. Default: [pgx-interpret.py path] + "/pgxGenes.txt".
--out_dir: Directory name for output files. Default: "pgx-interpret_out".
--discovery_gt_dir: Directory of DISCOVERY mode (GATK Unified Genotyper) VCF files.
-d, --depth: Read depth below which a low reads warning is printed.

Example commands: 
python ~/pgx-interpret.py 1330-01.knownsitesGT.vcf
(Assumes that directory containing pgx-interpret.py also contains pgx_genotyper_input.vcf,
pgxRefSites.txt, and pgxGenes.txt)

python pgx-interpret.py 1330-01.knownsitesGT.vcf -v --ref_vcf genotyper_input.vcf 
--ref_site_filename refsites.txt --gene_filename genes1.txt --discovery_gt_dir ../discovery_gt -d 25
 
 
Contact
-------
Diana Kolbe (diana-kolbe@uiowa.edu)
Julie Wertz (julie-wertz@uiowa.edu)
