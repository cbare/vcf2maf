## specification for MAF file transcribed from
## https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification
## 
## MAF version 2.3, accessed 2013-01-04
############################################################

maf_version = 2.3

columns = [
{
  'name':'Hugo_Symbol',
  'description':'HUGO symbol for the gene (HUGO symbols are always in all caps). If no gene exists within 3kb enter "Unknown". Source: http://genenames.org',
  'example':'EGFR',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Set or Unknown'
},
{
  'name':'Entrez_Gene_Id',
  'description':'Entrez gene ID (an integer). If no gene exists within 3kb enter "0". Source: http://ncbi.nlm.nih.gov/sites/entrez?db=gene',
  'example':'1956',
  'case sensitive':'No',
  'null':'No',
  'enumerated':'Set'
},
{
  'name':'Center',
  'description':'Genome sequencing center reporting the variant. If multiple institutions report the same mutation separate list using semicolons.',
  'example':'hgsc.bcm.edu;genome.wustl.edu',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Set'
},
{
  'name':'NCBI_Build',
  'description':'Any TGCA accepted genome identifier.  Can be string, integer or a float. (formerly: NCBI human genome build number.  Can be an integer or a float.)',
  'example':'hg18, hg19, GRCh37, GRCh37-lite, 36.1, 37, etc.',
  'case sensitive':'No',
  'null':'No',
  'enumerated':'Set'
},
{
  'name':'Chromosome',
  'description':'Chromosome number without "chr" prefix that contains the gene.',
  'example':'X, Y, M, 1, 2, etc.',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Set'
},
{
  'name':'Start_Position',
  'description':'Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate (1-based coordinate system).',
  'example':'999',
  'case sensitive':'No',
  'null':'No',
  'enumerated':'Set'
},
{
  'name':'End_Position',
  'description':'Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate (inclusive, 1-based coordinate system).',
  'example':'1000',
  'case sensitive':'No',
  'null':'No',
  'enumerated':'Set'
},
{
  'name':'Strand',
  'description':'Genomic strand of the reported allele. Variants should always be reported on the positive genomic strand. (Currently, only the positive strand is an accepted value).',
  'example':'+',
  'case sensitive':'No',
  'null':'No',
  'enumerated':'+'
},
{
  'name':'Variant_Classification',
  'description':'Translational effect of variant allele.',
  'example':'Missense_Mutation',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Silent', 'Splice_Site', 'Translation_Start_Site', 'Nonstop_Mutation', '3\'UTR', '3\'Flank', '5\'UTR', '5\'Flank', 'IGR', 'Intron', 'RNA', 'Targeted_Region', 'De_novo_Start_InFrame', 'De_novo_Start_OutOfFrame']
},
{
  'name':'Variant_Type',
  'description':'Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP but for 3 consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of 4 or more.',
  'example':'INS',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':['SNP', 'DNP', 'TNP', 'ONP', 'INS', 'DEL', 'Consolidated']
},
{
  'name':'Reference_Allele',
  'description':'The plus strand reference allele at this position. Include the sequence deleted for a deletion, or "-" for an insertion.',
  'example':'A',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':['A','C','G','T','-']
},
{
  'name':'Tumor_Seq_Allele1',
  'description':'Primary data genotype. Tumor sequencing (discovery) allele 1. " -" for a deletion represent a variant. "-" for an insertion represents wild-type allele. Novel inserted sequence for insertion should not include flanking reference bases.',
  'example':'C',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':['A','C','G','T','-']
},
{
  'name':'Tumor_Seq_Allele2',
  'description':'Primary data genotype. Tumor sequencing (discovery) allele 2. " -" for a deletion represents a variant. "-" for an insertion represents wild-type allele. Novel inserted sequence for insertion should not include flanking reference bases.',
  'example':'G',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'No'
},
{
  'name':'dbSNP_RS',
  'description':'Latest dbSNP rs ID (dbSNP_ID) or "novel" if there is no dbSNP record. source: http://ncbi.nlm.nih.gov/projects/SNP/',
  'example':'rs12345',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':'Set or "novel"'
},
{
  'name':'dbSNP_Val_Status',
  'description':'dbSNP validation status. Semicolon- separated list of validation statuses.',
  'example':'by2Hit2Allele;byCluster',
  'case sensitive':'No',
  'null':'Yes',
  'enumerated':['by1000genomes', 'by2Hit2Allele', 'byCluster', 'byFrequency', 'byHapMap', 'byOtherPop', 'bySubmitter', 'alternate_allele3']
},
{
  'name':'Tumor_Sample_Barcode',
  'description':'BCR aliquot barcode for the tumor sample including the two additional fields indicating plate and well position. i.e. TCGA-SiteID-PatientID-SampleID-PortionID-PlateID-CenterID. The full TCGA Aliquot ID.',
  'example':'TCGA-02-0021-01A-01D-0002-04',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Set'
},
{
  'name':'Matched_Norm_Sample_Barcode',
  'description':'BCR aliquot barcode for the matched normal sample including the two additional fields indicating plate and well position. i.e. TCGA-SiteID-PatientID-SampleID-PortionID-PlateID-CenterID. The full TCGA Aliquot ID; e.g. TCGA-02-0021-10A-01D-0002-04 (compare portion ID \'10A\' normal sample, to \'01A\' tumor sample).',
  'example':'TCGA-02-0021-10A-01D-0002-04',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Set'
},
{
  'name':'Match_Norm_Seq_Allele1',
  'description':'Primary data. Matched normal sequencing allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'T',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':['A','C','G','T','-']
},
{
  'name':'Match_Norm_Seq_Allele2',
  'description':'Primary data. Matched normal sequencing allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'ACGT',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':['A','C','G','T','-']
},
{
  'name':'Tumor_Validation_Allele1',
  'description':'Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'-',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':['A','C','G','T','-']
},
{
  'name':'Tumor_Validation_Allele2',
  'description':'Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'A',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':['A','C','G','T','-']
},
{
  'name':'Match_Norm_Validation_Allele1',
  'description':'Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'C',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':['A','C','G','T','-']
},
{
  'name':'Match_Norm_Validation_Allele2',
  'description':'Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'G',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':['A','C','G','T','-']
},
{
  'name':'Verification_Status',
  'description':'Second pass results from independent attempt using same methods as primary data source. Generally reserved for 3730 Sanger Sequencing.',
  'example':'Verified',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':['Verified', 'Unknown']
},
{
  'name':'Validation_Status',
  'description':'Second pass results from orthogonal technology.',
  'example':'Valid',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':['Valid', 'Unknown', 'Wildtype']
},
{
  'name':'Mutation_Status',
  'description':'Updated to reflect validation or verification status.',
  'example':'Somatic',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Somatic'
},
{
  'name':'Sequencing_Phase',
  'description':'TCGA sequencing phase. Phase should change under any circumstance that the targets under consideration change.',
  'example':'Phase_I',
  'case sensitive':'No',
  'null':'Yes',
  'enumerated':'No'
},
{
  'name':'Sequence_Source',
  'description':'Molecular assay type used to produce the analytes used for sequencing.',
  'example':'PCR;Capture',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':['PCR', 'Capture', 'WGS']
},
{
  'name':'Validation_Method',
  'description':'The assay platforms used for the validation call. Examples: Sanger_PCR_WGA, Sanger_PCR_gDNA, 454_PCR_WGA, 454_PCR_gDNA; separate multiple entries using semicolons.',
  'example':'Sanger_PCR_WGA;Sanger_PCR_gDNA',
  'case sensitive':'No',
  'null':'Yes',
  'enumerated':'No'
},
{
  'name':'Score',
  'description':'Not in use.',
  'example':'NA',
  'case sensitive':'No',
  'null':'Yes',
  'enumerated':'No'
},
{
  'name':'BAM_File',
  'description':'Not in use.',
  'example':'NA',
  'case sensitive':'No',
  'null':'Yes',
  'enumerated':'No'
},
{
  ## the MAF 2.3 specification allows these values:
  ## Illumina GAIIx, Illumina HiSeq, SOLID, 454, ABI 3730xl

  ## But, according to the validator, Sequencer must be one or more of:
  ## Illumina GAIIx, Illumina HiSeq, SOLID, 454, ABI 3730xl,
  ## Ion Torrent PGM, Ion Torrent Proton, PacBio RS, Illumina MiSeq,
  ## Illumina HiSeq 2500, 454 GS FLX Titanium, AB SOLiD 4 System
  ## separate multiple entries with semicolon

  'name':'Sequencer',
  'description':'Instrument used to produce primary data. Separate multiple entries using semicolons.',
  'example':'Illumina GAIIx;SOLID',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':['Illumina GAIIx', 'Illumina HiSeq', 'SOLID', '454', 'ABI 3730xl', 
                'Ion Torrent PGM', 'Ion Torrent Proton', 'PacBio RS', 'Illumina MiSeq',
                'Illumina HiSeq 2500', '454 GS FLX Titanium', 'AB SOLiD 4 System']
},
{
  'name':'Tumor_Sample_UUID',
  'description':'BCR aliquot UUID for tumor sample ',
  'example':'550e8400-e29b-41d4-a716-446655440000',
  'case sensitive':'Yes',
  'null':'No (for UUID-based files)',
  'enumerated':'No'
},
{
  'name':'Matched_Norm_Sample_UUID',
  'description':'BCR aliquot UUID for matched normal',
  'example':'567e8487-e29b-32d4-a716-446655443246',
  'case sensitive':'Yes',
  'null':'No (for UUID-based files)',
  'enumerated':'No'
}]


columns_by_name = {}
for column in columns:
  columns_by_name[column['name']] = column


