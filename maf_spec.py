## specification for MAF file transcribed from
## https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification
## 
## MAF version 2.3, accessed 2013-01-04
############################################################

maf_version = 2.3

maf_column_names = ['Hugo_Symbol',
                    'Entrez_Gene_Id',
                    'Center',
                    'NCBI_Build',
                    'Chromosome',
                    'Start_Position',
                    'End_Position',
                    'Strand',
                    'Variant_Classification',
                    'Variant_Type',
                    'Reference_Allele',
                    'Tumor_Seq_Allele1',
                    'Tumor_Seq_Allele2',
                    'dbSNP_RS',
                    'dbSNP_Val_Status',
                    'Tumor_Sample_Barcode',
                    'Matched_Norm_Sample_Barcode',
                    'Match_Norm_Seq_Allele1',
                    'Match_Norm_Seq_Allele2',
                    'Tumor_Validation_Allele1',
                    'Tumor_Validation_Allele2',
                    'Match_Norm_Validation_Allele1',
                    'Match_Norm_Validation_Allele2',
                    'Verification_Status',
                    'Validation_Status',
                    'Mutation_Status',
                    'Sequencing_Phase',
                    'Sequence_Source',
                    'Validation_Method',
                    'Score',
                    'BAM_File',
                    'Sequencer',
                    'Tumor_Sample_UUID',
                    'Matched_Norm_Sample_UUID']

maf_columns = {
'Hugo_Symbol' : {
  'index':1,
  'name':'Hugo_Symbol',
  'description':'HUGO symbol for the gene (HUGO symbols are always in all caps). If no gene exists within 3kb enter "Unknown". Source: http://genenames.org',
  'example':'EGFR',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Set or Unknown'
},
'Entrez_Gene_Id' : {
  'index':2,
  'name':'Entrez_Gene_Id',
  'description':'Entrez gene ID (an integer). If no gene exists within 3kb enter "0". Source: http://ncbi.nlm.nih.gov/sites/entrez?db=gene',
  'example':'1956',
  'case sensitive':'No',
  'null':'No',
  'enumerated':'Set'
},
'Center' : {
  'index':3,
  'name':'Center',
  'description':'Genome sequencing center reporting the variant. If multiple institutions report the same mutation separate list using semicolons.',
  'example':'hgsc.bcm.edu;genome.wustl.edu',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Set'
},
'NCBI_Build' : {
  'index':4,
  'name':'NCBI_Build',
  'description':'Any TGCA accepted genome identifier.  Can be string, integer or a float. (formerly: NCBI human genome build number.  Can be an integer or a float.)',
  'example':'hg18, hg19, GRCh37, GRCh37-lite, 36.1, 37, etc.',
  'case sensitive':'No',
  'null':'No',
  'enumerated':'Set'
},
'Chromosome' : {
  'index':5,
  'name':'Chromosome',
  'description':'Chromosome number without "chr" prefix that contains the gene.',
  'example':'X, Y, M, 1, 2, etc.',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Set'
},
'Start_Position' : {
  'index':6,
  'name':'Start_Position',
  'description':'Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate (1-based coordinate system).',
  'example':'999',
  'case sensitive':'No',
  'null':'No',
  'enumerated':'Set'
},
'End_Position' : {
  'index':7,
  'name':'End_Position',
  'description':'Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate (inclusive, 1-based coordinate system).',
  'example':'1000',
  'case sensitive':'No',
  'null':'No',
  'enumerated':'Set'
},
'Strand' : {
  'index':8,
  'name':'Strand',
  'description':'Genomic strand of the reported allele. Variants should always be reported on the positive genomic strand. (Currently, only the positive strand is an accepted value).',
  'example':'+',
  'case sensitive':'No',
  'null':'No',
  'enumerated':'+'
},
'Variant_Classification' : {
  'index':9,
  'name':'Variant_Classification',
  'description':'Translational effect of variant allele.',
  'example':'Missense_Mutation',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, 3'UTR, 3'Flank, 5'UTR, 5'Flank, IGR1 , Intron, RNA, Targeted_Region, De_novo_Start_InFrame, or De_novo_Start_OutOfFrame'
},
'Variant_Type' : {
  'index':10,
  'name':'Variant_Type',
  'description':'Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP but for 3 consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of 4 or more.',
  'example':'INS',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'SNP, DNP, TNP, ONP, INS, DEL, or Consolidated2'
},
'Reference_Allele' : {
  'index':11,
  'name':'Reference_Allele',
  'description':'The plus strand reference allele at this position. Include the sequence deleted for a deletion, or "-" for an insertion.',
  'example':'A',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'A,C,G,T, and/or -'
},
'Tumor_Seq_Allele1' : {
  'index':12,
  'name':'Tumor_Seq_Allele1',
  'description':'Primary data genotype. Tumor sequencing (discovery) allele 1. " -" for a deletion represent a variant. "-" for an insertion represents wild-type allele. Novel inserted sequence for insertion should not include flanking reference bases.',
  'example':'C',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'A,C,G,T, and/or -'
},
'Tumor_Seq_Allele2' : {
  'index':13,
  'name':'Tumor_Seq_Allele2',
  'description':'Primary data genotype. Tumor sequencing (discovery) allele 2. " -" for a deletion represents a variant. "-" for an insertion represents wild-type allele. Novel inserted sequence for insertion should not include flanking reference bases.',
  'example':'G',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'No'
},
'dbSNP_RS' : {
  'index':14,
  'name':'dbSNP_RS',
  'description':'Latest dbSNP rs ID (dbSNP_ID) or "novel" if there is no dbSNP record. source: http://ncbi.nlm.nih.gov/projects/SNP/',
  'example':'rs12345',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':'Set or "novel"'
},
'dbSNP_Val_Status' : {
  'index':15,
  'name':'dbSNP_Val_Status',
  'description':'dbSNP validation status. Semicolon- separated list of validation statuses.',
  'example':'by2Hit2Allele;byCluster',
  'case sensitive':'No',
  'null':'Yes',
  'enumerated':'by1000genomes;by2Hit2Allele; byCluster; byFrequency; byHapMap; byOtherPop; bySubmitter; alternate_allele3'
},
'Tumor_Sample_Barcode' : {
  'index':16,
  'name':'Tumor_Sample_Barcode',
  'description':'BCR aliquot barcode for the tumor sample including the two additional fields indicating plate and well position. i.e. TCGA-SiteID-PatientID-SampleID-PortionID-PlateID-CenterID. The full TCGA Aliquot ID.',
  'example':'TCGA-02-0021-01A-01D-0002-04',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Set'
},
'Matched_Norm_Sample_Barcode' : {
  'index':17,
  'name':'Matched_Norm_Sample_Barcode',
  'description':'BCR aliquot barcode for the matched normal sample including the two additional fields indicating plate and well position. i.e. TCGA-SiteID-PatientID-SampleID-PortionID-PlateID-CenterID. The full TCGA Aliquot ID; e.g. TCGA-02-0021-10A-01D-0002-04 (compare portion ID '10A' normal sample, to '01A' tumor sample).',
  'example':'TCGA-02-0021-10A-01D-0002-04',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Set'
},
'Match_Norm_Seq_Allele1' : {
  'index':18,
  'name':'Match_Norm_Seq_Allele1',
  'description':'Primary data. Matched normal sequencing allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'T',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':'A,C,G,T, and/or -'
},
'Match_Norm_Seq_Allele2' : {
  'index':19,
  'name':'Match_Norm_Seq_Allele2',
  'description':'Primary data. Matched normal sequencing allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'ACGT',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':'A,C,G,T, and/or -'
},
'Tumor_Validation_Allele1' : {
  'index':20,
  'name':'Tumor_Validation_Allele1',
  'description':'Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'-',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':'A,C,G,T, and/or -'
},
'Tumor_Validation_Allele2' : {
  'index':21,
  'name':'Tumor_Validation_Allele2',
  'description':'Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'A',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':'A,C,G,T, and/or -'
},
'Match_Norm_Validation_Allele1' : {
  'index':22,
  'name':'Match_Norm_Validation_Allele1',
  'description':'Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'C',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':'A,C,G,T, and/or -'
},
'Match_Norm_Validation_Allele2' : {
  'index':23,
  'name':'Match_Norm_Validation_Allele2',
  'description':'Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.',
  'example':'G',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':'A,C,G,T, and/or -'
},
'Verification_Status' : {
  'index':24,
  'name':'Verification_Status',
  'description':'Second pass results from independent attempt using same methods as primary data source. Generally reserved for 3730 Sanger Sequencing.',
  'example':'Verified',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':'Verified, Unknown'
},
'Validation_Status' : {
  'index':25,
  'name':'Validation_Status',
  'description':'Second pass results from orthogonal technology.',
  'example':'Valid',
  'case sensitive':'Yes',
  'null':'Yes',
  'enumerated':'Valid, Unknown, Wildtype'
},
'Mutation_Status' : {
  'index':26,
  'name':'Mutation_Status',
  'description':'Updated to reflect validation or verification status.',
  'example':'Somatic',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Somatic'
},
'Sequencing_Phase' : {
  'index':27,
  'name':'Sequencing_Phase',
  'description':'TCGA sequencing phase. Phase should change under any circumstance that the targets under consideration change.',
  'example':'Phase_I',
  'case sensitive':'No',
  'null':'Yes',
  'enumerated':'No'
},
'Sequence_Source' : {
  'index':28,
  'name':'Sequence_Source',
  'description':'Molecular assay type used to produce the analytes used for sequencing.',
  'example':'PCR;Capture',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'PCR, Capture, WGS'
},
'Validation_Method' : {
  'index':29,
  'name':'Validation_Method',
  'description':'The assay platforms used for the validation call. Examples: Sanger_PCR_WGA, Sanger_PCR_gDNA, 454_PCR_WGA, 454_PCR_gDNA; separate multiple entries using semicolons.',
  'example':'Sanger_PCR_WGA;Sanger_PCR_gDNA',
  'case sensitive':'No',
  'null':'Yes',
  'enumerated':'No'
},
'Score' : {
  'index':30,
  'name':'Score',
  'description':'Not in use.',
  'example':'NA',
  'case sensitive':'No',
  'null':'Yes',
  'enumerated':'No'
},
'BAM_File' : {
  'index':31,
  'name':'BAM_File',
  'description':'Not in use.',
  'example':'NA',
  'case sensitive':'No',
  'null':'Yes',
  'enumerated':'No'
},
'Sequencer' : {
  'index':32,
  'name':'Sequencer',
  'description':'Instrument used to produce primary data. Separate multiple entries using semicolons.',
  'example':'Illumina GAIIx;SOLID',
  'case sensitive':'Yes',
  'null':'No',
  'enumerated':'Illumina GAIIx, Illumina HiSeq, SOLID, 454, ABI 3730xl'
},
'Tumor_Sample_UUID' : {
  'index':33,
  'name':'Tumor_Sample_UUID',
  'description':'BCR aliquot UUID for tumor sample ',
  'example':'550e8400-e29b-41d4-a716-446655440000',
  'case sensitive':'Yes',
  'null':'No (for UUID-based files)',
  'enumerated':'No'
},
'Matched_Norm_Sample_UUID' : {
  'index':34,
  'name':'Matched_Norm_Sample_UUID',
  'description':'BCR aliquot UUID for matched normal',
  'example':'567e8487-e29b-32d4-a716-446655443246',
  'case sensitive':'Yes',
  'null':'No (for UUID-based files)',
  'enumerated':'No'
}}
