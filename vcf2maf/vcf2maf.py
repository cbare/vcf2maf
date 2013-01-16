## vcf2maf
## chris.bare@sagebase.org
##
############################################################

## example usage:
##
## python vcf2maf.py -v --vcf examples/example1.vcf --maf out.maf
##
## python vcf2maf.py -v --list-samples --vcf examples/example1.vcf
##
## gunzip --stdout TCGA-AA-A024_W_IlluminaGA-DNASeq_exome.vcf.gz | python vcf2maf.py -v > out.maf
##

import re
import sys
import argparse
import vcf
from time import time
import maf_spec


## output stream for logging
_log = sys.stderr

## regex matching TCGA barcodes
_barcode_pattern = re.compile(r'(TCGA-\w{2}-\w{4}-\w{3}-\w{3}-\w{4}-\w{2})')

## regex matching dbSNP rsIDs
_rsid_pattern = re.compile(r'rs\d+')

## genotypes may be delimited by | or slash characters
split_on_slash_or_bar = re.compile(r"[\|\/\\]")


## snpEFF Functional class {NONE, SILENT, MISSENSE, NONSENSE}.
## snpEFF Effects:
## TODO: should be checked
snpeff_variant_classification = {
  "UPSTREAM" : "5'Flank",
  "DOWNSTREAM" : "3'Flank",
  "INTERGENIC" : "IGR",
  "INTERGENIC_CONSERVED" : "IGR",
  "INTRAGENIC" : "Silent",
  "UTR_5_PRIME" : "5'UTR",
  "UTR_5_DELETED" : "5'UTR",
  "UTR_3_PRIME" : "3'UTR",
  "UTR_3_DELETED" : "3'UTR",
  "START_GAINED" : "Translation_Start_Site",
  "START_LOST" : "Translation_Start_Site",
  "SYNONYMOUS_START" : "Translation_Start_Site",
  "NON_SYNONYMOUS_START" : "Translation_Start_Site",
  "STOP_GAINED" : "Missense_Mutation",
  "SYNONYMOUS_STOP" : "Silent",
  "STOP_LOST" : "Nonstop_Mutation",
  "SPLICE_SITE_ACCEPTOR" : "Splice_Site",
  "SPLICE_SITE_DONOR" : "Splice_Site",
  "NON_SYNONYMOUS_CODING" : None,
  "SYNONYMOUS_CODING" : "Silent",
  "FRAME_SHIFT" : None,
  "CODON_CHANGE" : "Missense_Mutation",
  "CODON_INSERTION" : "In_Frame_Ins",
  "CODON_CHANGE_PLUS_CODON_INSERTION" : "In_Frame_Ins",
  "CODON_DELETION" : "In_Frame_Del",
  "CODON_CHANGE_PLUS_CODON_DELETION" : "In_Frame_Del",
  "TRANSCRIPT" : None,
  "CDS" : None,
  "GENE" : None,
  "EXON" : None,
  "EXON_DELETED" : "In_Frame_Del",
  "INTRON" : "Intron",
  "INTRON_CONSERVED" : "Intron",
  "RARE_AMINO_ACID" : "Missense_Mutation" }

## snpEff fields
EFFECT = 0
FUNCLASS = 2

## EFF=FRAME_SHIFT & (INFO:VT=DEL or INFO:SVTYPE=DEL) => "Frame_Shift_Del",
## EFF=FRAME_SHIFT & (INFO:VT=INS or INFO:SVTYPE=INS) => "Frame_Shift_Ins",
## EFF=CODON_DELETION => "In_Frame_Del",
## EFF=CODON_INSERTION =>  "In_Frame_Ins",
## EFF=NON_SYNONYMOUS_CODING & Functional_class=MISSENSE => "Missense_Mutation",
## EFF=NON_SYNONYMOUS_CODING & Functional_class=NONSENSE => "Nonsense_Mutation",
## EFF=SYNONYMOUS_CODING => "Silent",
## EFF=SPLICE_SITE_ACCEPTOR or SPLICE_SITE_DONOR => "Splice_Site",
## EFF=START_GAINED or START_LOST or SYNONYMOUS_START(?) or NON_SYNONYMOUS_START => "Translation_Start_Site",
## EFF=STOP_LOST => "Nonstop_Mutation",
## EFF=UTR_3_PRIME or UTR_3_DELETED => "3'UTR",
## EFF=UPSTREAM(?) => "5'Flank",
## EFF=UTR_5_PRIME or UTR_5_DELETED => "5'UTR",
## EFF=DOWNSTREAM(?) => "3'Flank",
## EFF=INTERGENIC or INTERGENIC_CONSERVED => "IGR",
## EFF=INTRON or INTRON_CONSERVED => "Intron",
##  "RNA",
##  "Targeted_Region",
## EFF=START_GAINED ?? => "De_novo_Start_InFrame",
## EFF=START_GAINED ?? => "De_novo_Start_OutOfFrame"]

def snpeff_to_variant_classification(snpeff, record):
  """map snpEFF Effects to MAF variant classifications"""

  # TODO: this mapping is imperfect and should be looked at by
  # someone who understands the intent of the MAF classifications
  if snpeff[EFFECT]=='FRAME_SHIFT':
    if record.INFO.get('VT', None) == 'DEL' or record.INFO.get('SVTYPE', None) == 'DEL':
      variant_class = 'Frame_Shift_Del'
    elif record.INFO.get('VT', None) == 'INS' or record.INFO.get('SVTYPE', None) == 'INS':
      variant_class = 'Frame_Shift_Ins'
    else:
      variant_class = None
  if snpeff[EFFECT] in ['EXON', 'CDS', 'NON_SYNONYMOUS_CODING', 'TRANSCRIPT', 'GENE']:
    if snpeff[FUNCLASS]=='MISSENSE':
      variant_class = 'Missense_Mutation'
    elif snpeff[FUNCLASS]=='NONSENSE':
      variant_class = 'Nonsense_Mutation'
    else:
      variant_class = None
  else:
    variant_class = snpeff_variant_classification[snpeff[EFFECT]]

  return variant_class



def list_samples(vcf_file, out_file):
  """
  List the samples in a vcf file
  """
  vcf_reader = vcf.Reader(vcf_file)
  out_file.write('Sample IDs:\n')
  for sample in vcf_reader.metadata['SAMPLE']:
    out_file.write(' '+sample['ID']+'\n')


def _extract_tcga_barcode(sample):
  """
  Look for TCGA barcode in Sample metadata, under the key SampleTCGABarcode
  and then in the File field.
  """
  if 'SampleTCGABarcode' in sample:
    barcode = sample['SampleTCGABarcode']
  ## TODO: some vcf files have barcode in SampleName
  else:
    ## saw an example where the barcode was in the File field, like this:
    ## 'File': '"/inside/grotto/bambam/coad/exome/redo/TCGA-A6-2672-01A-01W-0833-10_SOLiD.bam"'
    m = _barcode_pattern.search(sample['File'])
    if m:
      barcode = m.group()
    else:
      barcode = ''
  return barcode


def vcf2maf(vcf_file, maf_file, verbose=False):
  """
  Read a VCf file and output a corresponding MAF file.
  """

  vcf_reader = vcf.Reader(vcf_file)

  centers = vcf_reader.metadata['center']

  ## get genome assembly identifier
  ## also look in ##contig=<ID={ID},length={length},assembly={assembly}
  ncbi_build = vcf_reader.metadata['reference']['ID'].split(' ')[0]

  platform = vcf_reader.metadata['SAMPLE'][0]['Platform']
  # TODO check all platforms and warn if different

  ## TCGA identifier for an individual
  individual = vcf_reader.metadata['INDIVIDUAL']

  ## In order to make a MAF file, we need to find normal and tumor samples.
  ## As far as I can tell, sample names in VCF are arbitrary, so we'll have
  ## to figure out by experience which samples map to normal and tumor.
  ## TODO need command line params to select samples
  sample_names = [sample['ID'] for sample in vcf_reader.metadata['SAMPLE']]

  ## find normal sample
  normal_samples = [sample for sample in vcf_reader.metadata['SAMPLE'] if sample['ID'].lower()=='normal']
  if len(normal_samples)==0:
    raise Exception('No normal samples found among: %s' % (", ".join(sample_names),))
  sample = normal_samples[0]
  normal_sample_id = sample['ID']
  if verbose: _log.write('Normal sample found: %s\n' % (normal_sample_id,))

  ## look for TCGA barcode
  normal_sample_barcode = _extract_tcga_barcode(sample)
  if verbose: _log.write('normal sample barcode = '+str(normal_sample_barcode)+'\n')

  ## look for UUID
  normal_sample_uuid = sample['SampleUUID'] if 'SampleUUID' in sample else ''

  ## find tumor sample
  tumor_sample_names = ['tumor', 'primary']
  tumor_samples = [sample for sample in vcf_reader.metadata['SAMPLE'] if sample['ID'].lower() in tumor_sample_names]
  if len(tumor_samples)==0:
    raise Exception('No tumor samples found among: %s' % (", ".join(sample_names),))
  sample = tumor_samples[0]
  tumor_sample_id = sample['ID']
  if verbose: _log.write('Tumor sample found: %s\n' % (tumor_sample_id,))

  ## look for TCGA barcode
  tumor_sample_barcode = _extract_tcga_barcode(sample)
  if verbose: _log.write('normal sample barcode = '+str(tumor_sample_barcode)+'\n')

  ## look for UUID
  tumor_sample_uuid = sample['SampleUUID'] if 'SampleUUID' in sample else ''

  ## always use plus strand by convention
  strand = '+'

  ## write Mutation Annotation Format (MAF) file header
  maf_file.write('\t'.join([col['name'] for col in maf_spec.columns]))
  maf_file.write('\n')

  try:
    i = 0
    for record in vcf_reader:
      i += 1
      if verbose and i % 10000 == 0:
        _log.write('processing record %d\n' % (i,))

      fields = []

      ## parse out snpEff fields:

      ## NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Atc/Gtc|I300V|749|NOC2L|protein_coding|CODING|ENST00000327044|)

      ## EFFECT
      ## IMPACT
      ## FUNCLASS
      ## CODON
      ## AA
      ## AA_LEN
      ## GENE
      ## BIOTYPE
      ## CODING
      ## TRID
      ## EXID

      variant_class = None
      hugo = None
      gene_id = None
      if 'EFF' in record.INFO:
        snpeffs = record.INFO['EFF'].split(',')

        # TODO: which snpeff to pick?
        if len(snpeffs) > 0:
          snpeff = re.split(r'[\(\)\|]', snpeffs[0])
          variant_class = snpeff_to_variant_classification(snpeff, record)
          #if not variant_class:
          #  _log.write('Unclassifiable variant: ' + str(snpeff) + '\n')
          hugo = snpeff[6]
          gene_id = snpeff[9]

      ## Hugo_Symbol
      # TODO: fake data
      # TODO: how to find HUGO gene name given position
      if 'HUGO' in record.INFO:
        hugo = record.INFO['HUGO']
      if not hugo:
        hugo = 'unknown'
      fields.append(hugo)

      ## Entrez_Gene_Id
      ## see Table 8a: Annotation fields added for RNA-Seq variants
      ## of the TCGA's VCF spec:
      ## https://wiki.nci.nih.gov/display/TCGA/TCGA+Variant+Call+Format+%28VCF%29+Specification
      # TODO: fake data:
      # TODO: where to find this, if not in optional GENE field?
      if 'GENE' in record.INFO:
        gene_id = record.INFO['GENE']
      if not gene_id:
        gene_id = 'unknown'
      fields.append(gene_id)

      ## Center
      fields.append(';'.join(centers))

      ## NCBI_Build
      fields.append(ncbi_build)

      ## Chromosome
      fields.append(record.CHROM)

      ## Start_Position
      fields.append(str(record.start))

      ## End_Position
      fields.append(str(record.end))

      ## Strand
      fields.append(strand)

      ## Variant_Classification
      if not variant_class:
        variant_class = 'unknown'
      fields.append(variant_class)

      ## Variant_Type
      fields.append(record.INFO['VT'])

      ## Reference_Allele
      fields.append(record.REF)

      ## extract tumor alleles
      record.genotype(tumor_sample_id)

      ## find tumor alleles
      ## note: gt_bases will return None if the VCF file 
      ##       see lines 389-391 of PyVCF/vcf/parser.py
      if record.genotype(tumor_sample_id).gt_bases:
        tumor_alleles = split_on_slash_or_bar.split(record.genotype(tumor_sample_id).gt_bases)
        if len(tumor_alleles) < 1 or len(tumor_alleles) > 2:
          ## warning
          _log.write('%d tumor alleles detected in record %d. Only haploid and diploid genotypes are supported.\n' % (len(tumor_alleles), i, ))
      else:
        ## warning
        _log.write('No tumor alleles detected. Skipping record %d.\n' % (i,))
        ## not sure what the right thing to do is here. MAF doesn't seem to
        ## allow blank alleles.
        continue

      ## Tumor_Seq_Allele1
      fields.append(tumor_alleles[0] if len(tumor_alleles) > 0 else '')

      ## Tumor_Seq_Allele2
      fields.append(tumor_alleles[1] if len(tumor_alleles) > 1 else '')

      ## dbSNP_RS
      if record.ID:
        rsids = [rec_id for rec_id in record.ID.split(';') if _rsid_pattern.match(rec_id)]
        if rsids:
          fields.append(rsids[0])
        else:
          fields.append('')

      ## dbSNP_Val_Status
      # TODO: look up dbSNP val status
      fields.append('')

      ## Tumor_Sample_Barcode
      fields.append(tumor_sample_barcode)

      ## Matched_Norm_Sample_Barcode
      fields.append(normal_sample_barcode)

      ## find normal alleles
      if record.genotype(normal_sample_id).gt_bases:
        normal_alleles = split_on_slash_or_bar.split(record.genotype(normal_sample_id).gt_bases)
        if len(normal_alleles) < 1 or len(normal_alleles) > 2:
          ## warning
          _log.write('%d normal alleles detected in record %d. Only haploid and diploid genotypes are supported.\n' % (len(normal_alleles), i, ))
      else:
        ## warning
        _log.write('No normal alleles detected. Skipping record %d.\n' % (i,))
        ## not sure what the right thing to do is here. MAF doesn't seem to
        ## allow blank alleles.
        continue

      ## Match_Norm_Seq_Allele1
      fields.append(normal_alleles[0] if len(normal_alleles) > 0 else '')

      ## Match_Norm_Seq_Allele2
      fields.append(normal_alleles[1] if len(normal_alleles) > 1 else '')

      ## Tumor_Validation_Allele1
      fields.append('')

      ## Tumor_Validation_Allele2
      fields.append('')

      ## Match_Norm_Validation_Allele1
      fields.append('')

      ## Match_Norm_Validation_Allele2
      fields.append('')

      ## Verification_Status
      fields.append('')

      ## Validation_Status
      ## TODO: see: Including validation status in VCF file in
      ## https://wiki.nci.nih.gov/display/TCGA/TCGA+Variant+Call+Format+%28VCF%29+Specification
      fields.append('')

      ## Mutation_Status
      fields.append('Somatic')

      ## Sequencing_Phase (ex. Phase_I)
      ## TODO: how does this relate to phase in VCF, if at all?
      fields.append('')

      ## Sequence_Source
      # TODO: fake data
      fields.append('PCR')

      ## Validation_Method
      # TODO
      fields.append('')

      ## Score (not in use)
      fields.append('')

      ## BAM_File (not in use)
      fields.append('')

      ## Sequencer
      fields.append(platform)

      ## Tumor_Sample_UUID
      fields.append(tumor_sample_uuid)

      ## Matched_Norm_Sample_UUID
      fields.append(normal_sample_uuid)

      ## write a row of the MAF file
      for field in fields:
        maf_file.write(field)
        maf_file.write('\t')
      maf_file.write('\n')

  except Exception:
    print "Error processing VCF record %d" % (i,)
    raise

  ## return the number of VCF records processed
  return i


def main():
  parser = argparse.ArgumentParser(description='Convert VCF files to MAF format.')
  parser.add_argument('-v', '--verbose', action='store_true', help='Generate verbose output')
  parser.add_argument('--vcf', dest='vcf_filenames', nargs='+', help='Filename of VCF file(s) to read')
  parser.add_argument('--maf', dest='maf_filename', help='Filename of MAF file to output')
  parser.add_argument('-o', '--out', dest='out_filename', help='Name of output file')
  parser.add_argument('-l', '--list-samples', action='store_true', help='List the sample IDs in a VCF file')
  parser.add_argument('vcf_filenames2', metavar='VCF_FILENAME', nargs='*', help='Filename of a VCF file to read')
  
  args = parser.parse_args()

  ## VCF files can be specified by --vcf switch or as positional args
  if args.vcf_filenames2:
    if args.vcf_filenames is None:
      args.vcf_filenames = args.vcf_filenames2
    else:
      args.vcf_filenames.extend(args.vcf_filenames2)

  ## read from VCF files or STDIN
  if args.vcf_filenames:
    vcf_files = ( open(filename, 'r') for filename in args.vcf_filenames )
  else:
    vcf_files = [ sys.stdin ]

  if args.list_samples:
    ## list samples in VCF file
    with open(args.out_filename, 'w') if args.out_filename else sys.stdout as out_file:
      for vcf_file in vcf_files:
        if args.verbose:
          _log.write('processing file: ' + vcf_file.name + '\n')
        list_samples(vcf_file, out_file)
        vcf_file.close()

  else:
    ## convert vcf(s) to maf
    if args.maf_filename:
      maf_file = open(args.maf_filename, 'w')
    elif args.out_filename:
      maf_file = open(args.out_filename, 'w')
    else:
      maf_file = sys.stdout

    t = time()
    n = 0

    for vcf_file in vcf_files:
      if args.verbose:
        _log.write('processing file: ' + vcf_file.name + '\n')

      n += vcf2maf(vcf_file, maf_file, verbose=args.verbose)

      vcf_file.close()
    maf_file.close()

    t = time() - t
    if args.verbose:
      _log.write('vcf2maf processed %d records in %0.1f seconds\n' % (n, t,))
      _log.write('...or %0.1f seconds per %d records\n' % (t/n*100000, 100000,))


## call main method if this file is run as a script
if __name__ == '__main__':
    main()
