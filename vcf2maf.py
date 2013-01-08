## vcf2maf
## chris.bare@sagebase.org
############################################################

import re
import sys
import argparse
import vcf
from time import time

## regex matching TCGA barcodes
_barcode_pattern = re.compile(r'(TCGA-\w{2}-\w{4}-\w{3}-\w{3}-\w{4}-\w{2})')

## genotypes may be delimited by | or slash characters
split_on_slash_or_bar = re.compile(r"[\|\/\\]")


def list_samples(vcf_file):
  """
  List the samples in a vcf file
  """
  vcf_reader = vcf.Reader(vcf_file)
  print "Sample IDs:"
  for sample in vcf_reader.metadata['SAMPLE']:
    print "  ",sample['ID']


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

def vcf2maf(vcf_file, maf_file, verbose=True):
  """
  Read a VCf file and output a corresponding MAF file.
  """

  vcf_reader = vcf.Reader(vcf_file)

  centers = vcf_reader.metadata['center']
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
  if verbose: sys.stderr.write('Normal sample found: %s\n' % (normal_sample_id,))

  ## look for TCGA barcode
  normal_sample_barcode = _extract_tcga_barcode(sample)
  if verbose: sys.stderr.write('normal sample barcode = '+str(normal_sample_barcode)+'\n')

  ## look for UUID
  normal_sample_uuid = sample['SampleUUID'] if 'SampleUUID' in sample else ''

  ## find tumor sample
  tumor_sample_names = ['tumor', 'primary']
  tumor_samples = [sample for sample in vcf_reader.metadata['SAMPLE'] if sample['ID'].lower() in tumor_sample_names]
  if len(tumor_samples)==0:
    raise Exception('No tumor samples found among: %s' % (", ".join(sample_names),))
  sample = tumor_samples[0]
  tumor_sample_id = sample['ID']
  if verbose: sys.stderr.write('Tumor sample found: %s\n' % (tumor_sample_id,))

  ## look for TCGA barcode
  tumor_sample_barcode = _extract_tcga_barcode(sample)
  if verbose: sys.stderr.write('normal sample barcode = '+str(tumor_sample_barcode)+'\n')

  ## look for UUID
  tumor_sample_uuid = sample['SampleUUID'] if 'SampleUUID' in sample else ''

  strand = '+'

  ## write Mutation Annotation Format (MAF) file

  try:
    i = 0
    for record in vcf_reader:
      i += 1
      if verbose and i % 10000 == 0:
        sys.stderr.write('processing record %d\n' % (i,))

      ## Hugo_Symbol
      maf_file.write('???')
      maf_file.write('\t')

      ## Entrez_Gene_Id
      ## TODO: look for INFO:GENE field
      maf_file.write('???')
      maf_file.write('\t')

      ## Center
      maf_file.write(';'.join(centers))
      maf_file.write('\t')

      ## NCBI_Build
      maf_file.write(ncbi_build)
      maf_file.write('\t')

      ## Chromosome
      maf_file.write(record.CHROM)
      maf_file.write('\t')

      ## Start_Position
      maf_file.write(str(record.start))
      maf_file.write('\t')

      ## End_Position
      maf_file.write(str(record.end))
      maf_file.write('\t')

      ## Strand
      maf_file.write(strand)
      maf_file.write('\t')

      ## Variant_Classification
      maf_file.write('???')
      maf_file.write('\t')

      ## Variant_Type
      maf_file.write(record.INFO['VT'])
      maf_file.write('\t')

      ## Reference_Allele
      maf_file.write(record.REF)
      maf_file.write('\t')

      ## extract tumor alleles
      record.genotype(tumor_sample_id)

      ## find tumor alleles
      tumor_alleles = split_on_slash_or_bar.split(record.genotype(tumor_sample_id).gt_bases)
      if len(tumor_alleles) < 1 or len(tumor_alleles) > 2:
        ## warning
        sys.stderr("%d tumor alleles detected in record %d. Only haploid and diploid genotypes are supported." % (len(tumor_alleles), i, ))

      ## Tumor_Seq_Allele1
      maf_file.write(tumor_alleles[0] if len(tumor_alleles) > 0 else '')
      maf_file.write('\t')

      ## Tumor_Seq_Allele2
      maf_file.write(tumor_alleles[1] if len(tumor_alleles) > 1 else '')
      maf_file.write('\t')

      ## dbSNP_RS
      maf_file.write('???')
      maf_file.write('\t')

      ## dbSNP_Val_Status
      maf_file.write('???')
      maf_file.write('\t')

      ## Tumor_Sample_Barcode
      maf_file.write(tumor_sample_barcode)
      maf_file.write('\t')

      ## Matched_Norm_Sample_Barcode
      maf_file.write(normal_sample_barcode)
      maf_file.write('\t')

      ## find normal alleles
      normal_alleles = split_on_slash_or_bar.split(record.genotype(normal_sample_id).gt_bases)
      if len(normal_alleles) < 1 or len(normal_alleles) > 2:
        ## warning
        sys.stderr("%d normal alleles detected in record %d. Only haploid and diploid genotypes are supported." % (len(normal_alleles), i, ))

      ## Match_Norm_Seq_Allele1
      maf_file.write(normal_alleles[0] if len(normal_alleles) > 0 else '')
      maf_file.write('\t')

      ## Match_Norm_Seq_Allele2
      maf_file.write(normal_alleles[1] if len(normal_alleles) > 1 else '')
      maf_file.write('\t')

      ## Tumor_Validation_Allele1
      maf_file.write('')
      maf_file.write('\t')

      ## Tumor_Validation_Allele2
      maf_file.write('')
      maf_file.write('\t')

      ## Match_Norm_Validation_Allele1
      maf_file.write('')
      maf_file.write('\t')

      ## Match_Norm_Validation_Allele2
      maf_file.write('')
      maf_file.write('\t')

      ## Verification_Status
      maf_file.write('')
      maf_file.write('\t')

      ## Validation_Status
      ## TODO: see: Including validation status in VCF file in
      ## https://wiki.nci.nih.gov/display/TCGA/TCGA+Variant+Call+Format+%28VCF%29+Specification
      maf_file.write('')
      maf_file.write('\t')

      ## Mutation_Status
      maf_file.write('Somatic')
      maf_file.write('\t')

      ## Sequencing_Phase (ex. Phase_I)
      ## TODO: how does this relate to phase in VCF, if at all?
      maf_file.write('')
      maf_file.write('\t')

      ## Sequence_Source
      maf_file.write('???')
      maf_file.write('\t')

      ## Validation_Method
      maf_file.write('???')
      maf_file.write('\t')

      ## Score
      maf_file.write('')
      maf_file.write('\t')

      ## BAM_File
      maf_file.write('')
      maf_file.write('\t')

      ## Sequencer
      maf_file.write(platform)
      maf_file.write('\t')

      ## Tumor_Sample_UUID
      maf_file.write(tumor_sample_uuid)
      maf_file.write('\t')

      ## Matched_Norm_Sample_UUID
      maf_file.write(normal_sample_uuid)

      maf_file.write('\n')

  except Exception:
    print "Error processing VCF record %d" % (i,)
    raise

  ## return the number of VCF records processed
  return i


def main():
  t = time()
  parser = argparse.ArgumentParser(description='Convert VCF files to MAF format.')
  parser.add_argument('-v', '--verbose', action='store_true', help='Generate verbose output')
  parser.add_argument('--vcf', dest='vcf_filename', help='Filename of a VCF file to read')
  parser.add_argument('--maf', dest='maf_filename', help='Filename of MAF file to output')
  parser.add_argument('-l', '--list-samples', action='store_true', help='List the sample IDs in a VCF file')

  args = parser.parse_args()

  if args.vcf_filename:
    vcf_file = open(args.vcf_filename, 'r')
  else:
    vcf_file = sys.stdin

  if args.list_samples:
    list_samples(vcf_file)
    return(0)

  if args.maf_filename:
    maf_file = open(args.maf_filename, 'w')
  else:
    maf_file = sys.stdout

  n = vcf2maf(vcf_file, maf_file)

  t = time() - t

  if args.verbose:
    sys.stderr.write("vcf2maf processed %d records in %0.1f seconds\n" % (n, t,))
    sys.stderr.write("...or %0.1f seconds per %d records\n" % (t/n*100000, 100000,))

293178

## call main method if this file is run as a script
if __name__ == "__main__":
    main()

