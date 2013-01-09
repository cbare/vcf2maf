## vcf2maf
## chris.bare@sagebase.org
############################################################

import re
import sys
import argparse
import vcf
from time import time

## output stream for logging
_log = sys.stderr

## regex matching TCGA barcodes
_barcode_pattern = re.compile(r'(TCGA-\w{2}-\w{4}-\w{3}-\w{3}-\w{4}-\w{2})')

## genotypes may be delimited by | or slash characters
split_on_slash_or_bar = re.compile(r"[\|\/\\]")


def list_samples(vcf_file):
  """
  List the samples in a vcf file
  """
  vcf_reader = vcf.Reader(vcf_file)
  print 'Sample IDs:'
  for sample in vcf_reader.metadata['SAMPLE']:
    print ' ',sample['ID']


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

  strand = '+'

  ## write Mutation Annotation Format (MAF) file

  try:
    i = 0
    for record in vcf_reader:
      i += 1
      if verbose and i % 10000 == 0:
        _log.write('processing record %d\n' % (i,))

      fields = []

      ## Hugo_Symbol
      fields.append('???')

      ## Entrez_Gene_Id
      ## TODO: look for INFO:GENE field
      fields.append('???')

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
      fields.append('???')

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
          _log.write('%d tumor alleles detected in record %d. Only haploid and diploid genotypes are supported.' % (len(tumor_alleles), i, ))
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
      fields.append('???')

      ## dbSNP_Val_Status
      fields.append('???')

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
      fields.append('???')

      ## Validation_Method
      fields.append('???')

      ## Score
      fields.append('')

      ## BAM_File
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

  n = vcf2maf(vcf_file, maf_file, verbose=args.verbose)

  t = time() - t

  if args.verbose:
    _log.write('vcf2maf processed %d records in %0.1f seconds\n' % (n, t,))
    _log.write('...or %0.1f seconds per %d records\n' % (t/n*100000, 100000,))


## call main method if this file is run as a script
if __name__ == '__main__':
    main()

