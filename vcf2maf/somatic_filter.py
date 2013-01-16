## a script to filter out non-somatic (germ-line)
## mutations from a VCF file, based on the SOMATIC
## flag appearing in the INFO column.
############################################################

import argparse
from itertools import ifilter
import re
import sys


def somatic_filter(vcf_file, out_file):

  match_somatic = re.compile('somatic', re.IGNORECASE)

  ## pass through header lines and lines whose
  ## info column contains the SOMATIC flag
  for line in vcf_file:
    if line.startswith('#'):
      out_file.write(line)
    else:
      fields = line.split('\t')
      if any( ifilter( match_somatic.match, fields[7].split(';') ) ):
        out_file.write(line)


def main():
  parser = argparse.ArgumentParser(description='Filter out non-somatic mutations.')
  parser.add_argument('-v', '--verbose', action='store_true', help='Generate verbose output')
  parser.add_argument('--vcf', dest='vcf_filenames', nargs='+', help='Filename of VCF file(s) to read')
  parser.add_argument('-o', '--out', dest='out_filename', help='Name of output file')
  parser.add_argument('vcf_filenames2', metavar='VCF_FILENAME', nargs='*', help='Filename of VCF file(s) to read')

  args = parser.parse_args()

  ## VCF files can be specified by --vcf switch or as positional args
  if args.vcf_filenames2:
    if args.vcf_filenames is None:
      args.vcf_filenames = args.vcf_filenames2
    else:
      args.vcf_filenames.extend(args.vcf_filenames2)

  ## filter a single VCF file which might come from STDIN,
  ## to a single output file or to STDOUT
  if len(args.vcf_filenames) < 2:
    with open(args.vcf_filenames[0], 'r') if args.vcf_filenames else sys.stdin as vcf_file:
      with open(args.out_filename, 'w') if args.out_filename else sys.stdout as out_file:
        somatic_filter(vcf_file, out_file)

  ## for multiple VCF files, we'll output each filtered file XXX.vcf
  ## to an output file XXX.somatic.vcf
  else:
    for vcf_filename in args.vcf_filenames:
      out_filename = re.sub(r'.vcf$', r'.somatic.vcf', vcf_filename)
      with open(vcf_filename, 'r') as vcf_file:
        with open(out_file, 'w') as out_file:
          somatic_filter(vcf_file, out_file)


## call main method if this file is run as a script
if __name__ == '__main__':
    main()
