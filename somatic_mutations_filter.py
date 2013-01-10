## Somatic mutations filter
##   A subclass of PyVCF's vcf.filters.Base that passes
##   only somatic mutations.
##
##   example usage:
##   vcf_filter.py --local-script somatic_mutations_filter.py --no-filtered TCGA-A6-2672_W_SOLiD.vcf somatic > tmp.vcf
##
##   We could do almost the same this with grep, with the caveat that it wouldn't
##   properly update the metadata, like so:
##   grep -E "(^#)|(SOMATIC)" TCGA-A6-2672_W_SOLiD.vcf
############################################################

import vcf.filters

class SomaticMutationsFilter(vcf.filters.Base):
  'Filter that passes only somatic mutations'

  ## ID used in the FILTER metadata
  name = 'somatic'

  @classmethod
  def customize_parser(self, parser):
    parser.add_argument('--somatic', action='store_true', help='Filter sites below this quality')

  def filter_name(self):
    return self.name

  def __init__(self, args):
    self.apply_filter = args.somatic

  ## for some knuckledheaded reason, the return values they want are:
  ##   None => the variant passes the filter
  ##   any other value => the variant does NOT pass the filter
  def __call__(self, record):
    return None if 'SOMATIC' in record.INFO else False
