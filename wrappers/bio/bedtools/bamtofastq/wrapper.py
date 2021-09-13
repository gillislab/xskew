__author__ = "John Hover"
__copyright__ = "Copyright 2021, John Hover"
__email__ = "jhover@pobox.com, hover@cshl.edu"
__license__ = "MIT"

from snakemake.shell import shell

## Extract arguments
extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

shell( " bedtools bamtofastq -i {snakemake.input} -fq {snakemake.output} " )
