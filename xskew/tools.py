import datetime as dt
import logging
import os
import pandas as pd
import shutil
import subprocess
import sys
import traceback

gitpath = os.path.expanduser("~/git/cshlwork")
sys.path.append(gitpath)

from cshlwork.utils import *

STARSUBDIRS = ['_STARgenome', '_STARpass1','_STARtmp' ]


class NonZeroReturnException(Exception):
    """
    Thrown when a command has non-zero return code. 
    """

def setup_logging(level):
    """ 
    Setup logging 

    """
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(name)s %(filename)s:%(lineno)d %(funcName)s(): %(message)s'
    logging.basicConfig()
    logger = logging.getLogger()
    streamHandler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(FORMAT)
    streamHandler.setFormatter(formatter)
    logger.addHandler(streamHandler)
    logger.setLevel(level)

def string_modulo(instring, divisor):
    """
    Takes instring. Converts to bytes. Takes hex() value of bytes. converts to integer. 
    returns final integer % modbase
    
    """
    encoded = instring.encode('utf-8')
    hstring = encoded.hex()
    intval = int(hstring, 16)
    return intval % divisor


def modulo_filter(inlist, divisor, remainder):
    """
    Takes a list, returns list containing items in inlist that 
    have the given remainder modulo divisor. 
    """
    newlist = []
    for e in inlist:
        if string_modulo(e, divisor) == remainder:
            newlist.append(e)
    logging.debug(f'inlist len={len(inlist)}, {divisor} servers, {remainder} server idx. outlist len={len(newlist)}')
    return newlist


def list_sample(infile):
    """
    log a directory listing of all files starting with <sample>. for infile. 
    e.g. /home/hover/work/xskew1/SRR1480384.X.rg.bam
    log a directory listing, with sizes, for all files starting with SRR1480384 in /home/hover/work/xskew1/
    
    """
    listing = ""    
    filename= os.path.basename(infile)
    dirname=os.path.dirname(infile)
    sampleid= filename.split('.')[0]
    logging.debug(f'listing {dirname}')
    dirlist = os.listdir(dirname)
    for x in dirlist:
        if x.startswith(sampleid):
            st = os.stat(f'{dirname}/{x}')
            filesize = st.st_size
            listing += f"{x}\t{filesize}\n"
    logging.debug(f"directory listing for sample {sampleid}:\n{listing}")


        
def fasterq_dump(infile, outdir, nthreads, tempdir, outfile=None ):
    '''
    for single-ended reads, output is sometimes <sample>.fastq and 
    sometimes <sample>_1.fastq
    so we need to fix filename after run. 
    
    '''
    
    filepath = os.path.abspath(infile)    
    dirname = os.path.dirname(filepath)
    filename = os.path.basename(filepath)
    (sample, ext) = os.path.splitext(filename)
    logging.debug(f'handling sample={sample} dirname={dirname}')
    
    cmd = ['fasterq-dump', 
    '--split-files',
    '--include-technical',
    '--force', 
    '--threads', nthreads,
    '--outdir', outdir, 
    '-t', tempdir,
    '--log-level', 'debug', 
    infile]
    try:
        cp = run_command(cmd)
        logging.info(f'got rc={cp.returncode}')
        if outfile is not None:
            logging.info(f'outfile={outfile}')
            # implies single-ended output
            logging.info(f'checking for {outdir}/{sample}.fastq')
            if os.path.exists(f'{outdir}/{sample}.fastq'):
                os.rename(f'{outdir}/{sample}.fastq', outfile)
                logging.info(f'renamed output {outdir}/{sample}.fastq -> {outfile} ')
            logging.info(f'checking for {outdir}/{sample}_1.fastq')
            if os.path.exists(f'{outdir}/{sample}_1.fastq'):
                os.rename(f'{outdir}/{sample}_1.fastq',outfile)            
                logging.info(f'renamed output {outdir}/{sample}_1.fastq -> {outfile} ')        

    #except NonZeroReturnException as nzre:
    #    logging.error(f'problem with {infile}')
    #    logging.error(traceback.format_exc(None))
    except Exception as e:
        logging.error(traceback.format_exc(None))
        
        
        
    
#def star_nowasp(end1, end2, outprefix, outtemp, nthreads, genomedir):

def parse_assembly_report(reportfile):
    '''
    parse genbank assembly report. 
            Sequence-Name Sequence-Role Assigned-Molecule  Assigned-Molecule-Location/Type  
            GenBank-Accn    Relationship    RefSeq-Accn    Assembly-Unit    Sequence-Length    
            UCSC-style-name
            
     
    
    '''


def star_genome(genomedir, nthreads, gtffile, infile ): 
    cmd = ['STAR',
           '--runMode', 'genomeGenerate',
           '--runThreadN', nthreads,
           '--genomeDir', genomedir,
           '--sjdbGTFfile', gtffile, 
           #'--sjdbOverhang','100',
           '--genomeFastaFiles', infile  
       ]
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {end1}/{end2} input files.')
        logging.error(traceback.format_exc(None))
        raise    


def make_chr_label(reportfile, outfile, chr='chrX' ):
    """
    reads NCBI assembly report and extracts selected scaffold/assembly label. 
    """
    if os.path.exists(reportfile):
        colnames = ['Sequence-Name','Sequence-Role','Assigned-Molecule',
                    'Assigned-Molecule-Location/Type','GenBank-Accn',
                    'Relationship','RefSeq-Accn','Assembly-Unit',
                    'Sequence-Length','UCSC-style-name']   
        df = pd.read_csv(reportfile, comment="#", sep='\t')
        df.columns = colnames
        if chr.startswith('chr'):
            chrnum = chr[3:]
        else:
            chrnum = chr
        tagval = chrnum
        label = df[ df['Assigned-Molecule'] == tagval]['RefSeq-Accn'].values[0]
        logging.debug(f'extracted label {label} for {tagval} in {reportfile}')
    else:
        label = chr

    f = open(outfile, 'w')
    f.write(f'{label}\n')
    f.close()


def make_chr_index(infile, genomedir, chr, outfile):
    region = get_chr_label(genomedir, chr)
    samtools_faidx_region(infile, outfile, region)  


def get_chr_label(genomedir, chr='chrX'):
    labelfile = f"{genomedir}/{chr}label.txt"
    f = open(labelfile, 'r')
    label = f.read().strip()
    logging.debug(f"retrieved label {label} for {chr} in {labelfile}")
    return label

def get_label(regionfile):
    f = open(regionfile, 'r')
    label = f.read().strip()
    logging.debug(f"retrieved label {label} for {chr} in {regionfile}")
    return label


def bedtools_bamtofastq(infile, end1, end2):
    cmd = ['bedtools',
           'bamtofastq',
           '-i', infile, 
           '-fq', end1,
           '-fq2', end2 
       ]
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise       


def star_nowasp(end1, outprefix, nthreads, genomedir, end2=None):
    
    if end2 is None:
        end2 = ''
            
    cmd = ['STAR',
       '--readFilesIn', end1, end2, 
       '--outFileNamePrefix', outprefix,
#       '--outTmpDir', outtemp, 
       '--runThreadN', nthreads,
       '--genomeDir', genomedir, 
       '--twopassMode Basic',
       '--twopass1readsN','-1',
       '--outSAMtype BAM Unsorted', 
       '--quantMode GeneCounts'
       ]
    try:
        cp = run_command(cmd)
        logging.info(f'got rc={cp.returncode}')
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {end1}/{end2} input files.')
        logging.error(traceback.format_exc(None))
    except Exception as e:
        logging.error(traceback.format_exc(None))
    
    finally:
        for ext in STARSUBDIRS:
            dirpath = f'{outprefix}{ext}'
            if os.path.exists(dirpath) and os.path.isdir(dirpath):
                logging.debug(f'removing temp STAR dir {dirpath} ...')
                shutil.rmtree(dirpath)    

#def star_wasp(end1, end2, vcf, outprefix, outtemp, nthreads, genomedir):
def star_wasp(end1, vcf, outprefix, nthreads, genomedir, end2=None):
    " STAR  --genomeDir {params.gdir} --readFilesIn {input.end1} {input.end2} "
    " --runThreadN {threads} --twopassMode Basic --twopass1readsN -1 " 
    " --outSAMtype BAM Unsorted --quantMode GeneCounts "
    " --waspOutputMode SAMtag --varVCFfile {input.sfvcf} && "            
    
    if end2 is None:
        end2 = ''
    
    cmd = ['STAR',
       '--readFilesIn', end1, end2, 
       '--outFileNamePrefix', outprefix,
 #      '--outTmpDir', outtemp, 
       '--runThreadN', nthreads,
       '--genomeDir', genomedir, 
       '--twopassMode Basic',
       '--twopass1readsN -1',
       '--outSAMtype BAM Unsorted', 
       '--quantMode GeneCounts',
       '--waspOutputMode', 'SAMtag',
       '--varVCFfile', vcf 
       ]
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {end1}/{end2} input files.')
        logging.error(traceback.format_exc(None))
    
    except Exception as e:
        logging.error(traceback.format_exc(None))   
    
    finally:
        for ext in STARSUBDIRS:
            dirpath = f'{outprefix}{ext}'
            if os.path.exists(dirpath) and os.path.isdir(dirpath):
                shutil.rmtree(dirpath)


def samtools_faidx(infile, outfile):
    cmd = ['samtools',
           'faidx',
           '-o', outfile,
           infile, 
       ]
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise            

def samtools_faidx_region(infile, outfile, region):
    cmd = ['samtools',
           'faidx',
           '-o', outfile,
           infile,
           region 
       ]
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise   

        
def samtools_dict(infile, outfile):
    cmd = ['samtools',
           'dict',
           '-o', outfile, 
           infile, 
       ]
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))    
        raise    

def samtools_sort_readname(infile, outfile, memory, nthreads):
    
    nthreads = str(nthreads)
    mstring = f'{memory}M'

    oldcmd = ['samtools',
           'sort',
           '-n',
           '-m', f'{memory}M',
           '-o', outfile, 
           '-O', 'bam', 
           '-@', f'{nthreads}',
           infile
       ]
    
    cmd = ['samtools',
           'sort',
           '-n',
           '-m', mstring,
           '-o', outfile, 
           '-O', 'bam', 
           '-@', nthreads,
           infile
       ]
    try:
        run_command_shell(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise 


def samtools_sort(infile, outfile, memory, nthreads):
    cmd = ['samtools',
           'sort',
           '-m', f'{memory}M',
           '-o', outfile, 
           '-O', 'bam', 
           '-@', f'{nthreads}',
           infile
       ]
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise    

def samtools_index(infile, nthreads):
    cmd = ['samtools',
           'index',
           '-@', f'{nthreads}',
           infile,
       ]
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise    


def samtools_view_region(infile, outfile, regionfile):
    region = get_label(regionfile)
    cmd = ['samtools',
           'view',
           '-b', infile,
           '-o', outfile,
           region
       ]
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))    
        raise        
    
def samtools_view_quality(infile, outfile, quality, tag=None):
    cmd = ['samtools',
           'view',
           '-q', quality, 
           '-b', infile,
           '-o', outfile,
       ]
    if tag is not None:
        cmd.append('-d')
        cmd.append(tag)
    
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))        
        raise    


def gatk_arrg(infile, outfile):
    " gatk AddOrReplaceReadGroups -I={input.xfiltbam} -O={params.chrom}.rg.bam -SO=coordinate "
    "-RGID=id -RGLB=library -RGPL=platform -RGPU=machine -RGSM=sample && "
    cmd = [ 'gatk',
            'AddOrReplaceReadGroups',
            '-I', infile,
            '-O' , outfile,
            '-SO', 'coordinate',
            '-RGID', 'id',
            '-RGLB', 'library',
            '-RGPL', 'platform',
            '-RGPU', 'machine',
            '-RGSM', 'sample'
        ]
    try:
        run_command_shell(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))    
        raise

def gatk_md(infile, outfile, metrics):
    " gatk MarkDuplicates -I={params.chrom}.rg.bam -O={params.chrom}.dedupped.bam -CREATE_INDEX=true " 
    " -VALIDATION_STRINGENCY=SILENT -M=output.metrics && "
    
    cmd = [ 'gatk',
           'MarkDuplicates',
           '-I', infile, 
           '-O', outfile,
           '-CREATE_INDEX', 'true',
           '-VALIDATION_STRINGENCY', 'SILENT',
           '-M', metrics,
           ]
    try:
        run_command_shell(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise    

def gatk_sncr(infile, outfile, genome ):
    " gatk SplitNCigarReads -R {params.gdir}/GRCh38.p7.genome.fa -I {params.chrom}.dedupped.bam " 
    " -O {params.chrom}.split.filtered.bam --java-options '-XXgcThreads:2 -XX:ConcGCThreads ' && "
    cmd = [ 'gatk', 'SplitNCigarReads',
           '-R', genome,
           '-I', infile, 
           '-O', outfile, 
           #'--java-options', "'-XXgcThreads:2 -XX:ConcGCThreads'" 
           ]
    try:
        run_command_shell(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))           
        raise    
    
def gatk_htc(infile, outfile, genome, regionfile):
    " gatk HaplotypeCaller -R {params.gdir}/GRCh38.p7.genome.fa -L chr{params.chrom} " 
    " -I {params.chrom}.split.filtered.bam --dont-use-soft-clipped-bases -stand-call-conf 0.0 "
    " -O {params.chrom}.filtered.vcf && "    
    interval = get_label(regionfile)
    
    cmd = [ 'gatk',
           'HaplotypeCaller',
           '-L', interval, 
           '-R', genome,
           '-I', infile, 
           '-O', outfile,  
           '--dont-use-soft-clipped-bases',
           '-stand-call-conf', '0.0',
           ]
    try:
        run_command_shell(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise    


def gatk_sv(infile, outfile, genome, regionfile):
    " gatk SelectVariants -R {params.gdir}/GRCh38.p7.genome.fa -L chr{params.chrom} "
    "-V {params.chrom}.filtered.vcf -O {params.chrom}.snps.vcf -select-type SNP &&"
    interval = get_label(regionfile)
    cmd = [ 'gatk', 
            'SelectVariants',
            '-L', interval,             
            '-R', genome,
            '-V', infile, 
            '-O', outfile, 
            '-select-type', 'SNP'
            ]
    try:
        run_command_shell(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise    
    

def gatk_vf(infile, outfile, genome, regionfile):
    '''
      maybe subprocess.run(" ".join(cmd), shell=True)
    
    '''
    " gatk VariantFiltration -R {params.gdir}/GRCh38.p7.genome.fa -L chr{params.chrom} "
    " -V {params.chrom}.snps.vcf -O {params.chrom}.snps_filtered.vcf "
    ' -filter "QD < 2.0" --filter-name "QD2" '
    ' -filter "QUAL < 30.0" --filter-name "QUAL30" '
    ' -filter "SOR > 3.0" --filter-name "SOR3" '
    ' -filter "FS > 60.0" --filter-name "FS60" '
    ' -filter "MQ < 40.0" --filter-name "MQ40" '
    ' -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" '
    ' -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" && '
    interval = get_label(regionfile)
    
    cmd = [ 'gatk', 
           'VariantFiltration',
            '-R', genome,
            '-L', interval, 
            '-V', infile, 
            '-O', outfile,            
            '-filter', '"QD < 2.0"', '--filter-name', '"QD2"', 
            '-filter', '"QUAL < 30.0"', '--filter-name', '"QUAL30"',
            '-filter', '"SOR > 3.0"', '--filter-name', '"SOR3"',
            '-filter', '"FS > 60.0"', '--filter-name', '"FS60"',
            '-filter', '"MQ < 40.0"', '--filter-name', '"MQ40"',
            '-filter', '"MQRankSum < -12.5"', '--filter-name', '"MQRankSum-12.5"',
            '-filter', '"ReadPosRankSum < -8.0"', '--filter-name',  '"ReadPosRankSum-8"'        
        ]
    try:
        run_command_shell(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise    
    
    
def igvtools_count(infile, outfile, reg_fa):
    """
    
    http://software.broadinstitute.org/software/igv/genome_ids
    
    igvtools count -z 0 -w 1 --bases --strands read {input.sfbam} 
    mv -v tmp.{wildcards.sample}.wig {output.splitwig}
    """
    cmd = [ 'igvtools',
            'count',
            '-z','0',
            '-w','1',
            '--bases',
            '--strands',
            'read',
            infile,
            outfile,
            reg_fa
            ]
    try:
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))         
        raise    
    
    
                  