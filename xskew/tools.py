import datetime as dt
import logging
import os
import pandas as pd
import shutil
import subprocess
import sys
import traceback

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


def run_command(cmd):
    """
    cmd should be standard list of tokens...  ['cmd','arg1','arg2'] with cmd on shell PATH.
    
    """
    cmdstr = " ".join(cmd)
    logging.info(f"running command: {cmdstr} ")
    start = dt.datetime.now()
    cp = subprocess.run(cmd, 
                    text=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT)
    end = dt.datetime.now()
    elapsed =  end - start
    logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {elapsed.seconds} seconds.")
    
    if cp.stderr is not None:
        logging.debug(f"got stderr: {cp.stderr}")
    if cp.stdout is not None:
        logging.debug(f"got stdout: {cp.stdout}")
    
    if str(cp.returncode) == '0':
        logging.info(f'successfully ran {cmdstr}')
        return(cp.stderr, cp.stdout,cp.returncode)
    else:
        logging.warn(f'non-zero return code for cmd {cmdstr}')


def run_command_shell(cmd):
    """
    maybe subprocess.run(" ".join(cmd), shell=True)
    cmd should be standard list of tokens...  ['cmd','arg1','arg2'] with cmd on shell PATH.
    
    """
    cmdstr = " ".join(cmd)
    logging.info(f"running command: {cmdstr} ")
    start = dt.datetime.now()
    cp = subprocess.run(" ".join(cmd), 
                    shell=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT)
    #cp = subprocess.run(cmd, 
    #                shell=True, 
    #                stdout=subprocess.PIPE, 
    #                stderr=subprocess.STDOUT)
    end = dt.datetime.now()
    elapsed =  end - start
    logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {elapsed.seconds} seconds.")
    
    if cp.stderr is not None:
        logging.debug(f"got stderr: {cp.stderr}")
    if cp.stdout is not None:
        logging.debug(f"got stdout: {cp.stdout}")
    
    if str(cp.returncode) == '0':
        logging.info(f'successfully ran {cmdstr}')
        return(cp.stderr, cp.stdout,cp.returncode)
    else:
        logging.error(f'non-zero return code for cmd {cmdstr}')
        raise NonZeroReturnException(f'For cmd {cmdstr}')

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


        
def fasterq_dump(infile, outdir, nthreads, tempdir ):
    
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
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {infile}')
        logging.error(traceback.format_exc(None))
        raise    
    
#def star_nowasp(end1, end2, outprefix, outtemp, nthreads, genomedir):

def star_genome(genomedir, nthreads, gtffile, infile ): 
    cmd = ['STAR',
           '--runMode', 'genomeGenerate',
           '--runThreadN', nthreads,
           '--genomeDir', genomedir,
           '--sjdbGTFfile', gtffile, 
           '--sjdbOverhang','100',
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
    colnames = ['Sequence-Name','Sequence-Role','Assigned-Molecule',
                'Assigned-Molecule-Location/Type','GenBank-Accn',
                'Relationship','RefSeq-Accn','Assembly-Unit',
                'Sequence-Length','UCSC-style-name']   
    df = pd.read_csv(reportfile, comment="#", sep='\t')
    df.columns = colnames
    #tagval = f'{chr}'
    #label = df[ df['UCSC-style-name'] == tagval]['RefSeq-Accn'].values[0]
    chrnum = chr[3:]
    tagval = chrnum
    label = df[ df['Assigned-Molecule'] == tagval]['RefSeq-Accn'].values[0]
    logging.debug(f'extracted label {label} for {tagval} in {reportfile}')
    f = open(outfile, 'w')
    f.write(f'{label}\n')
    f.close()

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


def make_chr_index(infile, genomedir, chr, outfile):
    region = get_chr_label(genomedir, chr)
    samtools_faidx_region(infile, outfile, region)  


def star_nowasp(end1, end2, outprefix, nthreads, genomedir):
            
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
        run_command(cmd)
    except NonZeroReturnException as nzre:
        logging.error(f'problem with {end1}/{end2} input files.')
        logging.error(traceback.format_exc(None))
    finally:
        for ext in STARSUBDIRS:
            dirpath = f'{outprefix}{ext}'
            if os.path.exists(dirpath) and os.path.isdir(dirpath):
                shutil.rmtree(dirpath)    

#def star_wasp(end1, end2, vcf, outprefix, outtemp, nthreads, genomedir):
def star_wasp(end1, end2, vcf, outprefix, nthreads, genomedir):
    " STAR  --genomeDir {params.gdir} --readFilesIn {input.end1} {input.end2} "
    " --runThreadN {threads} --twopassMode Basic --twopass1readsN -1 " 
    " --outSAMtype BAM Unsorted --quantMode GeneCounts "
    " --waspOutputMode SAMtag --varVCFfile {input.sfvcf} && "            
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
        raise    
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
    
    
                  