import datetime as dt
import logging
import shutil
import subprocess
import sys

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
    logging.info(f"command: {cmdstr} running...")
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
            cmdstring = ' '.join(cmd)
            logging.debug(f'command string is {cmdstring}')
            try:
                run_command(cmd)
            except NonZeroReturnException as nzre:
                logging.error(f'problem with {infile}')
                logging.error(traceback.format_exc(None))

def star_nowasp(end1, end2, outprefix, outtemp, nthreads, genomedir):
            
            cmd = ['STAR',
               '--readFilesIn', end1, end2, 
               '--outFileNamePrefix', outprefix,
               '--outTmpDir', outtemp, 
               '--runThreadN', nthreads,
               '--genomeDir', genomedir, 
               '--twopassMode Basic',
               '--twopass1readsN -1',
               '--outSAMtype BAM Unsorted', 
               '--quantMode GeneCounts'
               ]
            cmdstring = ' '.join(cmd)
            logging.debug(f'command string is {cmdstring}')
            try:
                run_command(cmd)
            except NonZeroReturnException as nzre:
                logging.error(f'problem with {infile}')
                logging.error(traceback.format_exc(None))
            finally:
                shutil.rmtree(f'{outprefix}._STARgenome')
                shutil.rmtree(f'{outprefix}._STARpass1')    
