# 
conda create -n xskew
conda activate xskew
#conda install -c conda-forge -c bioconda snakemake
#conda install -c bioconda star gatk samtools bamtools sra-tools=2.11.0


conda install -c conda-forge -c bioconda python snakemake=6.8.0 star gatk4 samtools bamtools sra-tools=2.11.0 igvtools openjdk

# Move jdk11 to igvtoolsdir. 
cd ~/src 
wget https://data.broadinstitute.org/igv/projects/downloads/2.5/IGV_Linux_2.5.3.zip
unzip IGV_Linux_2.5.3.zip 
mv IGV_Linux_2.5.3/jdk-11 $CONDA_PREFIX/share/igvtools-2.5.3-1/ 
