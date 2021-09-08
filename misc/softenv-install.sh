#!/bin/bash -l
#  For werner1 project.   labelled gatk. 
#
#  samtools-1.11
#  bedtools-2.29.2
#  gatk-4.1.9.0
#  igvtools-2.8.9
#  STAR 2.7.2a  (tied to genome dir version)
#  genome 2.4.2a?
set -x

# Conda
conda create -y -n gatk python=3.8
conda activate gatk

# Conda available packages
conda install -y pandas numpy scipy seaborn matplotlib plotly h5py 
conda install -c conda-forge -c bioconda snakemake

# Tarball installs
#  Samtools 1.11
#  samtools 1.9 for compatiblity?
# ver=1.11
ver=1.9
wget https://github.com/samtools/samtools/releases/download/$ver/samtools-$ver.tar.bz2
tar -xvjf samtools-$ver.tar.bz2
cd samtools-$ver
./configure --prefix=$CONDA_PREFIX
make
make install
cd ..
rm -rf samtools-$ver.tar.bz2 samtools-$ver

#  bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
cp bin/* $CONDA_PREFIX/bin/
cd ..
rm -rf bedtools2 bedtools-2.29.1.tar.gz

#  gatk-4.1.9.0 
# ? 4.1.4.1 for compatibility?
# ver=4.1.9.0 
ver=4.1.4.1
wget https://github.com/broadinstitute/gatk/releases/download/$ver/gatk-$ver.zip
unzip gatk-$ver.zip
mv gatk-$ver $CONDA_PREFIX/
cd $CONDA_PREFIX/bin
ln -s ../gatk-$ver/gatk ./
cd -
rm gatk-$ver.zip

# igvtools-2.8.9
#
ver=2.8.10
pkg="IGV_${ver}"
wget https://data.broadinstitute.org/igv/projects/downloads/2.8/$pkg.zip
unzip $pkg.zip
mv $pkg $CONDA_PREFIX/
rm $pkg.zip

dpkg="IGV_Linux_${ver}_WithJava"
lpkg="IGV_Linux_${ver}"
wget https://data.broadinstitute.org/igv/projects/downloads/2.8/$dpkg.zip
unzip $dpkg.zip
rm $lpkg.zip

mv $lpkg  $CONDA_PREFIX/
mv $CONDA_PREFIX/$lpkg/jdk-11 $CONDA_PREFIX/$pkg/
cd $CONDA_PREFIX/bin/
ln -s $CONDA_PREFIX/$pkg/igvtools ./
cd -

# STAR
# may need to brew install gcc llvm libomp on mac. 
wget https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz
tar -xvzf STAR_2.4.2a.tar.gz
cd STAR-STAR-2.4.2a/source
#make STARforMacStatic CXX=/usr/local/Cellar/gcc/10.2.0/bin/g++-10 
make STAR
cp STAR $CONDA_PREFIX/bin/
cd -
rm -rf STAR-2.4.2a 2.4.2a.tar.gz

#
# Need newest Java in order to pass -XX:ActiveProcessorCount=8
# https://www.java.com/en/download/manual.jsp 
wget https://javadl.oracle.com/webapps/download/AutoDL?BundleId=243727_61ae65e088624f5aaa0b1d2d801acb16 -O jre-8u271-linux-x64.tar.gz
tar -xvzf jre-8u271-linux-x64.tar.gz
mv jre1.8.0_271 jre
cd $CONDA_PREFIX/bin/
ln -s cd $CONDA_PREFIX/jre/bin/* ./
cd -

	




