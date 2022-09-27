# vgan

vgan is a suite of tools for mitochondrial pangenomics. We currently support two subcommands: Haplocart (for modern human mtDNA haplogroup classification)
and duprm (for PCR duplicate removal from GAM files). The underlying data structure is the VG graph (see https://github.com/vgteam/vg).

## Installation:

vgan is supported for use on Linux systems.


### Release build

The easiest way to run vgan is to download the static binary. 

Step 1: Download the static binary

```
wget -nc -P $HOME/bin/ [URL WHERE WE HOST THE EXECUTABLE]
```

If you have root access, simply install the executable by running

```
sudo cp vgan /usr/bin
```

Otherwise, just leave the executable where it is or copy it in a bin/ directory in your home directory:

```
mkdir -P $HOME/bin/
sudo cp vgan $HOME/bin/
```

Step 2: Download HaploCart graph files. If you have root access please run

```
sudo mkdir -p /usr/bin/share/hcfiles/
sudo wget -nc -l1 --recursive --no-parent -P /usr/bin/share/ ftp://ftp.healthtech.dtu.dk:/public/haplocart/hcfiles/
```

If you do not have root access, you can download them in the directory of your choice but for ease, you can download them to a share directory in your home folder:

```
mkdir -P $HOME/share/
mkdir -P $HOME/share/hcfiles/
wget -nc -l1 --recursive --no-directories --no-parent -P $HOME/share/hcfiles/ ftp://ftp.healthtech.dtu.dk:/public/haplocart/hcfiles/
```

Alternatively, you can install the HaploCart graph files elsewhere and specify them using:
```
vgan haplocart --hc-files 
```

After this please add the ```$HOME/bin/``` directory to your path.  

```
export PATH="$HOME/bin:$PATH"
```
or add this line at the bottom of your $HOME/.bashrc to permanently add ```$HOME/bin/```  to your path.


### Bioconda

Install conda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and bioconda (https://bioconda.github.io/) and type the following:

```
 conda  install -c bioconda  vgan
````

### Building vgan from source

Please ensure that git is installed on your system.
You can check by typing "git --version".

If you wish to build from source you will need the build-essential package. If it is not already installed, please run

```
sudo apt-get install build-essential
```

You will also need wget. If you do not have it please run

```
sudo apt-get  install wget
```

vgan requires vg https://github.com/vgteam/vg. vg requires several packages to be installed, please refer to the README in vg. 

```
git clone --depth 1 https://github.com/grenaud/vgan.git
cd vgan/src && make
```

Please be patient, this will take some time.

If you are building on a WSL subsystem you may receive an error due to VG not building entirely.
This is likely to be a permissions issue. Please run the above command as a root user. (i.e. cd vgan/src && sudo make)

If the build has been successful, the executable will be found in the bin folder.


### Docker container

First pull the latest image, for example if the latest image is 1.0.0 you would do:
```
docker pull gabrielreno/vganv1.0.0:latest
```
Then run vgan as such:
```
docker run gabrielreno/vganv1.0.0:latest ../bin/vgan
```

### Potential Issues

vgan requires VG as a dependency. If there are issues with making vg please see https://github.com/vgteam/vg for further documentation. vg requires several packages to install. 

If you are on a WSL subsystem you may get the following error message:

```
make: *** [Makefile:65: vgmade] Error 2
```

In this case please make vg manually:

```
cd ../dep/vg && ./source_me.sh && make 
``` 
And then try again making vgan.


If your system is complaining that it cannot find package XYZ in the pkg-config search path, you 
need to modify your PKG_CONFIG_PATH environment variable. First find the location of your pkg-config files by running 

```
find / -iname 'XYZ.pc' 2>/dev/null
```

Then, set 

```
export PKG_CONFIG_PATH=[PATH TO PKG-CONFIG FILES]
```

Also, you may see the error "cc1plus: error: bad value (‘tigerlake’) for ‘-march=’ switch"

If this happens please run:

```
sudo apt install gcc-10 g++-10
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100 --slave /usr/bin/g++ g++ /usr/bin/g++-10 --slave /usr/bin/gcov gcov /usr/bin/gcov-10
sudo update-alternatives --config gcc
make clean && make
```
## HaploCart

HaploCart performs maximum likelihood estimation to predict the mitochondrial haplogroup for reads 
originating from an uncontaminated modern human sample. The program also optionally provides confidence estimation of the phylogenetic placement of the sample.

### Quick start

For consensus FASTA:

```
wget ftp://ftp.healthtech.dtu.dk:/public/haplocart/testdata/rCRS.fa.gz
vgan haplocart -f rCRS.fa.gz
```

For single-end FASTQ:

```
wget ftp://ftp.healthtech.dtu.dk:/public/haplocart/testdata/rCRS.fq.gz
vgan haplocart -fq1 rCRS.fq.gz
```


### Usage

For GAM input (note that the reads must have previously been aligned to OUR graph):

```
vgan haplocart -g [GAM file]
```

For consensus FASTA (may be gzipped):

```
vgan haplocart -f [FASTA file]
```

For single-end FASTQ (may be gzipped):

```
vgan haplocart -fq1 [FASTQ file]
```

For single-end interleaved (may be gzipped):

```
vgan haplocart -fq1 [FASTQ file] -i
```

For paired-end FASTQ (may be gzipped):

```
vgan haplocart -fq1 [FASTQ file] -fq2 [FASTQ file]
```

BAM/CRAM input:

Haplocart does not yet accept HTS-lib-compatible input directly. Please first unmap the reads and then follow the instructions for FASTQ input described above, e.g.

```
samtools bam2fq [FASTQ file] | vgan haplocart -fq1 /dev/stdin
```

### Haplocart Output File

The output file for Haplocart is a TSV file with the following columns:

```
| Sample | Predicted haplogroup	| # Mapped reads |
```

### Haplocart Posterior File

The output file for HaploCart's optional posterior probability estimation is a TSV file with the following columns:

```
| Sample | Considered clade | Posterior probability | Tree depth |
```
  
### HaploCart Notes

- HaploCart will make a prediction on any sample for which at least one read maps. It is up to the user to apply appropriate filters to HaploCart results.
If your input data is especially sparse please check the clade-level posterior probabilities and treat the predictions judiciously. 

- Although the HaploCart inference algorithm is deterministic, we rely upon VG Giraffe for mapping and we have noticed a degree of stochasticity in the mapping procedure. 
  Therefore in some edge cases there may be a (very slight) difference in output with different runs on the same sample. Unfortunately Giraffe has no random seed, so we cannot provide one either.
 

## duprm

vgan duprm removes PCR duplicate reads from a (SORTED!) GAM file. Usage is simply

```
vgan duprm [sorted_in.gam] > [out.gam]
```

## Support:

vgan is actively supported and maintained.
If you discover an issue with the program, please submit a bug report on Github or send an email. 

Contact:

  Haplocart: jdru@dtu.dk or gabriel.reno@gmail.com



