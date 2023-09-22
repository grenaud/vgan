[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/vgan/README.html)
# vgan

vgan is a suite of tools for pangenomics. We currently support two main subcommands: Haplocart (for modern human mtDNA haplogroup classification) and euka (for bilaterian abundance estimation of ancient environmental DNA). The underlying data structure is the VG graph (see https://github.com/vgteam/vg).

## Installation:

vgan is supported for use on Linux systems.


### Release build

The easiest way to run vgan is to download the static binary. 

Step 1: Download the static binary. The list of releases is here: https://github.com/grenaud/vgan/tags Each release comes with a static binary. Find the URL of the static binary by right-clicking and selecting "copy link"

```
wget [URL TO RELEASE BINARY]
```

Where you paste the URL of the binary.  If you have root access, simply install the executable by running

```
sudo cp vgan /usr/bin
```

Otherwise, just leave the executable where it is or copy it in a bin/ directory in your home directory:

```
mkdir -P $HOME/bin/
cp vgan $HOME/bin/
```

Step 2: Mark the binary executable:

```
chmod +x vgan
```

Step 3: Download the required graph files. If you have root access please run

For HaploCart:
```
sudo mkdir -p /usr/bin/share/vgan/hcfiles/
sudo wget -nc -l1 --recursive --no-parent -P /usr/bin/share/vgan/ ftp://ftp.healthtech.dtu.dk:/public/haplocart/hcfiles/
```

For euka:
```
sudo mkdir -p /usr/bin/share/vgan/euka_dir/
sudo wget -nc -l1 --recursive --no-parent -P /usr/bin/share/vgan/euka_dir/ ftp://ftp.healthtech.dtu.dk:/public/euka_files/
sudo mkdir -p /usr/bin/share/vgan/damageProfiles/
sudo wget -O $/usr/bin/share/vgan/damageProfiles/none.prof https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/damageProfiles/none.prof
```
If euka is to be used by multiple users, please ensure that the file ```/usr/bin/share/vgan/euka_dir/euka_db.dist``` has writing permission with the following command:
```
sudo chmod +w /usr/bin/share/vgan/euka_dir/euka_db.dist
```
If you do not have root access, you can download them in the directory of your choice but for ease, you can download them to a share directory in your home folder:

For HaploCart:
```
mkdir -p $HOME/share/vgan/
mkdir -p $HOME/share/vgan/hcfiles/
wget -nc -l1 --recursive --no-directories --no-parent -P $HOME/share/vgan/hcfiles/ ftp://ftp.healthtech.dtu.dk:/public/haplocart/hcfiles/
```
For euka:
```
mkdir -p $HOME/share/vgan/euka_dir/
wget -nc -l1 --recursive --no-directories --no-parent -P $HOME/share/vgan/euka_dir/ ftp://ftp.healthtech.dtu.dk:/public/euka_files/
mkdir -p $HOME/share/vgan/damageProfiles/
wget -O $HOME/share/vgan/damageProfiles/none.prof https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/damageProfiles/none.prof
```
If euka is to be used by multiple users, please ensure that the file ```$HOME/share/vgan/euka_dir/euka_db.dist``` has writing permission with the following command:
```
sudo chmod +w  $HOME/share/vgan/euka_dir/euka_db.dist
```

Additionally, you can download eukas visualisation scripts:
```
wget  -O $HOME/bin/make_tree_from_output.py https://raw.githubusercontent.com/grenaud/vgan/main/bin/make_tree_from_output.py
wget  -O $HOME/bin/plot_taxon.R https://raw.githubusercontent.com/grenaud/vgan/main/bin/plot_taxon.R
wget  -O $HOME/bin/visualize_detected_taxa.sh https://raw.githubusercontent.com/grenaud/vgan/main/bin/visualize_detected_taxa.sh
wget  -O $HOME/bin/euka.yml https://raw.githubusercontent.com/grenaud/vgan/main/bin/euka.yml
```

Alternatively, you can install the graph files elsewhere and specify them using:

For HaploCart:
```
vgan haplocart --hc-files 
```
For euka:
```
vgan euka --euka_dir
```

After this please add the ```$HOME/bin/``` directory to your path.  

```
export PATH="$HOME/bin:$PATH"
```
or add this line at the bottom of your $HOME/.bashrc to permanently add ```$HOME/bin/```  to your path.


### Bioconda

Install Conda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and Bioconda (https://bioconda.github.io/) and type the following:

```
 conda  install -c bioconda vgan
```

To use euka from bioconda, it is necessary to download eukas graph files externally with the following command:
```
mkdir ./euka_dir
wget -nc -l1 --recursive --no-directories --no-parent -P euka_dir/ ftp://ftp.healthtech.dtu.dk:/public/euka_files/
```
When calling the euka command, please make sure to add ```--euka_dir``` option with the correct input path to the euka_dir directory with the graph file.
eukas associated graph files have a size of 11Gb; please make sure that enough space is available when downloading. 

If euka is to be used by multiple users, please ensure that the file ```euka_dir/euka_db.dist``` has writing permission with the following command:
```
sudo chmod +w euka_dir/euka_db.dist
```

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
git clone https://github.com/grenaud/vgan.git
cd vgan/src && make
```

Please be patient, this will take some time.

If you are building on a WSL subsystem you may receive an error due to VG not building entirely.
This is likely to be a permissions issue. Please run the above command as a root user. (i.e. cd vgan/src && sudo make)

If the build has been successful, the executable will be found in the bin folder.

If you build from source, you can download the requisite files for HaploCart with 

```
make hcfilesmade
```

and for euka:

```
make eukafilesmade
```
If euka is to be used by multiple users, please ensure that the file ```$HOME/share/vgan/euka_dir/euka_db.dist``` has writing permission with the following command:
```
sudo chmod +w  $HOME/share/vgan/euka_dir/euka_db.dist
```

from the src directory

### Docker container

First pull the latest image, for example if the latest image is 1.0.2 you would do:
```
docker pull gabrielreno/vganv1.0.2:latest
```
Then run vgan as such:
```
docker run gabrielreno/vganv1.0.2:latest ../bin/vgan
```

### Potential Issues

vgan requires VG as a dependency. If there are issues with making vg please see https://github.com/vgteam/vg for further documentation. vg requires several packages to install. 

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

## euka 

euka is a tool to characterize ancient environmental DNA samples for arthropodic and tetrapodic mitochondrial DNA. The foundation of this tool is a curated database of complete mitochondrial genomes sorted into taxonomic groups and built into a variation graph reference structure. euka takes FASTQ input files and maps the file against our reference graph structure.  We use an ancient-aware maximum-likelihood framework to analyse our mapping results on a per fragment basis. Additionally, we filter on mapping quality and coverage evenness across the pangenome graph of a taxa. 

### Usage

For merged or single-end FASTQ (gzipped okay):
```
./vgan euka -fq1 [FASTQ file]
```
Please note that if you would like to input multiple FASTQ file, we encourage the use of file descriptors:
```
./vgan euka -fq1 <(zcat [FASTQ file1] [FASTQ file2] [FASTQ file3] ...)
```
For paired-end FASTQ (gzipped okay):
```
./vgan euka -fq1 [FASTQ file] -fq2 [FASTQ file]
```

For incorporation of ancient damage profiles during the maximum likelihood estimation: 
```
./vgan euka -fq1 [FASTQ file] --deam5p [5end.prof] --deam3p [3end.prof]
```

For user specific MCMC runs:
```
./vgan euka -fq1 [FASTQ file] -iter [INT] -burnin [INT]
```

For additional output options:
```
./vgan euka -fq1 [FASTQ file] --outFrag --outGroup [Taxa name]
```
Please note that the Taxa name provided to the ```--outGroup``` argument must be a taxon defined in euka's database. 

If you specify the ```--outFrag``` option you will be provided with a list of all read names that belong to the detected taxa (one taxon per row). To extract these reads for further downstream analysis you can use the following command:
```
less -S [output file prefix]_FragNames.tsv | sed '[row number]!d' | sed 's/\t/\n/g' | seqtk subseq [FASTQ input file] - > output.fq
```
### euka filter options:
euka has four filter options to modify the stringency of taxa detection. We always recommend to have a first look at your samples with the default parameters. These parameters have been thoroughly tested to provide confident abundance estimations. However, to detect more divergent taxa (for example Formicidae), it may be necessary to adjust the filter parameters. Furthermore, we want to highlight that the reference genomes for many arthropodic species have low-complexity reference genomes and are more prone to spurious alignments. Even with our standard parameter filters, we can see more false-positive detections for these taxa, and results should always be evaluated carefully. Our ~/vgan/share/vgan/euka_dir/euka_db.bins file lists every taxa in our database with their respective bins (for our coverage estimation). The file shows the Node ID range for each bin and the calculated entropy score for this bin.


Example for incorporation of lower-entropy regions in the mitogenome: 
```
./vgan euka -fq1 [FASTQ file] --entropy 1.13 --minBins 3 
```

Example for more conservative parameter settings: 
```
./vgan euka -fq1 [FASTQ file] --minMQ 40 --minFrag 20
```


### euka Options: 
A list of all of eukas options:

<pre>
  Input options:
                --euka_dir [STR]        euka database location (default: current working directory)
                --dbprefix [STR]        database prefix name (defualt: euka_db)
                -fq1 [STR]              Input FASTQ file (for merged and single-end reads)
                -fq2 [STR]              Second input FASTQ file (for paired-end reads)
                -i                      Paired-end reads are interleaved (default: false)
                -o [STR]                Output file prefix (default: euka_output)
                -t                      Number of threads (-1 for all available)
                -Z                      Temporary directory (default: /tmp)

Filter options:
                --minMQ [INT]           Set the mapping quality minimum for a fragment (default: 29)
                --minFrag [INT]         Minimum amount of fragments that need to map to a group (default: 10)
                --entropy [double]      Minimum entropy score for a bin to be considered (default: 1.17)
                --minBins [INT]         Minimum number of bins that need to be available for a group (default: 6)

Damage options:
                --deam5p                [.prof] 5p deamination frequency for eukaryotic species (default: no damage)
                --deam3p                [.prof] 3p deamination frequency for eukaryotic species (default: no damage)
                -l [INT]                Set length for substitution matrix (default: 5)
                --out_dir [STR]         Path for output prof-file (default: current working directory)

Markov chain Monte Carlo options:
                --no-mcmc               The MCMC does not run (default: false)
                --iter [INT]            Define the number of iterartions for the MCMC (default: 10000)
                --burnin [INT]          Define the burnin period for the MCMC (default: 100)

Additional output option:
                --outFrag               Outputs a file with all read names per taxonomic group (default: false)
                --outGroup [string]     Outputs all information about a taxonmic group of interest (default: empty)
</pre>  


### euka output files:
euka's default output consists of four TSV files:

The abundance and the detected TSV file consists of the following columns:
``` 
| Taxa | detected | Number of reads | proportion estimate | 85% CI lower bound | 85% CI upper bound | 95% CI lower bound | 95% CI upper bound |
```
The abundance file will list all availabe taxa, while the detected file will only list taxa that have passed all of euka's detection filters.

There are multiple ways to visualise eukas output. To make sure the provided scripts work, we recommend installing the provided conda environment:
```
conda env create -f ~/vgan/bin/euka.yml
```
Or make sure the following packages are installed for python ete3, csv and argparse; for R, you will need the libraries ggplot2 and ggpubr. 

You have the option to visualize these two files with a taxonomic tree using the provided Python script (make sure Python is in your path):
``` 
python ~/vgan/bin/make_tree_from_ouput.py [output file prefix]_abundance.tsv
python ~/vgan/bin/make_tree_from_ouput.py [output file prefix]_detected.tsv 
```

The coverage file consists of all detected taxa and their corresponding number of bins. It shows the number of reads sorted into each bin. The coverage TSV file is used by our plotting script to show for a barplot.

The inSize file provides a list of all fragment sizes for each taxa. It is used by our plotting script to create a histogram of the fragment length distribution.

For each detected taxa we estimate a damage profile in a .prof file. These files are used as input for the damage estimation plots in our plotting script. 
Additionally, we provide an average damage profile for the 5’ end and the 3’ end that can be directly inputted into our damage model options ``` --deam5p ``` and ``` --deam3p ```.

To visualize detected taxa you can use the following provided script:
``` 
~/vgan/bin/visualize_detected_taxa.sh [output file prefix] 
```
To visualize only a specific detected taxon please use: 
```
Rscript plot_taxon.R [output file prefix] [taxon name]
```

### euka notes:
- We recommend to run euka first with standard parameters and no incorporated damage profiles.
- Further, we encourage the use of all available input files for one ecological site. Especially with ancient environmental DNA samples where read counts/abundances are extremely low euka will have difficulties “detecting” a taxa if, for example, only one sequencing lane is provided. 

## Unit tests:

We provide comprehensive unit tests for our main subcommands (HaploCart and Euka). These are built on the Boost Test Suite. Running the tests requires building vgan from source because it involves compiling a separate executable. To run the unit tests, please run

```
make test
./../bin/test --run_test
```

You can also run subcommand-specific tests:

```
make test
./../bin/test --run_test=haplocart
```

or

```
make test
./../bin/test --run_test=euka
```

## General notes:
- Please be aware that the multithreading for paired-end input files does not work with the vg version we use. However, from our experience, vg giraffe often reverts to mapping single-end, which is multi-threaded. If you have paired-end reads please consider merging them or concatenate them with a file descriptor:
```
vgan euka -fq1 <(zcat reads1.fq.gz reads2.fq.gz) -t 30 

```

## Support:

Vgan is actively supported and maintained.
If you discover an issue with the program, please submit a bug report on Github or send an email. 

Contact:<br>

  Haplocart: jdru@dtu.dk or gabriel.reno@gmail.com<br>
  euka: navo@dtu.dk or gabriel.reno@gmail.com
