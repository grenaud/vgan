[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?styl=flat)](http://bioconda.github.io/recipes/vgan/README.html)
# vgan

vgan is a suite of tools for pangenomics. We currently support four main subcommands: Haplocart (for modern human mtDNA haplogroup classification); euka (for bilaterian abundance estimation of ancient environmental DNA);
soibean (for sedaDNA mitochondrial mixtures of eukaryotes); and TrailMix (for hominin ancient mtDNA samples, possibly mixtures).
The underlying data structure is the VG graph (see https://github.com/vgteam/vg).

# Installation:

vgan is supported for use on Linux systems.


## Release build

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
For soibean:
```
sudo mkdir -p /usr/bin/share/vgan/soibean_dir/
sudo wget -nc --recursive --no-parent -P /usr/bin/share/vgan/soibean_dir/ ftp://ftp.healthtech.dtu.dk:/public/soibean_files/
sudo mkdir -p /usr/bin/share/vgan/damageProfiles/
sudo wget -O $/usr/bin/share/vgan/damageProfiles/none.prof https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/damageProfiles/none.prof
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
wget -nc -l1 --recursive --no-parent -P $HOME/share/vgan/euka_dir/ ftp://ftp.healthtech.dtu.dk:/public/euka_files/
mkdir -p $HOME/share/vgan/damageProfiles/
wget -O $HOME/share/vgan/damageProfiles/none.prof https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/damageProfiles/none.prof
```
If euka is to be used by multiple users, please ensure that the file ```$HOME/share/vgan/euka_dir/euka_db.dist``` has writing permission with the following command:
```
sudo chmod +w  $HOME/share/vgan/euka_dir/euka_db.dist
```
Additionally, you can download eukas visualisation scripts:
```
wget -O $HOME/share/vgan/plottingScripts/make_tree_from_output.py https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/make_tree_from_output.py
wget -O $HOME/share/vgan/plottingScripts/plot_taxon.R https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/plot_taxon.R
wget -O $HOME/share/vgan/plottingScripts/visualize_detected_taxa.sh https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/visualize_detected_taxa.sh
wget -O $HOME/share/vgan/plottingScripts/euka.yml https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/euka.yml
```
For soibean:
```
mkdir -p $HOME/share/vgan/soibean_dir/
wget -nc --recursive --no-parent -P $HOME/share/vgan/soibean_dir/ ftp://ftp.healthtech.dtu.dk:/public/soibean_files/
mkdir -p $HOME/share/vgan/damageProfiles/
wget -O $HOME/share/vgan/damageProfiles/none.prof https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/damageProfiles/none.prof
```
The tree files need to be unzipped:
```
cd $HOME/share/vgan/soibean_dir/tree_dir/
unzip tree.zip
```
To download all additional scripts for soibean, please use the following commands:
```
wegt -O $HOME/share/vgan/soibean_dir/make_graph_files.sh https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/soibean_dir/make_graph_files.sh
chmod +x $HOME/share/vgan/soibean_dir/make_graph_files.sh
wget -O $HOME/share/vgan/plottingScripts/soibeanPlotTrace.R https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/soibeanPlotTrace.R
wget -O $HOME/share/vgan/plottingScripts/soibeanPlots.R https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/soibeanPlots.R
wget -O $HOME/share/vgan/plottingScripts/soibeanPlotk.R https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/soibeanPlotk.R
wget -O $HOME/share/vgan/plottingScripts/soibean.yml https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/soibean.yml
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
For soibean:
```
vgan soibean --soibean_dir
```

After this please add the ```$HOME/bin/``` directory to your path.  

```
export PATH="$HOME/bin:$PATH"
```
or add this line at the bottom of your $HOME/.bashrc to permanently add ```$HOME/bin/```  to your path.


## Bioconda

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

## Building vgan from source

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

and for soibean:

```
make soibeanfilesmade
```

## Docker container

First pull the latest image, for example if the latest image is 1.0.2 you would do:
```
docker pull gabrielreno/vganv1.0.2:latest
```
Then run vgan as such:
```
docker run gabrielreno/vganv1.0.2:latest ../bin/vgan
```

## Potential Issues

vgan requires VG as a dependency. If there are issues with making vg please see https://github.com/vgteam/vg for further documentation. vg requires several packages to install. 

# HaploCart

HaploCart performs maximum likelihood estimation to predict the mitochondrial haplogroup for reads 
originating from an uncontaminated modern human sample. The program also optionally provides confidence estimation of the phylogenetic placement of the sample.

## Quick start

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


## Usage

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

## Haplocart Output File

The output file for Haplocart is a TSV file with the following columns:

```
| Sample | Predicted haplogroup	| # Mapped reads |
```

## Haplocart Posterior File

The output file for HaploCart's optional posterior probability estimation is a TSV file with the following columns:

```
| Sample | Considered clade | Posterior probability | Tree depth |
```
  
## HaploCart Notes

- HaploCart will make a prediction on any sample for which at least one read maps. It is up to the user to apply appropriate filters to HaploCart results.
If your input data is especially sparse please check the clade-level posterior probabilities and treat the predictions judiciously. 

- Although the HaploCart inference algorithm is deterministic, we rely upon VG Giraffe for mapping and we have noticed a degree of stochasticity in the mapping procedure. 
  Therefore in some edge cases there may be a (very slight) difference in output with different runs on the same sample. Unfortunately Giraffe has no random seed, so we cannot provide one either.

- Please remember, if the graph file and its indexes are used by multiple users the graph.dist file needs to be writeable for all. You can use the following command:
```
chmod +w [graph.dist]
```
 

# duprm

vgan duprm removes PCR duplicate reads from a (SORTED!) GAM file. Usage is simply

```
vgan duprm [sorted_in.gam] > [out.gam]
```

# euka 

euka is a tool for characterizing ancient environmental DNA samples for arthropodic and tetrapodic mitochondrial DNA. Its foundation is a curated database of complete mitochondrial genomes sorted into taxonomic groups and built into a variation graph reference structure. euka takes FASTQ input files and maps them against our reference graph structure.  We use an ancient-aware maximum-likelihood framework to analyse our mapping results on a per-fragment basis. Additionally, we filter on mapping quality and coverage evenness across the pangenome graph of a taxa. 

## Quick start
For a metagenomic fastq file:
```
wget https://github.com/grenaud/vgan/raw/main/test/input_files/euka/three_dhigh_100.fq.gz
vgan euka -fq1 three_dhigh_100.fq.gz -t 20
```
For a full workflow example, please see our tab: Example workflow.

## Usage

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

For user-specific MCMC runs:
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
## euka filter options:
euka has four filter options to modify the stringency of taxa detection. We always recommend to have a first look at your samples with the default parameters. These parameters have been thoroughly tested to provide confident abundance estimations. However, to detect more divergent taxa (for example Formicidae), it may be necessary to adjust the filter parameters. Furthermore, we want to highlight that the reference genomes for many arthropodic species have low-complexity reference genomes and are more prone to spurious alignments. Even with our standard parameter filters, we can see more false-positive detections for these taxa, and results should always be evaluated carefully. Our ~/vgan/share/vgan/euka_dir/euka_db.bins file lists every taxa in our database with their respective bins (for our coverage estimation). The file shows the Node ID range for each bin and the calculated entropy score for this bin.


Example for incorporation of lower-entropy regions in the mitogenome: 
```
./vgan euka -fq1 [FASTQ file] --entropy 1.13 --minBins 3 
```

Example for more conservative parameter settings: 
```
./vgan euka -fq1 [FASTQ file] --minMQ 40 --minFrag 20
```


## euka Options: 
A list of all of euka's options:

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


## euka output files:
euka's default output consists of four TSV files:

The abundance and the detected TSV file consists of the following columns:
``` 
| Taxa | detected | Number of reads | proportion estimate | 85% CI lower bound | 85% CI upper bound | 95% CI lower bound | 95% CI upper bound |
```
The abundance file will list all availabe taxa, while the detected file will only list taxa that have passed all of euka's detection filters.

There are multiple ways to visualise eukas output. To make sure the provided scripts work, we recommend installing the provided conda environment:
```
conda env create -f ~/vgan/share/vgan/plottingScripts/euka.yml
```
Or make sure the following packages are installed for python ete3, csv and argparse; for R, you will need the libraries ggplot2 and ggpubr. 

You have the option to visualize these two files with a taxonomic tree using the provided Python script (make sure Python is in your path):
``` 
python ~/vgan/share/vgan/plottingScripts/make_tree_from_ouput.py [output file prefix]_abundance.tsv
python ~/vgan/share/vgan/plottingScripts/make_tree_from_ouput.py [output file prefix]_detected.tsv 
```

The coverage file consists of all detected taxa and their corresponding number of bins. It shows the number of reads sorted into each bin. The coverage TSV file is used by our plotting script to show for a barplot.

The inSize file provides a list of all fragment sizes for each taxa. It is used by our plotting script to create a histogram of the fragment length distribution.

For each detected taxa, we estimate a damage profile in a .prof file. These files are used as input for the damage estimation plots in our plotting script. 
Additionally, we provide an average damage profile for the 5’ end and the 3’ end that can be directly inputted into our damage model options ``` --deam5p ``` and ``` --deam3p ```.

To visualize detected taxa, you can use the following provided script:
``` 
~/vgan/share/vgan/plottingScripts/visualize_detected_taxa.sh [output file prefix] 
```
To visualize only a specific detected taxon, please use: 
```
Rscript plot_taxon.R [output file prefix] [taxon name]
```

## euka notes:
- We recommend running euka first with standard parameters and no incorporated damage profiles.
- Further, we encourage the use of all available input files for one ecological site. Especially with ancient environmental DNA samples where read counts/abundances are extremely low, euka will have difficulties “detecting” a taxon if, for example, only one sequencing lane is provided.
- Please remember, if the graph file and its indexes are used by multiple users, the euka graph.dist file needs to be writeable for all. You can use the following command:
```
chmod +w euka_db.dist
```

# soibean
soibean is a tool crafted for identifying and placing ancient environmental DNA from **pre-analysed FASTQ files** onto a phylogenetic tree. Pre-analysed FASTQ files are all FASTQ reads that have been **previously classified to a lower taxonomic resolution** (e.g. family level). It leverages mitochondrial pangenome graphs across 335 arthropodic and tetrapodic taxa, drawing on the euka database enhanced with ancestral states and complemented by HKY model phylogenetic trees. Soibean maps input FASTQs to a specific subgraph, determining the likely species or species mix within. It uses an MCMC algorithm to estimate the species' positions on the tree and their proportions in the sample.

## Quickstart
For fastq mapped to the taxon of bears (Ursidae):
```
wget https://github.com/grenaud/vgan/raw/main/test/input_files/soibean/k1.fq.gz
cd $HOME/vgan/share/vgan/soibean_dir/tree_dir
unzip trees.zip
./make_graph_files.sh taxa
./make_graph_files.sh Ursidae
vgan soibean -fq1 k1.fq.gz --dbprefix Ursidae -t 20
```

## Example workflow:
A possible workflow for your ancient environmental DNA sample could be:
```
wget https://github.com/grenaud/vgan/raw/main/test/input_files/euka/test.fq.gz

vgan euka -fq1 test.fq.gz --outFrag 
```
The taxon of bears (Ursidae) has been detected. We plot the taxon to verify that it looks ancient:
```
./visualize_detected_taxa.sh euka_output
```

After the verification, we extract the reads mapped to the bear taxon to do species detection:
```
less -S euka_output_FragNames.tsv | sed '3!d' | sed 's/\t/\n/g' | seqtk subseq test.fq.gz - | gzip > Ursidae.fq.gz
```
We extract the subgraph for the bears from our soibean graph:
``` 
./make_graph_files.sh Ursidae 
```
The script extracts the bear subgraph and creates all necessary files for mapping and analysing with soibean:
```
vgan soibean -fq1 Ursidae.fq.gz --dbprefix Ursidae -t 20 -o soibeanOut 
```

Once soibean is done we can plot our results. First, we check the k-curve, which shows the number of sources best describing our data:
```
Rscript soibeanPlotk.R soibeanOut
```
Afterwards, we can choose the best run of our MCMC and plot the distribution on the phylogenetic tree (the tree files can be found: ```~/vgan/share/vgan/soibean_dir/tree_dir```.:
```
Rscript soibeanPlots.R Ursidae.new.dnd soibeanOutResult10.mcmc 
```

### Usage
Extract a subgraph of interest:
```
./make_graph_files.sh [taxon name] [optional: path to vg binary]
```
For example, to extract the graph for the taxon of bears (Ursidae):
```
./make_graph_files.sh Ursidae
```
You will get all the necessary graph files for mapping with soibean as output. You only have to create these files once per taxon.

For merged or single-end FASTQ (gzipped okay):
```
./vgan soibean -fq1 [FASTQ file] --dbprefix [taxon name]
```
Please note that if you would like to input multiple FASTQ files, we encourage the use of file descriptors:
```
./vgan soibean -fq1 <(cat [FASTQ file1] [FASTQ file2] [FASTQ file3] ...) --dbprefix [taxon name]
```
For paired-end FASTQ (gzipped okay):
```
./vgan soibean -fq1 [FASTQ file] -fq2 [FASTQ file] --dbprefix [taxon name]
```
For incorporation of ancient damage profiles during the maximum likelihood estimation: 
```
./vgan soibean -fq1 [FASTQ file] --dbprefix [taxon name] --deam5p [5end.prof] --deam3p [3end.prof]
```
An example of a damage file can be found: ```~/vgan/share/vgan/damageProfiles/dhigh5p.prof``` and ```~/vgan/share/vgan/damageProfiles/dhigh3p.prof```. 
euka provides damage profiles for each detected taxon ending in ```.prof```. The taxon-specific damage profiles can be adjusted to be used as soibean input:
```
line=$(grep -n "^A>C" [euka taxon output].prof | head -2 | tail -1 | cut -d: -f1) && awk -v line=$line -v OFS='\t' 'BEGIN{FS=OFS} NR < line {NF--; print}' [euka taxon output].prof > [taxon 5end].prof && (head -1 [euka taxon output].prof && awk -v line=$line -v OFS='\t' 'BEGIN{FS=OFS} NR >= line {if(NR > line) {NF--; print}}' [euka taxon output].prof | tac) > [taxon 3end].prof
```

For user-defined options in the MCMC, e.g., the number of iterations, the number of burn-in iterations, or the maximum number of sources (k):
```
./vgan soibean -fq1 [FASTQ file] --dbprefix [taxon name] --iter [INT] --burnin [INT] -k [INT]
```
or

```
./vgan soibean -fq1 [FASTQ file] --dbprefix [taxon name] --iter [INT] --burnin [INT] --randStart true
```

## soibean Options: 
A list of all of soibeans options:

<pre>
Options:
  --soibean_dir <directory>   Specify the directory containing the soibean files
  --tree_dir <directory>      Specify the directory containing the HKY trees
  --dbprefix <prefix>         Specify the prefix for the database files
  -o [STR]                    Output file prefix (default: beanOut)
  -fq1 <filename>             Specify the input FASTQ file (single-end or first pair)
  -fq2 <filename>             Specify the input FASTQ file (second pair)
  -t <threads>                Specify the number of threads to use (default: 1)
  -z <directory>              Specify the temporary directory for intermediate files (default: /tmp/)
  -i                          Enable interleaved input mode
Damage options:
  --deam5p                    [.prof] 5p deamination frequency for eukaryotic species (default: no damage)
  --deam3p                    [.prof] 3p deamination frequency for eukaryotic species (default: no damage)
Markov chain Monte Carlo options:
  --no-mcmc                   The MCMC does not run (default: false)
  --chains [INT]              Define the number of chains for the MCMC (default: 4)
  --iter [INT]                Define the number of iterations for the MCMC (default: 500.000)
  --burnin [INT]              Define the burn-in period for the MCMC (default: 75.000)
  --randStart [bool]          Set to get random starting nodes in the tree instead of the signature nodes (default: false)
  -k [INT]                    User defined value of k (k = number of expected sources) (default: not defined)
</pre>

## soibean output files:
soibean's default output consists of a results file and a trace file per k per chain, as well as three diagnostics files per k. All files are tab-separated.

soibean's output can be visualised with multiple R scripts. To ensure they run probably please create an activate soibean's conda environment.
```
conda env create -f ~/vgan/share/vgan/plottingScripts/soibean.yml
conda actiavte soibean
```
The first visualisation step is plot the curve for k:
```
Rscript soibeanPlotk.R [output prefix]
```
The resulting plot shows the highest log-likelihood per k and can give an estimation of the most likely number of sources.

Further evidence for the most likely k can be extracted from the diagnostic files. The main diagnostics file for each k has the following columns:
``` 
Source | Highest log-likelihood | for chain | Rhat for the proportion estimate | Rhat for the branch position estimate
```
For each of the two parameters the MCMC is estimating soibean outputs a seperate diagnostics file. The diagnostics for the proportion estimate:
```
Source | Chain | Mean Proportion Estimate | 5% CI  | Median Proportion Estimate | 95% CI | Effective Sample Size | Autocorrelation Variance
```
and the diagnostics for the branch position estimate:
```
Source | Chain | Mean Branch Position | 5% CI | Median Branch Position | 95% CI | Effective Sample Size | Autocorrelation Variance | Effective Sample Size for the source estimation
```
It is recommended that each of the parameters have an Effective Sample Size of > 200. 

After deciding on the most likely k and the best chain the trace file can be plotted using:
```
Rscript soibeanPlotTrace.R [Tracefile.mcmc] 
```
The results can be plotted with the following command:
```
Rscript soibeanPlots.R [tree file] [Resultsfile.mcmc]
```
The tree file can be found in ```$HOME/share/vgan/soibean_dir/tree_dir```.

## soibean notes:
- Please remember, if the graph file and its indexes are used by multiple users, the graph.dist file needs to be writeable for all. You can use the following command:
```
chmod +w [graph.dist]

```

# TrailMix

## TrailMix notes

## TrailMix Options

<pre>

Algorithm parameters:
        -fq1 [STR]                           FASTQ input file
        -fq2 [STR]                           FASTQ second input file (for paired-end reads)
        -g [STR]                             GAM input file
        -i                                   Input FASTQ (-fq1) is interleaved
        -k                                   Number of distinct contributing haplogroups
        -o [STR]                             Output file (default: stdout)
Non-algorithm parameters:
        -s [STR]                             Sample name
        --dbprefix <prefix>                  Specify the prefix for the database files
        -t                                   Number of threads
        -v                                   Verbose mode
        -z                                   Temporary directory (default: /tmp/)
Markov chain Monte Carlo options:
        --chains [INT]                       Define the number of chains for the MCMC (default: 4)
        --iter [INT]                         Define the number of iterations for the MCMC (default: 1.000.000)
        --randStart [bool]                   Set to get random starting nodes in the tree instead of the signature nodes (default: false)
        --burnin [INT]                       Define the burn-in period for the MCMC (default: 100.000)
Initialization options:
        --mu [INT,INT,...]                   Define the fragment length mean per source (for read count proportion estimation)
        --library-type [STR]                 Strand-specific library type (fr: read1 forward, rf: read1 reverse) (default: unstranded)

</pre>

# Unit tests:

We provide comprehensive unit tests for our main subcommands (HaploCart, euka and soibean). These are built on the Boost Test Suite. Running the tests requires building vgan from source because it involves compiling a separate executable. To run the unit tests, please run

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

or

```
make test
./../bin/test --run_test=Soibean
```

# General notes:
- Please be aware that the multithreading for paired-end input files does not work with the vg version we use. However, from our experience, vg giraffe often reverts to mapping single-end, which is multi-threaded. If you have paired-end reads please consider merging them or concatenate them with a file descriptor:
```
vgan euka -fq1 <(zcat reads1.fq.gz reads2.fq.gz) -t 30 

```

# Support:

Vgan is actively supported and maintained.
If you discover an issue with the program, please submit a bug report on GitHub or send an email. 

Contact:<br>

  Haplocart: jdru@dtu.dk or gabriel.reno@gmail.com
  euka: navo@dtu.dk or gabriel.reno@gmail.com
  soibean: navo@dtu.dk or gabriel.reno@gmail.com
<br>
