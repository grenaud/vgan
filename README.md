[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?styl=flat)](http://bioconda.github.io/recipes/vgan/README.html)
# vgan

vgan is a suite of tools for pangenomics. We currently support four main subcommands: [HaploCart](#haplocart) (for modern human mtDNA haplogroup classification), [euka](#euka) (for bilaterian taxa detection and abundance estimation of ancient environmental DNA), [soibean](#soibean) (for species identification and multiple source estimation of filtered bilateria data from ancient environmental DNA) and [TrailMix](#trailmix) (for the deconvolution and phylogenetic placement of hominin sedaDNA mixtures). The underlying data structure is the VG graph (see https://github.com/vgteam/vg).

**vgan is supported for use on Linux systems.**

## Contents

1. [Installation](#installation)
2. [euka](#euka)
3. [HaploCart](#haplocart)
4. [soibean](#soibean)
5. [TrailMix](#trailmix)
6. [aeDNA example workflow](#aedna-example-workflow)
7. [Unit tests](#unit-tests)
8. [Support](#support)

---

## Installation

Paths below assume the vgan folder structure in `$HOME/share/vgan/` (created in step 3 of the release build). If you built vgan from source, the same folders are in `share/vgan/` of your clone.

### Release build (recommended)

The easiest way to run vgan is to download the static binary.

**Step 1**: Download the static binary.

The list of releases is here: https://github.com/grenaud/vgan/tags. Each release comes with a static binary. Find the URL of the static binary by right-clicking and selecting "copy link"

```
wget [URL TO RELEASE BINARY]
```

Where you paste the URL of the binary, copy the executable in a bin/ directory in your home directory:

```
mkdir -p $HOME/bin/
cp vgan $HOME/bin/
```

**Step 2**: Mark the binary executable:

```
chmod +x $HOME/bin/vgan
```

**Step 3**: Create the vgan folder structure:

```
mkdir -p $HOME/share/vgan/euka_dir/
mkdir -p $HOME/share/vgan/hcfiles/
mkdir -p $HOME/share/vgan/soibean_dir/
mkdir -p $HOME/share/vgan/tmfiles/
mkdir -p $HOME/share/vgan/damageProfiles/
mkdir -p $HOME/share/vgan/plottingScripts/
```

Once the folder structure is created, download the database files of the subcommands you want to use, see [Database files](#database-files).

### Bioconda

Make sure conda is installed or follow the conda install guide: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html. The Bioconda channel must be installed in your environment.

To install vgan, use:

```
conda install -c bioconda vgan
```

When using vgan via conda, you also must download the database files, see [Database files](#database-files). You can find detailed instructions for each subcommand in its section below.

### Database files

The database files are hosted on Google Drive as `.tar.gz` archives. Each archive unpacks into one folder named like the archive. That folder belongs in `$HOME/share/vgan/`, or anywhere else if you point the subcommand to it with `--euka_dir`, `--hc-files`, `--soibean_dir` or `--tm-files`.

| Archive | Size | Used by | Link |
|---|---|---|---|
| `euka_dir.tar.gz` | 1.3 GB | [euka](#euka) | [download](https://drive.google.com/file/d/1cHDZujMQaMiJxZHhSD56BqlziJWthhPk/view?usp=sharing) |
| `hcfiles.tar.gz` | 91 MB | [HaploCart](#haplocart) | [download](https://drive.google.com/file/d/1gPue0IKdEz3UpY2KMbAHnXZXKvEcCdmQ/view?usp=sharing) |
| `soibean_dir.tar.gz` | 459 MB | [soibean](#soibean) | [download](https://drive.google.com/file/d/1Hzjs4Bc8V71OYp8vY7waI2Q2x4vi8fSP/view?usp=sharing) |
| `tmfiles.tar.gz` | 558 MB | [TrailMix](#trailmix) | [download](https://drive.google.com/file/d/1jPs3artDM_tpKSo0WQgDrqnV9WIk3_Nh/view?usp=sharing) |
| `publication_tmfiles.tar.gz` | 2.4 MB | [TrailMix](#trailmix) (database of the publication) | [download](https://drive.google.com/file/d/1VrFP9ImM7ocN5HR0ppVp2Kp2aVmUo7fX/view?usp=sharing) |

You can download in a browser from these links. To download on the command line, each section below has a `curl` command for its archive; `gdown [file id]` (`pip install gdown`) works as well. The `.tar.gz` can be deleted after unpacking.

### Building vgan from source

Please ensure that git is installed on your system.
You can check by typing "git --version".

If you wish to build from source, you will need the build-essential package. If it is not already installed, please run

```
sudo apt-get install build-essential
```

You will also need wget. If you do not have it please run

```
sudo apt-get install wget
```

vgan requires vg https://github.com/vgteam/vg. vg requires several packages to be installed; please refer to the README in vg.

```
git clone https://github.com/grenaud/vgan.git
cd vgan/src && make
```

Please be patient, this will take some time.

If you are building on a WSL subsystem you may receive an error due to VG not building entirely.
This is likely to be a permissions issue. Please run the above command as a root user. (i.e. `cd vgan/src && sudo make`)

If the build has been successful, the executable will be found in the bin folder.

If you build from source, you can download the requisite files for HaploCart (from the `src` directory) with

```
make hcfilesmade
```

and for euka:

```
make eukafilesmade
```

If euka is to be used by multiple users, please ensure that the file ```$HOME/share/vgan/euka_dir/euka_db.dist``` has writing permission with the following command (from the `src` directory):

```
sudo chmod +w $HOME/share/vgan/euka_dir/euka_db.dist
```

and for soibean:

```
make soibeanfilesmade
```

### Potential issues

vgan requires VG as a dependency. If there are issues with making vg, please see https://github.com/vgteam/vg for further documentation. vg requires several packages to install.

---

## euka

### Introduction

euka is a tool for characterizing ancient environmental DNA samples for arthropodic and tetrapodic mitochondrial DNA. Its foundation is a curated database of complete mitochondrial genomes sorted into taxonomic groups and built into a variation graph reference structure. euka takes FASTQ input files and maps them against our reference graph structure. We use an ancient-aware maximum-likelihood framework to analyze our mapping results on a per-fragment basis. Additionally, we filter on mapping quality and coverage evenness across the pangenome graph of a taxon.

### euka download instructions

Please make sure that you created the vgan folder structure (see [Installation](#installation)). This is necessary for vgan to automatically find the graph input files. If you would like to store the graph input files in a different folder structure, please use the ```--euka_dir``` option to let the program know where your graph input files are.

```
cd $HOME/share/vgan
curl -L -o euka_dir.tar.gz "https://drive.usercontent.google.com/download?id=1cHDZujMQaMiJxZHhSD56BqlziJWthhPk&export=download&confirm=t"
tar -xzf euka_dir.tar.gz && rm euka_dir.tar.gz
wget -O $HOME/share/vgan/damageProfiles/none.prof https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/damageProfiles/none.prof
```

If euka is to be used by multiple users, please ensure that the file ```$HOME/share/vgan/euka_dir/euka_db.dist``` has writing permission with the following command:

```
sudo chmod +w $HOME/share/vgan/euka_dir/euka_db.dist
```

Additionally, you can download euka's visualization scripts:

```
wget -O $HOME/share/vgan/plottingScripts/make_tree_from_output.py https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/make_tree_from_output.py
wget -O $HOME/share/vgan/plottingScripts/plot_taxon.R https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/plot_taxon.R
wget -O $HOME/share/vgan/plottingScripts/visualize_detected_taxa.sh https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/visualize_detected_taxa.sh
wget -O $HOME/share/vgan/plottingScripts/euka.yml https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/euka.yml
```

You can use the ```euka.yml``` file to create a conda environment for plotting. The conda environment includes all necessary packages to create euka's output plots.

### euka quick start

For a metagenomic FASTQ file:

```
wget https://github.com/grenaud/vgan/raw/main/test/input_files/euka/three_dhigh_100.fq.gz
vgan euka -fq1 three_dhigh_100.fq.gz -t 20
```

For a full workflow example, please see the [aeDNA example workflow](#aedna-example-workflow).

### euka usage

For merged or single-end FASTQ (gzipped okay):

```
vgan euka -fq1 [FASTQ file]
```

Please note that if you would like to input multiple FASTQ files, we encourage the use of file descriptors:

```
vgan euka -fq1 <(zcat [FASTQ file1] [FASTQ file2] [FASTQ file3] ...)
```

For paired-end FASTQ (gzipped okay):

```
vgan euka -fq1 [FASTQ file] -fq2 [FASTQ file]
```

euka can take different minimizer files for mapping, which increase/decrease sensitivity and specificity (for a detailed explanation, please see [SAFARI, Rubin, Joshua, et al.](https://doi.org/10.1101/2024.08.12.607489)).
We generally recommend the ```euka_db.rec.min``` minimizer file, which provides the best ratio of sensitivity and specificity. With very noisy datasets, we recommend the ```euka_db.con.min``` minimizer file, which is more specific, but less sensitive. The ```euka_db.sen.min``` file is very sensitive, but can be very unspecific:

```
vgan euka -fq1 [FASTQ file] -M [minimizer file]
```

For incorporation of ancient damage profiles during the maximum likelihood estimation:

```
vgan euka -fq1 [FASTQ file] --deam5p [5end.prof] --deam3p [3end.prof]
```

For user-specific MCMC runs:

```
vgan euka -fq1 [FASTQ file] --iter [INT] --burnin [INT]
```

For additional output options:

```
vgan euka -fq1 [FASTQ file] --outFrag --outGroup [Taxa name]
```

Please note that the taxon name provided to the ```--outGroup``` argument must be a taxon defined in euka's database.

If you specify the ```--outFrag``` option you will be provided with a list of all read names that belong to the detected taxa (one taxon per row). To extract these reads for further downstream analysis you can use the following command:

```
less -S [output file prefix]_FragNames.tsv | sed '[row number]!d' | sed 's/\t/\n/g' | seqtk subseq [FASTQ input file] - > output.fq
```

To make euka less strict in detecting taxa you can use the ```--maxBins``` flag. The default is 0, meaning no bin in the coverage estimation is allowed to have no reads mapping to it. If you set it higher, it corresponds to the number of bins that are allowed to be empty (have no reads mapping to them).

```
vgan euka -fq1 [FASTQ file] --maxBins 2
```

### euka filter options

euka has four filter options to modify the stringency of taxa detection. We always recommend having a first look at your samples with the default parameters. These parameters have been thoroughly tested to provide confident abundance estimations. However, to detect more divergent taxa (for example Formicidae), it may be necessary to adjust the filter parameters. Furthermore, we want to highlight that the reference genomes for many arthropodic species have low-complexity reference genomes and are more prone to spurious alignments. Even with our standard parameter filters, we can see more false-positive detections for these taxa, and results should always be evaluated carefully. The file `$HOME/share/vgan/euka_dir/euka_db.bins` lists every taxon in our database with its respective bins (for our coverage estimation). The file shows the Node ID range for each bin and the calculated entropy score for this bin.

Example for incorporation of lower-entropy regions in the mitogenome:

```
vgan euka -fq1 [FASTQ file] --entropy 1.13 --minBins 3
```

Example for more conservative parameter settings:

```
vgan euka -fq1 [FASTQ file] --minMQ 40 --minFrag 20
```

### euka options

A list of all of euka's options:

```
Input options:
        --euka_dir [STR]        euka database location (default: current working directory)
        --dbprefix [STR]        database prefix name (default: euka_db)
        -fq1 [STR]              Input FASTQ file (for merged and single-end reads)
        -fq2 [STR]              Second input FASTQ file (for paired-end reads)
        -i                      Paired-end reads are interleaved (default: false)
        -g [STR]                GAM input file
        -M [STR]                Alternative minimizer prefix input (default: euka_db)
        -o [STR]                Output file prefix (default: euka_output)
        -t                      Number of threads (-1 for all available)
        -Z                      Temporary directory (default: /tmp)

Filter options:
        --minMQ [INT]           Set the mapping quality minimum for a fragment (default: 29)
        --minFrag [INT]         Minimum amount of fragments that need to map to a group (default: 10)
        --entropy [double]      Minimum entropy score for a bin to be considered (default: 1.17)
        --minBins [INT]         Minimum number of bins with an entropy higher than the threshold (default: 6)
        --maxBins [INT]         Maximum number of bins without coverage (default: 0)

Damage options:
        --deam5p [.prof]        5p deamination frequency for eukaryotic species (default: no damage)
        --deam3p [.prof]        3p deamination frequency for eukaryotic species (default: no damage)
        -l [INT]                Set length for substitution matrix (default: 5)
        --out_dir [STR]         Path for output prof-file (default: current working directory)

Markov chain Monte Carlo options:
        --no-mcmc               The MCMC does not run (default: false)
        --iter [INT]            Define the number of iterations for the MCMC (default: 10000)
        --burnin [INT]          Define the burn-in period for the MCMC (default: 100)

Additional output option:
        --outFrag               Outputs a file with all read names per taxonomic group (default: false)
        --outGroup [string]     Outputs all information about a taxonomic group of interest (default: empty)
```

### euka output files

euka's default output consists of four TSV files:

The abundance and the detected TSV file consists of the following columns:

```
| Taxa | detected | Number of reads | proportion estimate | 85% CI lower bound | 85% CI upper bound | 95% CI lower bound | 95% CI upper bound |
```

The abundance file will list all available taxa, while the detected file will only list taxa that have passed all of euka's detection filters.

There are multiple ways to visualize euka's output. To make sure the provided scripts work, we recommend installing the provided conda environment:

```
conda env create -f $HOME/share/vgan/plottingScripts/euka.yml
```

Or make sure the following packages are installed for Python: ete3, csv and argparse; for R, you will need the libraries ggplot2 and ggpubr.

You have the option to visualize these two files with a taxonomic tree using the provided Python script (make sure Python is in your path):

```
python $HOME/share/vgan/plottingScripts/make_tree_from_output.py [output file prefix]_abundance.tsv
python $HOME/share/vgan/plottingScripts/make_tree_from_output.py [output file prefix]_detected.tsv
```

The coverage file consists of all detected taxa and their corresponding number of bins. It shows the number of reads sorted into each bin. The coverage TSV file is used by our plotting script to show a barplot.

The inSize file provides a list of all fragment sizes for each taxon. It is used by our plotting script to create a histogram of the fragment length distribution.

For each detected taxon, we estimate a damage profile in a .prof file. These files are used as input for the damage estimation plots in our plotting script.
Additionally, we provide an average damage profile for the 5' end and the 3' end that can be directly inputted into our damage model options ```--deam5p``` and ```--deam3p```.

To visualize detected taxa, you can use the following provided script:

```
$HOME/share/vgan/plottingScripts/visualize_detected_taxa.sh [output file prefix]
```

To visualize only a specific detected taxon, please use:

```
Rscript $HOME/share/vgan/plottingScripts/plot_taxon.R [output file prefix] [taxon name]
```

### euka notes

- We recommend running euka first with standard parameters and no incorporated damage profiles.
- Further, we encourage the use of all available input files for one ecological site. Especially with ancient environmental DNA samples where read counts/abundances are extremely low, euka will have difficulties "detecting" a taxon if, for example, only one sequencing lane is provided.
- Please remember, if the graph file and its indexes are used by multiple users, the euka graph.dist file needs to be writeable for all. You can use the following command:

```
chmod +w euka_db.dist
```

---

## HaploCart

### Introduction

HaploCart performs maximum likelihood estimation to predict the mitochondrial haplogroup for reads
originating from an uncontaminated modern human sample. The program also optionally provides confidence estimation of the phylogenetic placement of the sample.

### HaploCart download instructions

Please download HaploCart's necessary graph input files. The files will be automatically read in by vgan if you constructed the vgan folder structure. If not, please remember to use ```--hc-files``` and specify the path to the input files.

```
cd $HOME/share/vgan
curl -L -o hcfiles.tar.gz "https://drive.usercontent.google.com/download?id=1gPue0IKdEz3UpY2KMbAHnXZXKvEcCdmQ&export=download&confirm=t"
tar -xzf hcfiles.tar.gz && rm hcfiles.tar.gz
```

### HaploCart quick start

For consensus FASTA:

```
wget https://github.com/grenaud/vgan/raw/main/test/input_files/rCRS.fa
vgan haplocart -f rCRS.fa
```

For single-end FASTQ:

```
wget https://github.com/grenaud/vgan/raw/main/test/input_files/rCRS.fq
vgan haplocart -fq1 rCRS.fq
```

### HaploCart usage

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

HaploCart does not yet accept HTS-lib-compatible input directly. Please first unmap the reads and then follow the instructions for FASTQ input described above, e.g.

```
samtools bam2fq [FASTQ file] | vgan haplocart -fq1 /dev/stdin
```

VCF input:

We generally advise against using VCF as input as this will rely on a model that was used for calling the variants. Furthermore, a lot of the information is lost and we therefore recommend to use it on the full dataset to gain maximum advantage from the pangenome graph. However, if you wish to use it anyway, we provide a helper script to convert to fasta. If the rCRS (with header J01415.2 in the fasta) was used to call the variants, then use:

```
share/vgan/hc_scripts/vcf2fasta.py input.vcf.gz rCRS J01415.2 | gzip > input.fa.gz
vgan haplocart -f input.fa.gz
```

### HaploCart output file

The output file for HaploCart is a TSV file with the following columns:

```
| Sample | Predicted haplogroup | # Mapped reads |
```

### HaploCart posterior file

The output file for HaploCart's optional posterior probability estimation is a TSV file with the following columns:

```
| Sample | Considered clade | Posterior probability | Tree depth |
```

### HaploCart notes

- HaploCart will make a prediction on any sample for which at least one read maps. It is up to the user to apply appropriate filters to HaploCart results.
If your input data is especially sparse please check the clade-level posterior probabilities and treat the predictions judiciously.

- Although the HaploCart inference algorithm is deterministic, we rely upon VG Giraffe for mapping and we have noticed a degree of stochasticity in the mapping procedure.
  Therefore in some edge cases there may be a (very slight) difference in output with different runs on the same sample. Unfortunately Giraffe has no random seed, so we cannot provide one either.

- Please remember, if the graph file and its indexes are used by multiple users the graph.dist file needs to be writeable for all. You can use the following command:

```
chmod +w [graph.dist]
```

---

## soibean

### Introduction

soibean is a tool crafted for identifying and placing ancient environmental DNA from **pre-analyzed FASTQ files** onto a phylogenetic tree. Pre-analyzed FASTQ files are all FASTQ reads that have been **previously classified to a lower taxonomic resolution** (e.g. family level). It leverages mitochondrial pangenome graphs across 335 arthropodic and tetrapodic taxa, drawing on the euka database enhanced with ancestral states and complemented by HKY model phylogenetic trees. soibean maps input FASTQs to a specific subgraph, determining the likely species or species mix within. It uses an MCMC algorithm to estimate the species' positions on the tree and their proportions in the sample.

### soibean download instructions

We recommend following the vgan folder structure explained in the [Installation](#installation) section. However, if you prefer to store the soibean parent graph and its subgraphs in different directories, use the ```--soibean_dir``` option to specify the path to the graph files and the ```--tree_dir``` option to specify the path to the directory storing the tree files. Please remember that the ```--dbprefix``` option for the subgraph is mandatory for soibean.

```
cd $HOME/share/vgan
curl -L -o soibean_dir.tar.gz "https://drive.usercontent.google.com/download?id=1Hzjs4Bc8V71OYp8vY7waI2Q2x4vi8fSP&export=download&confirm=t"
tar -xzf soibean_dir.tar.gz && rm soibean_dir.tar.gz
wget -O $HOME/share/vgan/damageProfiles/none.prof https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/damageProfiles/none.prof
```

The archive includes the tree files in `tree_dir/` (already unzipped; `trees.zip` is the same set) and the subgraph files of the bears (`Ursidae.*`).

To download all additional scripts for soibean, please use the following commands:

```
wget -O $HOME/share/vgan/soibean_dir/make_graph_files.sh https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/soibean_dir/make_graph_files.sh
chmod +x $HOME/share/vgan/soibean_dir/make_graph_files.sh
wget -O $HOME/share/vgan/plottingScripts/soibeanPlotTrace.R https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/soibeanPlotTrace.R
wget -O $HOME/share/vgan/plottingScripts/soibeanPlots.R https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/soibeanPlots.R
wget -O $HOME/share/vgan/plottingScripts/soibeanPlotk.R https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/soibeanPlotk.R
wget -O $HOME/share/vgan/plottingScripts/soibean.yml https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/soibean.yml
```

### soibean quick start

For FASTQ mapped to the taxon of bears (Ursidae):

```
wget https://github.com/grenaud/vgan/raw/main/test/input_files/soibean/k1.fq.gz
vgan soibean -fq1 k1.fq.gz --dbprefix Ursidae -t 20
```

The subgraph of the bears is part of the download. For other taxa, create the subgraph first with `make_graph_files.sh` (see below).

### soibean usage

Extract a subgraph of interest:

```
$HOME/share/vgan/soibean_dir/make_graph_files.sh [taxon name] [--vg-path <path to vg>]
```

For example, to extract the graph for the taxon of bears (Ursidae):

```
$HOME/share/vgan/soibean_dir/make_graph_files.sh Ursidae
```

You will get all the necessary graph files for mapping with soibean as output. You only have to create these files once per taxon.

For merged or single-end FASTQ (gzipped okay):

```
vgan soibean -fq1 [FASTQ file] --dbprefix [taxon name]
```

Please note that if you would like to input multiple FASTQ files, we encourage the use of file descriptors:

```
vgan soibean -fq1 <(cat [FASTQ file1] [FASTQ file2] [FASTQ file3] ...) --dbprefix [taxon name]
```

For paired-end FASTQ (gzipped okay):

```
vgan soibean -fq1 [FASTQ file] -fq2 [FASTQ file] --dbprefix [taxon name]
```

For incorporation of ancient damage profiles during the maximum likelihood estimation:

```
vgan soibean -fq1 [FASTQ file] --dbprefix [taxon name] --deam5p [5end.prof] --deam3p [3end.prof]
```

Examples of damage files are in the repository: ```share/vgan/damageProfiles/dhigh5p.prof``` and ```share/vgan/damageProfiles/dhigh3p.prof```.
euka provides damage profiles for each detected taxon ending in ```.prof```. The taxon-specific damage profiles can be adjusted to be used as soibean input:

```
line=$(grep -n "^A>C" [euka taxon output].prof | head -2 | tail -1 | cut -d: -f1) && awk -v line=$line -v OFS='\t' 'BEGIN{FS=OFS} NR < line {NF--; print}' [euka taxon output].prof > [taxon 5end].prof && (head -1 [euka taxon output].prof && awk -v line=$line -v OFS='\t' 'BEGIN{FS=OFS} NR >= line {if(NR > line) {NF--; print}}' [euka taxon output].prof | tac) > [taxon 3end].prof
```

For user-defined options in the MCMC, e.g., the number of iterations, the number of burn-in iterations, or the maximum number of sources (k):

```
vgan soibean -fq1 [FASTQ file] --dbprefix [taxon name] --iter [INT] --burnin [INT] -k [INT]
```

or

```
vgan soibean -fq1 [FASTQ file] --dbprefix [taxon name] --iter [INT] --burnin [INT] --randStart true
```

### soibean options

A list of all of soibean's options:

```
Options:
  --soibean_dir <directory>   Specify the directory containing the soibean files
  --tree_dir <directory>      Specify the directory containing the HKY trees
  --dbprefix <prefix>         Specify the prefix for the database files
  -o [STR]                    Output file prefix (default: beanOut)
  -fq1 <filename>             Specify the input FASTQ file (single-end or first pair)
  -fq2 <filename>             Specify the input FASTQ file (second pair)
  -g <filename>               Specify the input GAM file (SAFARI output)
  -t <threads>                Specify the number of threads to use (default: 1)
  -z <directory>              Specify the temporary directory for intermediate files (default: /tmp/)
  -i                          Enable interleaved input mode
  -M <filename>               Specify an alternative minimizer index (default: <prefix>.rec.min)
  -P [int]                    Penalty for unsupported paths (default: 7)
Damage options:
  --deam5p [.prof]            5p deamination frequency for eukaryotic species (default: no damage)
  --deam3p [.prof]            3p deamination frequency for eukaryotic species (default: no damage)
Markov chain Monte Carlo options:
  --no-mcmc                   The MCMC does not run (default: false)
  --chains [INT]              Define the number of chains for the MCMC (default: 4)
  --iter [INT]                Define the number of iterations for the MCMC (default: 500.000)
  --burnin [INT]              Define the burn-in period for the MCMC (default: 75.000)
  --randStart [bool]          Set to get random starting nodes in the tree instead of the signature nodes (default: false)
  -k [INT]                    User defined value of k (k = number of expected sources) (default: not defined)
Additional output options:
  --alignment-detail [bool]   Additional TSV file with all match alignment information to identify SNPs (default: false)
  --pathThres [INT]           Reports all matches in additional TSV file for nodes with a number of paths less or equal to the threshold set (default: all paths)
```

### soibean output files

soibean's default output consists of a results file and a trace file per k per chain, as well as three diagnostics files per k. All files are tab-separated.

soibean's output can be visualized with multiple R scripts. To ensure they run properly please create and activate soibean's conda environment.

```
conda env create -f $HOME/share/vgan/plottingScripts/soibean.yml
conda activate soibean
```

The first visualization step is to plot the curve for k:

```
Rscript $HOME/share/vgan/plottingScripts/soibeanPlotk.R [output prefix]
```

The resulting plot shows the highest log-likelihood per k and can give an estimation of the most likely number of sources.

Further evidence for the most likely k can be extracted from the diagnostic files. The main diagnostics file for each k has the following columns:

```
Source | Highest log-likelihood | for chain | Rhat for the proportion estimate | Rhat for the branch position estimate
```

For each of the two parameters the MCMC is estimating, soibean outputs a separate diagnostics file. The diagnostics for the proportion estimate:

```
Source | Chain | Mean Proportion Estimate | 5% CI | Median Proportion Estimate | 95% CI | Effective Sample Size | Autocorrelation Variance
```

and the diagnostics for the branch position estimate:

```
Source | Chain | Mean Branch Position | 5% CI | Median Branch Position | 95% CI | Effective Sample Size | Autocorrelation Variance | Effective Sample Size for the source estimation
```

It is recommended that each of the parameters have an Effective Sample Size of > 200.

After deciding on the most likely k and the best chain the trace file can be plotted using:

```
Rscript $HOME/share/vgan/plottingScripts/soibeanPlotTrace.R [Tracefile.mcmc]
```

The results can be plotted with the following command:

```
Rscript $HOME/share/vgan/plottingScripts/soibeanPlots.R [tree file] [Resultsfile.mcmc]
```

The tree file can be found in ```$HOME/share/vgan/soibean_dir/tree_dir```.

### soibean notes

- Please remember, if the graph file and its indexes are used by multiple users, the graph.dist file needs to be writeable for all. You can use the following command:

```
chmod +w [graph.dist]
```

---

## TrailMix

### Introduction

TrailMix deconvolves a mixture of hominin mitochondrial DNA from ancient sediment (sedaDNA) into k contributing sources. For each source it estimates the proportion in the mixture and its position on a phylogenetic tree of mitochondrial genomes, using an MCMC. It takes FASTQ or GAM input.

### TrailMix download instructions

The default database is `tmfiles`; vgan looks for it in `$HOME/share/vgan/tmfiles/` (or the folder given with `--tm-files`).

```
cd $HOME/share/vgan
curl -L -o tmfiles.tar.gz "https://drive.usercontent.google.com/download?id=1jPs3artDM_tpKSo0WQgDrqnV9WIk3_Nh&export=download&confirm=t"
tar -xzf tmfiles.tar.gz && rm tmfiles.tar.gz
```

The database of the publication is a separate archive. Its files start with `pub.graph`, so it needs `--dbprefix pub.graph` in addition to `--tm-files`:

```
cd $HOME/share/vgan
curl -L -o publication_tmfiles.tar.gz "https://drive.usercontent.google.com/download?id=1VrFP9ImM7ocN5HR0ppVp2Kp2aVmUo7fX&export=download&confirm=t"
tar -xzf publication_tmfiles.tar.gz && rm publication_tmfiles.tar.gz
```

To download the plotting script (it needs Python with numpy, pandas, biopython, matplotlib and scipy):

```
wget -O $HOME/share/vgan/plottingScripts/plot_trailmix_trees.py https://raw.githubusercontent.com/grenaud/vgan/main/share/vgan/plottingScripts/plot_trailmix_trees.py
```

### TrailMix usage

For merged or single-end FASTQ (gzipped okay), with two sources:

```
vgan trailmix -fq1 [FASTQ file] -k 2 -o [output prefix]
```

With the database of the publication, and several chains:

```
vgan trailmix -fq1 [FASTQ file] -k 2 --chains 5 --iter 100000 --tm-files $HOME/share/vgan/publication_tmfiles/ --dbprefix pub.graph -o [output prefix]
```

For paired-end FASTQ, add `-fq2 [FASTQ file]`.

### TrailMix options

```
Algorithm parameters:
        -fq1 [STR]              FASTQ input file
        -fq2 [STR]              FASTQ second input file (for paired-end reads)
        -g [STR]                GAM input file
        -k                      Number of distinct contributing haplogroups
        -o [STR]                Output file prefix (default: TM_out)
        --depth [INT]           Tree depth to whitelist the MCMC search space (-1 for no restriction)

Non-algorithm parameters:
        --dbprefix <prefix>     Prefix of the database files (default: graph)
        --tm-files [STR]        TrailMix file directory (default: "../share/vgan/tmfiles/")
        -t                      Number of threads
        -v                      Verbose mode
        -z                      Temporary directory (default: /tmp/)

Markov chain Monte Carlo options:
        --chains [INT]          Number of chains for the MCMC (default: 1)
        --iter [INT]            Number of iterations for the MCMC (default: 10000)
        --burnin [INT]          Burn-in period for the MCMC (default: 100)
        --randStart [bool]      Random starting nodes in the tree instead of the signature nodes (default: false)

Initialization options:
        --library-type [STR]    Strand-specific library type (fr: read1 forward, rf: read1 reverse) (default: unstranded)

Contamination mode options:
        --contamination-mode, --cont-mode
                                Requires exactly two sources. Accepted for compatibility; it no longer changes the likelihood (use --is-ancient)
        --is-ancient [BOOL_VECTOR]
                                One 0/1 per source (-k entries), comma-delimited, where 1 means ancient and 0 modern (e.g. 1,0). Default: all 1.
```

### TrailMix output files

`-o` is a filename prefix; the results are written next to it. For `k` sources and chain `N`:

- `[prefix]Result[k]_chain[N].mcmc.gz`: one row per MCMC iteration after the burn-in, with the columns `Source_i`, `Node_i`, `Loglikelihood`, `proportion_i` and `branch_position_derived_i` for each source `i`. The position is the fraction of the way from the parent to the node along the branch named by `Source_i`.
- `[prefix]Trace[k]_chain[N].detail.mcmc.gz`: every proposal of the chain, with its acceptance probability and move.
- `[prefix]BranchEstimate.txt` and `[prefix]BaseProportionEstimates.txt`: the estimates per source and chain, with effective sample sizes. `Source` is the most visited branch of the slot; the statistics cover every iteration after the burn-in; `5% CI` and `95% CI` are the 5th and 95th percentiles.
- `[prefix]Diagnostics.txt`: R-hat for the proportions and for the placements.
- `[prefix]RPVGHaplotypes.txt` and `[prefix]RPVGHaplotypeTranscripts.txt`: the haplotypes RPVG finds in the reads before the MCMC runs.
- `[prefix]summary.txt`: `key: value` lines with the sample, k, chains, iterations, burn-in and the database that was used.

The sources are exchangeable, so the slot numbers can swap between chains. Compare chains (R-hat, effective sample size and the plots) before trusting an estimate.

To plot the results on the tree:

```
python $HOME/share/vgan/plottingScripts/plot_trailmix_trees.py [output prefix]
```

It reads `[prefix]summary.txt` for the database and the chains, draws the posterior placements of every source on the tree and the proportions, and marks the best-mixing chain. The options are described at the top of the script.

---

## aeDNA example workflow

A possible workflow for your ancient environmental DNA sample could be:

```
wget https://github.com/grenaud/vgan/raw/main/test/input_files/euka/three_dhigh_100.fq.gz

vgan euka -fq1 three_dhigh_100.fq.gz --outFrag
```

The taxon of bears (Ursidae) has been detected. We plot the taxon to verify that it looks ancient:

```
$HOME/share/vgan/plottingScripts/visualize_detected_taxa.sh euka_output
```

After the verification, we extract the reads mapped to the bear taxon to do species detection (this step needs [seqtk](https://github.com/lh3/seqtk)). `[row number]` is the row of the bear taxon in `euka_output_FragNames.tsv`:

```
less -S euka_output_FragNames.tsv | sed '[row number]!d' | sed 's/\t/\n/g' | seqtk subseq three_dhigh_100.fq.gz - | gzip > Ursidae.fq.gz
```

We extract the subgraph for the bears from our soibean graph (you can skip this step for the bears, because the `Ursidae.*` files are part of the soibean download):

```
$HOME/share/vgan/soibean_dir/make_graph_files.sh Ursidae
```

The script extracts the bear subgraph and creates all necessary files for mapping and analyzing with soibean:

```
vgan soibean -fq1 Ursidae.fq.gz --dbprefix Ursidae -t 20 -o soibeanOut
```

Once soibean is done, we can plot our results. First, we check the k-curve, which shows the number of sources best describing our data:

```
Rscript $HOME/share/vgan/plottingScripts/soibeanPlotk.R soibeanOut
```

Afterwards, we can choose the best run of our MCMC and plot the distribution on the phylogenetic tree (the tree files can be found in ```$HOME/share/vgan/soibean_dir/tree_dir```):

```
Rscript $HOME/share/vgan/plottingScripts/soibeanPlots.R Ursidae.new.dnd soibeanOutResult10.mcmc
```

---

## Unit tests

We provide comprehensive unit tests for our main subcommands (HaploCart, euka and soibean). These are built on the Boost Test Suite. Running the tests requires building vgan from source because it involves compiling a separate executable. To run all unit tests, please run

```
cd vgan/src
make test
../bin/test
```

You can also run subcommand-specific tests:

```
../bin/test --run_test=haplocart
```

or

```
../bin/test --run_test=euka
```

or

```
../bin/test --run_test=Soibean
```

---

## Support

vgan is actively supported and maintained.
If you discover an issue with the program, please submit a bug report on GitHub or send an email.

Contact:

- HaploCart: jdru@dtu.dk or gabriel.reno@gmail.com
- euka: n.alexandra.vogel@gmail.com or gabriel.reno@gmail.com
- soibean: n.alexandra.vogel@gmail.com or gabriel.reno@gmail.com
