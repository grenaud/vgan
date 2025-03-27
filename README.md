[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?styl=flat)](http://bioconda.github.io/recipes/vgan/README.html)
# vgan

vgan is a suite of tools for pangenomics. We currently support three main subcommands: Haplocart (for modern human mtDNA haplogroup classification), euka (for bilaterian abundance estimation of ancient environmental DNA) and soibean (for species identification and multiple source estimation of filtered bilateria data from ancient environmental DNA). The underlying data structure is the VG graph (see https://github.com/vgteam/vg).

You can find the complete **installation guide** and detailed **manuals** for each subcommand in our [**wiki**](https://github.com/grenaud/vgan/wiki).

**vgan is supported for use on Linux systems.**

## Release build (recommended)

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
mkdir -p $HOME/share/vgan/damageProfiles/
mkdir -p $HOME/share/vgan/plottingScripts/
```
Once the folder structure is created, we can download the necessary files. 
Find download instructions as well as detailed manuals for **euka**, **HaploCart**, and **soibean** in our [**wiki**](https://github.com/grenaud/vgan/wiki). 

# Support:

vgan is actively supported and maintained.
If you discover an issue with the program, please submit a bug report on GitHub or send an email. 

Contact:<br>

  Haplocart: jdru@dtu.dk or gabriel.reno@gmail.com 
  euka: n.alexandra.vogel@gmail.com or gabriel.reno@gmail.com
  soibean: n.alexandra.vogel@gmail.com or gabriel.reno@gmail.com
<br>
