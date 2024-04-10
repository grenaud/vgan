#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: $0 [options] <input_name>"
    echo ""
    echo "Options:"
    echo "  -h, --help               Show this help message and exit."
    echo "  --vg-path <path>         Explicitly specify the path to the VG executable."
    echo "  taxa                     Overview of all available taxa in the soibean database."
    echo ""
    echo "Arguments:"
    echo "  input_name               Taxon of interest. This taxon will be extracted from the soibean database."
    echo ""
    echo "Examples:"
    echo "  $0 taxa                  Provides a list of all available taxa in the soibean database"
    echo "  $0 taxon_name --vg-path /path/to/vg  Extracts the taxon of interest from the soibean database and creates all index files needed for soibean."
    echo ""
}

# Set the script's working directory to $HOME/vgan/share/vgan/soibean_dir/
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"


# No arguments, show help
if [ $# -eq 0 ]; then
    show_help
    exit 0
fi

# Check for help option
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    exit 0
fi

# Initialize VG_EXECUTABLE_PATH variable
VG_EXECUTABLE_PATH=""

# Parse options
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
  --vg-path )
    shift
    # Check if the path ends with "vg" or "/", and adjust accordingly
    if [[ "${1}" == */vg ]]; then
      VG_EXECUTABLE_PATH="${1}"
    elif [[ "${1}" == */ ]]; then
      VG_EXECUTABLE_PATH="${1}vg"
    else
      VG_EXECUTABLE_PATH="${1}/vg"
    fi
    ;;
  * )
    show_help
    exit 1
    ;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi


input_name=$1

# Function to check for vg binary or download it
find_or_download_vg_binary() {
    if [[ -n "$VG_EXECUTABLE_PATH" ]]; then
        if [[ ! -x "$VG_EXECUTABLE_PATH" ]]; then
            # VG path specified but executable not found, download it
            echo "VG executable not found at $VG_EXECUTABLE_PATH. Downloading..."
            wget -O "$VG_EXECUTABLE_PATH" https://github.com/vgteam/vg/releases/download/v1.44.0/vg
            chmod +x "$VG_EXECUTABLE_PATH"
            echo "VG executable downloaded to: $VG_EXECUTABLE_PATH"
        else
            echo "Using specified VG executable path: $VG_EXECUTABLE_PATH"
        fi
        return 0
    fi

    # VG path not specified, search for it
    # (You can add default search locations and logic here)
}

# Call find_or_download_vg_binary to handle VG executable path logic
find_or_download_vg_binary

CLADE_FILE="${SCRIPT_DIR}/soibean_db.clade"
# Check for the soibean_db.clade file, download if not exists
if [ ! -f "$CLADE_FILE" ]; then
    echo "soibean_db.clade file not found. Downloading..."
    wget -nc --recursive --no-parent -P "${SCRIPT_DIR}" ftp://ftp.healthtech.dtu.dk:/public/soibean_files/soibean_db.clade
else
    echo ""
fi

# Subcommand for reading the second column
if [[ $1 == "taxa" ]]; then
  while IFS=' ' read -r -a array
  do
    echo "${array[1]}"
  done < "${SCRIPT_DIR}/soibean_db.clade"
  exit 0
fi


fifth_element=""
sixth_element=""

while IFS= read -r line
do
  if echo "$line" | grep -q "${input_name}"; then
    set -- $line
    fifth_element=$5
    sixth_element=$6
    echo "5th element: $fifth_element"
    echo "6th element: $sixth_element"
    break
  fi
done < "${SCRIPT_DIR}/soibean_db.clade"

# Check that both elements are numbers
if ! [[ $fifth_element =~ ^[0-9]+$ ]] || ! [[ $sixth_element =~ ^[0-9]+$ ]]; then
  echo "5th and 6th elements are not both numbers"
  exit 1
fi

if [[ -z "$fifth_element" || -z "$sixth_element" ]]; then
    echo "Name not found in the file or the line doesn't have at least 6 elements"
    exit 1
fi

$VG_EXECUTABLE_PATH chunk -r ${fifth_element}:${sixth_element} -x soibean_db.og > ${input_name}_pre.vg
echo "Subgraph created!"
echo "Building indexing files for subgraph ... "
$VG_EXECUTABLE_PATH mod -X 5 ${input_name}_pre.vg > ${input_name}.vg
$VG_EXECUTABLE_PATH snarls ${input_name}.vg > ${input_name}.snarls
$VG_EXECUTABLE_PATH view ${input_name}.vg > ${input_name}.gfa
$VG_EXECUTABLE_PATH gbwt -o ${input_name}.gbwt -g ${input_name}.gg -G ${input_name}.gfa
$VG_EXECUTABLE_PATH index -j ${input_name}.dist ${input_name}.vg
$VG_EXECUTABLE_PATH minimizer -g ${input_name}.gbwt -i ${input_name}.min -k 20 -w 10 ${input_name}.vg
$VG_EXECUTABLE_PATH convert -g ${input_name}.gfa -o > ${input_name}.og

chmod +w ${input_name}.dist

echo "... done!"
