#!/bin/bash

# Function to show help
show_help() {
    echo "Usage: $0 <taxon_name> [options]"
    echo ""
    echo "First argument must be the taxon name."
    echo "Then, options can be specified."
    echo ""
    echo "Options:"
    echo "  -h, --help           Show this help message"
    echo "  --vg-path <path>     Specify the VG executable path (optional)"
    echo ""
    echo "Example:"
    echo "$0 Ursidae --vg-path /usr/local/bin/"
}

# Function to check for vg binary or download it
find_or_download_vg_binary() {
    # If a VG path is specified, check if 'vg' exists there.
    if [[ -n "$VG_EXECUTABLE_PATH" ]]; then
        # Append "/vg" if the path does not end with it, assuming it's a directory.
        if [[ ! "$VG_EXECUTABLE_PATH" =~ vg$ ]]; then
            VG_EXECUTABLE_PATH="${VG_EXECUTABLE_PATH}/vg"
        fi
        # Check if the vg binary exists and is executable.
        if [[ ! -x "$VG_EXECUTABLE_PATH" ]]; then
            echo "VG executable not found at $VG_EXECUTABLE_PATH. Downloading..."
            # Download directly to the specified path (considering VG_EXECUTABLE_PATH includes '/vg').
            wget -O "$VG_EXECUTABLE_PATH" https://github.com/vgteam/vg/releases/download/v1.44.0/vg
            chmod +x "$VG_EXECUTABLE_PATH"
            echo "VG executable downloaded to: $VG_EXECUTABLE_PATH"
        else
            echo "Using VG executable found at $VG_EXECUTABLE_PATH"
        fi
    else
        # No VG path specified, check the current working directory.
        if [[ -x "./vg" ]]; then
            VG_EXECUTABLE_PATH="./vg"
            echo "Using VG executable found in the current directory."
        else
            echo "VG executable not found in the current directory. Downloading..."
            VG_EXECUTABLE_PATH="./vg"
            wget -O "$VG_EXECUTABLE_PATH" https://github.com/vgteam/vg/releases/download/v1.44.0/vg
            chmod +x "$VG_EXECUTABLE_PATH"
            echo "VG executable downloaded to the current directory."
        fi
    fi
}


# Check for help option
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    exit 0
fi

input_name="$1" # The first argument is always the input/taxa name.
shift # Remove the first argument so we can process the rest.

# Initialize VG_EXECUTABLE_PATH variable with an empty value indicating not set.
VG_EXECUTABLE_PATH=""

# Now process the rest of the arguments for options.
while [[ "$1" =~ ^- ]]; do
    case "$1" in
        --vg-path)
            shift # Move past the option to get its value.
            VG_EXECUTABLE_PATH="$1"
            # Ensure the path ends with "vg" or is a directory.
            if [[ ! "$VG_EXECUTABLE_PATH" == */vg && ! -d "$VG_EXECUTABLE_PATH" ]]; then
                VG_EXECUTABLE_PATH="$VG_EXECUTABLE_PATH/vg"
            fi
            ;;
        *)
            show_help
            exit 1
            ;;
    esac
    shift # Move to the next option.
done

# Call the function to check for the VG executable or download it.
find_or_download_vg_binary

echo "Processing input name: $input_name"

SCRIPT_DIR=$(readlink -f "$(dirname "$0")")


CLADE_FILE="${SCRIPT_DIR}/soibean_db.clade"
# Check for the soibean_db.clade file, download if not exists
if [ ! -f "$CLADE_FILE" ]; then
    echo "soibean_db.clade file not found. Downloading..."
    wget -nc --recursive --no-parent -P "${SCRIPT_DIR}" ftp://ftp.healthtech.dtu.dk:/public/soibean_files/soibean_db.clade
    mv ftp.healthtech.dtu.dk/public/soibean_files/soibean_db.clade ${SCRIPT_DIR} 
    rm -r ftp.healthtech.dtu.dk/public/soibean_files/
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

$VG_EXECUTABLE_PATH chunk -r ${fifth_element}:${sixth_element} -x ${SCRIPT_DIR}/soibean_db.og > ${SCRIPT_DIR}/${input_name}_pre.vg
echo "Subgraph created!"
echo "Building indexing files for subgraph ... "
$VG_EXECUTABLE_PATH mod -X 5 ${SCRIPT_DIR}/${input_name}_pre.vg > ${SCRIPT_DIR}/${input_name}.vg
$VG_EXECUTABLE_PATH snarls ${SCRIPT_DIR}/${input_name}.vg > ${SCRIPT_DIR}/${input_name}.snarls
$VG_EXECUTABLE_PATH view ${SCRIPT_DIR}/${input_name}.vg > ${SCRIPT_DIR}/${input_name}.gfa
$VG_EXECUTABLE_PATH gbwt -o ${SCRIPT_DIR}/${input_name}.gbwt -g ${SCRIPT_DIR}/${input_name}.gg -G ${SCRIPT_DIR}/${input_name}.gfa
$VG_EXECUTABLE_PATH index -j ${SCRIPT_DIR}/${input_name}.dist ${SCRIPT_DIR}/${input_name}.vg
$VG_EXECUTABLE_PATH minimizer -g ${SCRIPT_DIR}/${input_name}.gbwt -i ${SCRIPT_DIR}/${input_name}.min -k 20 -w 10 ${SCRIPT_DIR}/${input_name}.vg
$VG_EXECUTABLE_PATH convert -g ${SCRIPT_DIR}/${input_name}.gfa -o > ${SCRIPT_DIR}/${input_name}.og

chmod +w ${SCRIPT_DIR}/${input_name}.dist

echo "... done! You can find your graph files in ${SCRIPT_DIR}!"
