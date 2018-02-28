#!/bin/bash
# mmlong toolbox
# By Rasmus Kirkegaard and Søren Karst
# Version
VERSION=0.1.0

################################################################################
### Preparation ----------------------------------------------------------------
################################################################################

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [ name ...]

-- mmlong toolbox v. $VERSION: Misc. tools for mmlong workflow.

where:
    -h   Show this help text.
    name Name of tool.
    ...  Commands for tool.

Available tools:
* metaflow: De novo assembly of metagenome using nanopore data, read coverage
  estimation based on nanopore and illumina data, aswell as taxonomic classification
  and detection of essntial genes and 16S. For more info type: \`mmlong metaflow -h\`
* reassembly: Reassembly of metagenome bins using unicycler. For more info type: 
  \`mmlong reassembly -h\`
* miniasm-racon: Basic nanopore assembly with miniasm followed by racon polishing.
  For more info type: \`mmlong miniasm-racon -h\`
* checkassembly: align assembly to metagenome bin and map bin reads to
  metagenome bin. Output alignments and read mappings in SAM format for
  visualization in CLC. For more info type: \`mmlong checkassembly -h\`
"

### Terminal Arguments ---------------------------------------------------------


# Import user arguments
while getopts ':hz' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done
TOOL=$1
TOOL_ARG="${*:2}"

# Check missing arguments
if [ -z "${TOOL}" ]; then printf "\n No tool selected!\n\n"; echo "$USAGE"; exit 1; fi; 
if [ -z "${TOOL_ARG}" ]; then
  printf "\n No arguments provided to $TOOL!\n\n"
  TOOL_ARG="-h"
  echo "$USAGE"
  exit 1
fi

# Paths
SCRIPT_PATH=`readlink ${BASH_SOURCE[0]}`
export MMLONG_DIR=${SCRIPT_PATH%/*}

### Call tool ------------------------------------------------------------------

$MMLONG_DIR/mmlong_tools/mmlong_$TOOL $TOOL_ARG
