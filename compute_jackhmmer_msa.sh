#!/bin/bash -e

# Script to compute MSA for all entries in a FASTA file using jackhmmer
# Runs searches in parallel across multiple databases

JACKHMMER="jackhmmer-patch"
MMSEQS="mmseqs"
PDBDB="pdb100_230517"

# Usage function
usage() {
    echo "Usage: $0 -q <QUERY> -d <database_dir> -b <BASE> [-c <cores>] [-M] [-e] [-T]"
    echo "  -q    Input FASTA file with protein sequences"
    echo "  -d    Base directory containing sequence databases"
    echo "  -b    Output directory for MSA results"
    echo "  -c    Number of cores to use for parallel processing"
    echo "  -M    Compute MSA using jackhmmer"
    echo "  -e    Compute environmental sequences using jackhmmer"
    echo "  -T    Compute templates using MMseqs2"
    exit 1
}

COMPUTE_MSA="0"
COMPUTE_ENV="0"
COMPUTE_TEMPLATE="0"
CORES=1

# Split FASTA file into individual sequences
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

# Parse command line arguments
while getopts "q:c:d:b:MeTh" opt; do
    case $opt in
        q) QUERY="$OPTARG" ;;
        d) DBBASE="$OPTARG" ;;
        b) BASE="$OPTARG" ;;
        c) CORES="$OPTARG" ;;
        h) usage ;;
        M) COMPUTE_MSA="1" ;;
        e) COMPUTE_ENV="1" ;;
        T) COMPUTE_TEMPLATE="1" ;;
        *) usage ;;
    esac
done

NUMBER_DB=$((COMPUTE_MSA * 2 + COMPUTE_ENV))

# Check required arguments
if [ -z "$QUERY" ] || [ -z "$DBBASE" ] || [ -z "$BASE" ]; then
    usage
fi

if ! [[ -f "$QUERY" ]]; then
    echo "Error: FASTA file '$QUERY' does not exist."
    exit 1
fi

if [[ -f "$BASE" ]]; then
    echo "Error: Output path '$BASE' is a file."
    exit 1
fi

if ! [[ -d "$DBBASE" ]]; then
    echo "Error: Database directory '$DBBASE' does not exist."
    exit 1
fi

if [ "$CORES" -lt ${NUMBER_DB} ]; then
    echo "Error: At least ${NUMBER_DB} cores are required for the requested databases (compute MSA: ${COMPUTE_MSA}, compute ENV: ${COMPUTE_ENV})."
    exit 1
fi

# Database paths (modify these variables as needed)
UNIPROT_DB="${DBBASE}/uniprot.fasta"
UNIREF90_DB="${DBBASE}/uniref90.fasta"
MGNIFY_DB="${DBBASE}/mgy_clusters_2022_05.fasta"

# Create output directory
mkdir -p "$BASE"

awk '/^>/ {if (seq) {print seq > file; seq=""}; file=sprintf("'$TEMP_DIR'/seq_%04d.fasta", ++count); print > file; next} {seq = seq $0} END {if (seq) print seq > file}' "$QUERY"

# Function to run jackhmmer on a single sequence against one database
run_jackhmmer() {
    local seq_file=$1
    local db_path=$2
    local db_name=$3
    local seq_limit=$4
    local seq_id=$(grep "^>" "$seq_file" | head -1 | sed 's/^>protein|//;s/ .*//')
    local seq_output_dir="${BASE}/${seq_id}"
    
    mkdir -p "$seq_output_dir"
    ${JACKHMMER} -A "${seq_output_dir}/hits_${db_name}.a3m" --seq_limit "${seq_limit}" --F1 0.0005 --F2 0.00005 --F3 0.0000005 --incE 0.0001 -E 0.0001 --cpu $((CORES / NUMBER_DB)) -N 1 "$seq_file" "$db_path"
}

# Process all sequences in parallel
for file in "$TEMP_DIR"/seq_*.fasta; do
    if [ "${COMPUTE_MSA}" = "1" ]; then
        run_jackhmmer "$file" "$UNIPROT_DB" "uniprot" 50000 &
        run_jackhmmer "$file" "$UNIREF90_DB" "uniref90" 10000 &
    fi
    if [ "${COMPUTE_ENV}" = "1" ]; then
        run_jackhmmer "$file" "$MGNIFY_DB" "mgnify" 5000 &
    fi

    wait
done

echo "THIS IS BASE ${BASE}"
if [ "${COMPUTE_TEMPLATE}" = "1" ]; then
    module load MMseqs2
    echo "Searching for templates using MMseqs2"
    TEMPLATE_SEARCH_PARAM="--db-load-mode 2 --threads ${CORES} -s 7.5 -a -e 0.0001 --gpu 1"
    "${MMSEQS}" createdb "${QUERY}" "${BASE}/qdb"
    "${MMSEQS}" search "${BASE}/qdb" "${DBBASE}/${PDBDB}" "${BASE}/res_pdb" "${BASE}/tmp" $TEMPLATE_SEARCH_PARAM
    "${MMSEQS}" convertalis "${BASE}/qdb" "${DBBASE}/${PDBDB}" "${BASE}/res_pdb" "${BASE}/all_chain_templates.m8" --threads ${CORES} --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,cigar --db-load-mode 2 --db-output 1
    "${MMSEQS}" rmdb "${BASE}/res_pdb"
    "${MMSEQS}" rmdb "${BASE}/qdb"
    "${MMSEQS}" rmdb "${BASE}/qdb_h"
    "${MMSEQS}" rmdb "${BASE}/res"
    rm -f "${BASE}"/*.dbtype
    rm -f "${BASE}"/*.index
    rm -rf -- "${BASE}/tmp"
    
    sed -i "s/protein[|]//" "${BASE}/all_chain_templates.m8"
fi

echo "MSA computation complete. Results in: $BASE"