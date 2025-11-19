#!/bin/bash

VENV="/trainfs/data/comp_bio/apps/chai-lab-0.6.1/venv"
PROF="/trainfs/data/comp_bio/conda.profile"

module load MMseqs2

source ${PROF}
conda activate ${VENV}

MSA_SCRIPT="/trainfs/data/comp_bio/apps/chai-lab-0.6.1/compute_msa.sh"
BASE_DB="/trainfs/data/comp_bio/data/CHAIDB/"

help() {
echo "Usage: chai-lab fold [OPTIONS] FASTA_FILE OUTPUT_DIR"
echo ""
echo "  Run Chai-1 to fold a complex."
echo ""
echo "Arguments:"
echo "  FASTA_FILE  [required]"
echo "  OUTPUT_DIR  [required]"
echo ""
echo "Options:"
echo "  --use-esm-embeddings / --no-use-esm-embeddings"
echo "                                  [default: use-esm-embeddings]"
echo "  --use-msa"
echo "                                  [default: False]"
echo "  --msa-output-directory PATH"
echo "  --constraint-path PATH"
echo "  --use-templates"
echo "                                  [default: False]"
echo "  --template-hits-path PATH"
echo "  --recycle-msa-subsample INTEGER"
echo "                                  [default: 0]"
echo "  --num-trunk-recycles INTEGER    [default: 3]"
echo "  --num-diffn-timesteps INTEGER   [default: 200]"
echo "  --num-diffn-samples INTEGER     [default: 5]"
echo "  --num-trunk-samples INTEGER     [default: 1]"
echo "  --seed INTEGER"
echo "  --device TEXT"
echo "  --low-memory / --no-low-memory  [default: low-memory]"
echo "  --help                          Show this message and exit."
}

TMPMSADIR=`mktemp -d`
MSA_OUTPUT_DIRECTORY="${TMPMSADIR}"

FASTA_FILE=""
get_fasta_file() {
ARGS=("$@")
while [[ $# -gt 0 ]]; do
    case $1 in
        --use-esm-embeddings|--no-use-esm-embeddings)
            # Handle boolean flags
            shift
            ;;
        --use-templates|--use-msa)
            shift
            ;;
        --low-memory|--no-low-memory)
            shift
            ;;
        --constraint-path|--recycle-msa-subsample|--num-trunk-recycles|--num-diffn-timesteps|--num-diffn-samples|--num-trunk-samples|--seed|--device)
            # These options take a value
            shift 2
            ;;
        --msa-output-directory)
            MSA_OUTPUT_DIRECTORY="$2"
            shift 2  # Remove both --msa-directory and its value
            ;;
        --help)
            # Show help and exit
            shift
            ;;
        -*)
            echo "Unknown option: $1"
            exit 1
            ;;
        *)
            # Positional argument
            positional_args+=("$1")
            shift
            ;;
    esac
done

positional_args=($(remove_arg_from_params "fold" "0" "${positional_args[@]}"))

if [ ${#positional_args[@]} -eq 2 ]; then
    FASTA_FILE="${positional_args[0]}"
    OUTPUT_DIR="${positional_args[1]}"
else
    echo "Error: FASTA_FILE and OUTPUT_DIR positional arguments are required."
    exit 255
fi

if [[ -d "${OUTPUT_DIR}" ]]; then
    echo "the output directory '${OUTPUT_DIR}' already exists"
    exit 255
fi  

set -- "${ARGS[@]}"
}

remove_arg_from_params() {
    local arg_to_remove="$1"
    local do_skip="$2"
    shift 2
    local skip="0"

    local new_args=()
    for arg in "$@"; do
        if [ "$arg" != "$arg_to_remove" ] && [ "$skip" = "0" ]; then
            new_args+=("$arg")
        fi
        skip=0
        if [ "$arg" = "$arg_to_remove" ] && [ "$do_skip" = "1" ]; then
            skip=1
            continue
        fi
    done

    # Print the new argument list
    printf '%s\n' "${new_args[@]}"
}

check_arg_exists() {
    local arg_to_find="$1"
    shift
    for arg in "$@"; do
        if [ "$arg" = "$arg_to_find" ]; then
            return 0  # Found
        fi
    done
    return 1  # Not found
}

ADDITIONAL_ARGS=""
if [ "$1" = "fold" ]; then
    
    MSA_ARGS=""
    if check_arg_exists "--help" "$@"; then
        help
        exit 0
    fi
    
    get_fasta_file "$@"
    if ! [[ -f "$FASTA_FILE" ]]; then
        echo -e "FASTA file '${FASTA_FILE}' does not exist."
        exit 255
    fi

    if check_arg_exists "--use-msa" "$@"; then
        echo "MSA_DIR :: ${MSA_OUTPUT_DIRECTORY}"
        ADDITIONAL_ARGS="${ADDITIONAL_ARGS} --msa-directory ${MSA_OUTPUT_DIRECTORY}"
        MSA_ARGS="${MSA_ARGS} -e -M"
    fi

    if check_arg_exists "--use-templates" "$@"; then
        ADDITIONAL_ARGS="${ADDITIONAL_ARGS} --template-hits-path ${MSA_OUTPUT_DIRECTORY}/all_chain_templates.m8"
        MSA_ARGS="${MSA_ARGS} -T"
    fi

    if [ "$MSA_ARGS" != "" ]; then
        echo "EXECUTING :: ${MSA_SCRIPT}" -q "${FASTA_FILE}" -d "${BASE_DB}" -b "${MSA_OUTPUT_DIRECTORY}" -G ${MSA_ARGS}
        "${MSA_SCRIPT}" -q "${FASTA_FILE}" -d "${BASE_DB}" -b "${MSA_OUTPUT_DIRECTORY}" -G ${MSA_ARGS}
        for file in "${MSA_OUTPUT_DIRECTORY}"/*; do
            if [[ -d "${file}" ]]; then
                chai-lab a3m-to-pqt "${file}" --output-directory "${MSA_OUTPUT_DIRECTORY}"
            fi
        done

        if [[ -f "${MSA_OUTPUT_DIRECTORY}/all_chain_templates.m8" ]]; then
            sed -i 's/\x00//g' "${MSA_OUTPUT_DIRECTORY}/all_chain_templates.m8"
        fi
    fi

    ARGS=($(remove_arg_from_params "--msa-output-directory" "1" "$@"))
    ARGS=($(remove_arg_from_params "--use-msa" "0" "${ARGS[@]}"))
    ARGS=($(remove_arg_from_params "--use-templates" "0" "${ARGS[@]}"))

    echo "chai-lab ${ARGS[@]}" $ADDITIONAL_ARGS
fi

chai-lab "${ARGS[@]}" $ADDITIONAL_ARGS

trap "rm -rf '$TMPMSADIR'" EXIT
