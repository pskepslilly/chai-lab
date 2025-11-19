#!/bin/bash -e

# Exit handler to provide feedback
trap 'echo "Error: MSA Script failed at line $LINENO with exit code $?" >&2' ERR

export OMP_NUM_THREADS=100
MMSEQS="mmseqs"

SENSITIVITY=8

MGNIFY="mgnify_gpu_db"
UNIPROT="uniprot_gpu_db"
UNIREF90="uniref90_gpu_db"
PDB100="pdb100_230517"

FILTER=1
GPU=0

COMPUTE_ENV=0
COMPUTE_TEMPLATE=0
COMPUTE_MSA=0

# Default values
MMSEQS="mmseqs"
DBBASE=""
BASE=""
QUERY=""
THREADS=1
UNIREF90="uniref90_gpu_db"
PDB100="pdb100_230517"
MGNIFY="mgnify_gpu_db"
UNIPROT="uniprot_gpu_db"

EXPAND_EVAL=inf
ALIGN_EVAL=10
DIFF=3000
QSC=0.8
MAX_ACCEPT=1000000

# Parse command line arguments
while getopts "m:q:d:b:t:u:p:g:n:x:a:GD:Q:A:eTMf:h" opt; do
    case $opt in
        m) MMSEQS="$OPTARG" ;;
        q) QUERY="$OPTARG" ;;
        d) DBBASE="$OPTARG" ;;
        b) BASE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        u) UNIREF90="$OPTARG" ;;
        p) PDB100="$OPTARG" ;;
        g) MGNIFY="$OPTARG" ;;
        n) UNIPROT="$OPTARG" ;;
        x) EXPAND_EVAL="$OPTARG" ;;
        a) ALIGN_EVAL="$OPTARG" ;;
        G) GPU="1" ;;
        D) DIFF="$OPTARG" ;;
        Q) QSC="$OPTARG" ;;
        A) MAX_ACCEPT="$OPTARG" ;;
        e) COMPUTE_ENV="1" ;;
        T) COMPUTE_TEMPLATE="1" ;;
        M) COMPUTE_MSA="1" ;;
        f) FILTER="$OPTARG" ;;
        h)
            echo "Usage: $0 -q QUERY -d DBBASE -b BASE [OPTIONS]"
            echo "Required:"
            echo "  -q QUERY            Input query file"
            echo "  -d DBBASE           Database base path"
            echo "  -b BASE             Output base path"
            echo "Optional:"
            echo "  -m MMSEQS           MMseqs2 binary (default: mmseqs)"
            echo "  -t THREADS          Number of threads (default: 1)"
            echo "  -u UNIREF90         UniRef90 database (default: ${UNIREF90})"
            echo "  -p PDB100           PDB100 database (default: ${PDB100})"
            echo "  -g MGNIFY           MGnify database (default: ${MGNIFY})"
            echo "  -n UNIPROT          UniProt database (default: ${UNIPROT})"
            echo "  -x EXPAND_EVAL      Expand eval value (default: ${EXPAND_EVAL})"
            echo "  -a ALIGN_EVAL       Align eval value (default: ${ALIGN_EVAL})"
            echo "  -G GPU              Use GPU acceleration (default: ${GPU})"
            echo "  -D DIFF             Diff value (default: ${DIFF})"
            echo "  -Q QSC              QSC value (default: ${QSC})"
            echo "  -A MAX_ACCEPT       Max accept value (default: ${MAX_ACCEPT})"
            echo "  -e COMPUTE_ENV      Compute environment (default: ${COMPUTE_ENV})"
            echo "  -T COMPUTE_TEMPLATE Compute templates (default: ${COMPUTE_TEMPLATE})"
            echo "  -M COMPUTE_MSA      Compute MSA (default: ${COMPUTE_MSA})"
            echo "  -f FILTER           Filter value (default: ${FILTER})"
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

if [ "${GPU}" = "0" ]; then
    echo "Currently only supports GPU acceleration"
    exit 1
fi

add_taxid_to_a3m() {
    local raw_a3m=$1
    local taxid_tsv=$2
    local output_a3m=$3
    
    awk -v taxfile="${taxid_tsv}" '
    BEGIN {
        # Load TaxID mapping from file
        while ((getline line < taxfile) > 0) {
            split(line, fields, "\t")
            if (length(fields) >= 2) {
                taxid[fields[1]] = fields[2]
            }
        }
        close(taxfile)
    }
    /^>/ {
        # Extract sequence ID (first word after >)
        seqid = substr($1, 2)
        
        # Get TaxID, default to "" if not found
        tid = (seqid in taxid) ? taxid[seqid] : ""
        
        # Print header with TaxID

        if (tid != "") {
            printf ">%s TaxID=%s", seqid, tid   
        }
        else {
            printf ">%s", seqid
        }
        
        # Add rest of original header if present
        for (i = 2; i <= NF; i++) {
            printf " %s", $i
        }
        printf "\n"
        next
    }
    {
        # Print sequence lines as-is
        print
    }
    ' "${raw_a3m}" > "${output_a3m}"
}

extract_protein_sequences() {
    local input_fasta=$1
    local output_fasta=$2
    
    awk '
    /^>protein\|/ {
        printing = 1
    }
    /^>/ && !/^>protein\|/ {
        printing = 0
    }
    printing {
        print
    }
    ' "${input_fasta}" > "${output_fasta}"
}


# Check required arguments
if [ -z "$QUERY" ] || [ -z "$DBBASE" ] || [ -z "$BASE" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 -q QUERY -d DBBASE -b BASE [OPTIONS]"
    echo "Use -h for help"
    exit 1
fi

extract_protein_sequences "${QUERY}" "${TMPDIR}/query_proteins.fasta"
QUERY="${TMPDIR}/query_proteins.fasta"

# Check if query file has any sequences
if [ ! -s "${QUERY}" ]; then
    echo "Error: No protein sequences found in input file"
    exit 1
fi

export MMSEQS_CALL_DEPTH=1

SEARCH_PARAM="--threads ${THREADS} --num-iterations 3 --db-load-mode 2 -a -e 0.0001 -s ${SENSITIVITY} --gpu ${GPU}"
TEMPLATE_SEARCH_PARAM="--db-load-mode 2 --threads ${THREADS} -s 7.5 -a -e 0.0001 --gpu ${GPU}"
if [ "${GPU}" = "1" ]; then
    SEARCH_PARAM="${SEARCH_PARAM} --prefilter-mode 1"
    TEMPLATE_SEARCH_PARAM="${TEMPLATE_SEARCH_PARAM} --prefilter-mode 1"
fi
FILTER_PARAM="--filter-msa ${FILTER} --filter-min-enable 1000 --diff ${DIFF} --qid 0.0,0.2,0.4,0.6,0.8,1.0 --qsc 0 --max-seq-id 0.95"
EXPAND_PARAM="--expansion-mode 0 -e ${EXPAND_EVAL} --expand-filter-clusters ${FILTER} --max-seq-id 0.95"
mkdir -p "${BASE}"
"${MMSEQS}" createdb "${QUERY}" "${BASE}/qdb"

if [ "${COMPUTE_MSA}" = "1" ]; then
    "${MMSEQS}" search "${BASE}/qdb" "${DBBASE}/${UNIREF90}" "${BASE}/res" "${BASE}/tmp" $SEARCH_PARAM --max-seqs 10000
    "${MMSEQS}" mvdb "${BASE}/tmp/latest/profile_1" "${BASE}/prof_res"
    "${MMSEQS}" lndb "${BASE}/qdb_h" "${BASE}/prof_res_h"
    # "${MMSEQS}" expandaln "${BASE}/qdb" "${DBBASE}/${UNIREF90}.idx" "${BASE}/res" "${DBBASE}/${UNIREF90}.idx" "${BASE}/res_exp" --db-load-mode 2 --threads ${THREADS} ${EXPAND_PARAM}
    "${MMSEQS}" align "${BASE}/qdb" "${DBBASE}/${UNIREF90}" "${BASE}/res" "${BASE}/res_realign" --db-load-mode 2 --threads ${THREADS} -e ${ALIGN_EVAL} --max-accept ${MAX_ACCEPT} --alt-ali 10 -a
    "${MMSEQS}" filterresult "${BASE}/qdb" "${DBBASE}/${UNIREF90}" "${BASE}/res_realign" "${BASE}/res_realign_filter" --db-load-mode 2 --threads ${THREADS} --qid 0 --qsc $QSC --diff 0 --max-seq-id 1.0 --filter-min-enable 100
    "${MMSEQS}" result2msa "${BASE}/qdb" "${DBBASE}/${UNIREF90}" "${BASE}/res_realign_filter" "${BASE}/raw_hits_uniref90.a3m" --msa-format-mode 6 --db-load-mode 2 --threads ${THREADS} ${FILTER_PARAM}
    "${MMSEQS}" convertalis "${BASE}/qdb" "${DBBASE}/${UNIREF90}" "${BASE}/res_realign_filter" "${BASE}/uniref90_taxid_mapping.tsv" --format-output "target,taxid" --threads ${THREADS} --db-load-mode 2 --db-output 0
    add_taxid_to_a3m "${BASE}/raw_hits_uniref90.a3m" "${BASE}/uniref90_taxid_mapping.tsv" "${BASE}/hits_uniref90.a3m"
    "${MMSEQS}" rmdb "${BASE}/prof_res"
    "${MMSEQS}" rmdb "${BASE}/prof_res_h"
    "${MMSEQS}" rmdb "${BASE}/res_realign"
    # "${MMSEQS}" rmdb "${BASE}/res_exp"
    "${MMSEQS}" rmdb "${BASE}/res"
    "${MMSEQS}" rmdb "${BASE}/res_realign_filter"

    "${MMSEQS}" search "${BASE}/qdb" "${DBBASE}/${UNIPROT}" "${BASE}/res" "${BASE}/tmp" $SEARCH_PARAM --max-seqs 50000
    "${MMSEQS}" mvdb "${BASE}/tmp/latest/profile_1" "${BASE}/prof_res"
    "${MMSEQS}" lndb "${BASE}/qdb_h" "${BASE}/prof_res_h"
    # "${MMSEQS}" expandaln "${BASE}/qdb" "${DBBASE}/${UNIPROT}.idx" "${BASE}/res" "${DBBASE}/${UNIPROT}.idx" "${BASE}/res_exp" --db-load-mode 2 --threads ${THREADS} ${EXPAND_PARAM}
    "${MMSEQS}" align "${BASE}/qdb" "${DBBASE}/${UNIPROT}" "${BASE}/res" "${BASE}/res_realign" --db-load-mode 2 --threads ${THREADS} -e ${ALIGN_EVAL} --max-accept ${MAX_ACCEPT} --alt-ali 10 -a
    "${MMSEQS}" filterresult "${BASE}/qdb" "${DBBASE}/${UNIPROT}" "${BASE}/res_realign" "${BASE}/res_realign_filter" --db-load-mode 2 --threads ${THREADS} --qid 0 --qsc $QSC --diff 0 --max-seq-id 1.0 --filter-min-enable 100
    "${MMSEQS}" result2msa "${BASE}/qdb" "${DBBASE}/${UNIPROT}" "${BASE}/res_realign_filter" "${BASE}/raw_hits_uniprot.a3m" --msa-format-mode 6 --db-load-mode 2 --threads ${THREADS} ${FILTER_PARAM}
    "${MMSEQS}" convertalis "${BASE}/qdb" "${DBBASE}/${UNIPROT}" "${BASE}/res_realign_filter" "${BASE}/uniprot_taxid_mapping.tsv" --format-output "target,taxid" --threads ${THREADS} --db-load-mode 2 --db-output 0
    add_taxid_to_a3m "${BASE}/raw_hits_uniprot.a3m" "${BASE}/uniprot_taxid_mapping.tsv" "${BASE}/hits_uniprot.a3m"
    "${MMSEQS}" rmdb "${BASE}/res_realign"
    # "${MMSEQS}" rmdb "${BASE}/res_exp"
    "${MMSEQS}" rmdb "${BASE}/res"
    "${MMSEQS}" rmdb "${BASE}/res_realign_filter"
fi

if [ "${COMPUTE_TEMPLATE}" = "1" ]; then
    query="${BASE}/qdb"
    if [[ -f "${BASE}/prof_res" ]]; then
        query="${BASE}/prof_res"
    fi
  "${MMSEQS}" search "${query}" "${DBBASE}/${PDB100}" "${BASE}/res_pdb" "${BASE}/tmp" $TEMPLATE_SEARCH_PARAM
  "${MMSEQS}" convertalis "${query}" "${DBBASE}/${PDB100}" "${BASE}/res_pdb" "${BASE}/all_chain_templates.m8" --threads ${THREADS} --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,cigar --db-load-mode 2 --db-output 1
  "${MMSEQS}" rmdb "${BASE}/res_pdb"
fi

if [ "${COMPUTE_ENV}" = "1" ]; then
  "${MMSEQS}" search "${BASE}/qdb" "${DBBASE}/${MGNIFY}" "${BASE}/res_env" "${BASE}/tmp" $SEARCH_PARAM --max-seqs 5000
#   "${MMSEQS}" expandaln "${BASE}/prof_res" "${DBBASE}/${MGNIFY}.idx" "${BASE}/res_env" "${DBBASE}/${MGNIFY}.idx" "${BASE}/res_env_exp" -e ${EXPAND_EVAL} --threads ${THREADS} --expansion-mode 0 --db-load-mode 2
  "${MMSEQS}" align "${BASE}/tmp/latest/profile_1" "${DBBASE}/${MGNIFY}" "${BASE}/res_env" "${BASE}/res_env_realign" --db-load-mode 2 --threads ${THREADS} -e ${ALIGN_EVAL} --max-accept ${MAX_ACCEPT} --alt-ali 10 -a
  "${MMSEQS}" filterresult "${BASE}/qdb" "${DBBASE}/${MGNIFY}.idx" "${BASE}/res_env_realign" "${BASE}/res_env_realign_filter" --db-load-mode 2 --threads ${THREADS} --qid 0 --qsc $QSC --diff 0 --max-seq-id 1.0 --filter-min-enable 100
  "${MMSEQS}" result2msa "${BASE}/qdb" "${DBBASE}/${MGNIFY}.idx" "${BASE}/res_env_realign_filter" "${BASE}/hits_mgnify.a3m" --msa-format-mode 6 --db-load-mode 2 --threads ${THREADS} ${FILTER_PARAM}
  "${MMSEQS}" rmdb "${BASE}/res_env_realign_filter"
  "${MMSEQS}" rmdb "${BASE}/res_env_realign"
#   "${MMSEQS}" rmdb "${BASE}/res_env_exp"
  "${MMSEQS}" rmdb "${BASE}/res_env"
fi
"${MMSEQS}" rmdb "${BASE}/qdb"
"${MMSEQS}" rmdb "${BASE}/qdb_h"
"${MMSEQS}" rmdb "${BASE}/res"
rm -f -- "${BASE}/prof_res"*
rm -rf -- "${BASE}/tmp"


# POST PROCESSING:
# Separate a3m files into chunks based on protein headers
for a3m_file in "${BASE}/hits_mgnify.a3m" "${BASE}/hits_uniprot.a3m" "${BASE}/hits_uniref90.a3m"; do
    if [ -f "$a3m_file" ]; then
        basename=$(basename "$a3m_file")
        # Replace >protein| with > in the output files
        awk -v base="${BASE}" -v fname="$basename" '
            BEGIN { RS="\x00"; FS="\n" }
            {
                if (match($1, /^>protein\|(.*)/, arr)) {
                    if (arr[1] != "") {
                        dir = base "/" arr[1]
                        system("mkdir -p \"" dir "\"")
                        outfile = dir "/" fname
                        # Replace >protein| with > in the output
                        gsub(/^>protein\|/, ">", $1)
                        # make sure to remove TaxID if present
                        gsub(/TaxID=[0-9]+/, "", $1)
                        
                        # Add newline to every line in the record
                        for (i = 1; i <= NF; i++) {
                            print $i >> outfile
                        }
                    }
                }
            }
        ' "$a3m_file"
    fi
done

# need to strip the extra null characters left in as well
if [[ -f "${BASE}/all_chain_templates.m8" ]]; then
    sed -i "s/protein[|]//" "${BASE}/all_chain_templates.m8"
fi

# FINAL CLEANUP
rm -f "${BASE}"/*.dbtype
rm -f "${BASE}"/*.index
rm -f "${BASE}"/hits_mgnify.a3m "${BASE}"/hits_uniprot.a3m "${BASE}"/hits_uniref90.a3m
rm -f "${BASE}"/raw_hits_mgnify.a3m "${BASE}"/raw_hits_uniprot.a3m "${BASE}"/raw_hits_uniref90.a3m
rm -f "${BASE}"/uniprot_taxid_mapping.tsv  "${BASE}"/uniref90_taxid_mapping.tsv
