#!/bin/bash
#module load python/3.6.2
#module load py_packages/3.6
#

INPUT_FILE=$1
#BASE_DIR="/sc/orga/projects/Signatures/GlobalPVal"
BASE_DIR="/Users/calina01/PycharmProjects/GlobalPVal"
#BASE_DIR="/Users/anna/PycharmProjects/Pathway_Enrichment_Ranked_Expression"

if [ $# -eq 0 ]
  then
    echo "USAGE: sh analyze_expression_profiles.sh path/to/input"
    exit
fi

if [ ! -f ${INPUT_FILE} ]
then
    echo "Could not find $INPUT_FILE."
    exit
fi

# Read in database choice
read -p 'Database to use (k for KEGG, r for Reactome, s for HMDB/SMPDB): ' DB_CHOICE

if [ ${DB_CHOICE} == 'r' ] || [ ${DB_CHOICE} == 'R' ]; then
    DATABASE='reactome'
elif [ ${DB_CHOICE} == 'k' ] || [ ${DB_CHOICE} == 'K' ]; then
    DATABASE='KEGG'
elif [ ${DB_CHOICE} == 's' ] || [ ${DB_CHOICE} == 'S' ]; then
    DATABASE='HMDB_SMPDB'
else
    DATABASE='KEGG'
    echo 'Database choice not recognized. Defaulting to KEGG.'
fi

# Read in Gene Expression Choice
#read -p 'Genes of interest (u for upregulated, d for downregulated): ' GENE_EXPRESSION_CHOICE
#if [ ${GENE_EXPRESSION_CHOICE} == 'd' ] || [ ${GENE_EXPRESSION_CHOICE} == 'D' ]; then
#    ASCENDING='y'
#elif [ ${GENE_EXPRESSION_CHOICE} == 'u' ] || [ ${GENE_EXPRESSION_CHOICE} == 'U' ]; then
#    ASCENDING='n'
#else
#    ASCENDING='y'
#    echo 'Gene expression choice not recognized. Defaulting to downregulated.'
#fi

# Read in custom pathways
read -p 'Are you using a list of custom pathways (y/n): ' CUSTOM_PATHWAY
if [ ${CUSTOM_PATHWAY} == 'y' ] || [ ${CUSTOM_PATHWAY} == 'Y' ]; then
    read -p "Please name your custom pathway: ${DATABASE}_" CUSTOM_PW_NAME_PART

    CUSTOM_PATHWAY_NAME="${DATABASE}_${CUSTOM_PW_NAME_PART}"

    CUSTOM_PATHWAY_DB="${BASE_DIR}/databases/user_pathways/${CUSTOM_PATHWAY_NAME}.pkl"

    if [ -f ${CUSTOM_PATHWAY_DB} ]
    then
        read -p "${CUSTOM_PATHWAY_NAME} exists in the db already. Use existing (y or n)? " USE_EXISTING

        if [ ${USE_EXISTING} == 'y' ]
        then
            CUSTOM_PATHWAY_EXISTS='y'
            CUSTOM_PATHWAY_FILE=1
        else
            CUSTOM_PATHWAY_EXISTS='n'
            COPY=1
            CUSTOM_PATHWAY_NAME_UPDATED="${BASE_DIR}/databases/user_pathways/${CUSTOM_PATHWAY_NAME}_${COPY}.pkl"

            while [ -f ${CUSTOM_PATHWAY_NAME_UPDATED} ]
            do
                let "COPY++"
                CUSTOM_PATHWAY_NAME_UPDATED="${BASE_DIR}/databases/user_pathways/${CUSTOM_PATHWAY_NAME}_${COPY}.pkl"
            done

            CUSTOM_PATHWAY_NAME="${CUSTOM_PATHWAY_NAME}_${COPY}"
        fi
    else
        CUSTOM_PATHWAY_EXISTS='n'
    fi
else
    CUSTOM_PATHWAY_FILE=0
    CUSTOM_PATHWAY_NAME=0
    CUSTOM_PATHWAY_EXISTS='n'
fi

if [ ${CUSTOM_PATHWAY_EXISTS} == 'n' ] && [ ${CUSTOM_PATHWAY_NAME} != 0 ]
then
    read -e -p 'Path to custom pathway file: ' CUSTOM_PATHWAY_FILE

    if [ ! -f ${CUSTOM_PATHWAY_FILE} ]
    then
        echo "Could not find $CUSTOM_PATHWAY_FILE."
        exit
    fi
fi


#echo "#!/bin/bash
##BSUB -J ${INPUT_FILE}
##BSUB -W 15:00
##BSUB -n 1
##BSUB -q premium
##BSUB -P acc_Signatures
##BSUB -m mothra
##BSUB -o /sc/orga/scratch/calina01/globalpval.out
##BSUB -e /sc/orga/scratch/calina01/globalpval.err
##BSUB -u anna.calinawan@mssm.edu
##BSUB -R rusage[mem=32000]
##BSUB -L /bin/bash
#
#module load python/3.6.2
#module load py_packages/3.6
#
#python ${BASE_DIR}/main.py ${INPUT_FILE} ${DATABASE} ${ASCENDING} ${BASE_DIR} ${CUSTOM_PATHWAY_FILE} ${CUSTOM_PATHWAY_NAME} ${CUSTOM_PATHWAY_EXISTS}"  > "${BASE_DIR}/run_job.lsf"
#
#bsub < "${BASE_DIR}/run_job.lsf"
#ASCENDING='y'
#python ${BASE_DIR}/main.py ${INPUT_FILE} ${DATABASE} ${ASCENDING} ${BASE_DIR} ${CUSTOM_PATHWAY_FILE} ${CUSTOM_PATHWAY_NAME} ${CUSTOM_PATHWAY_EXISTS}


ASCENDING='n'
python ${BASE_DIR}/main.py ${INPUT_FILE} ${DATABASE} ${ASCENDING} ${BASE_DIR} ${CUSTOM_PATHWAY_FILE} ${CUSTOM_PATHWAY_NAME} ${CUSTOM_PATHWAY_EXISTS}
