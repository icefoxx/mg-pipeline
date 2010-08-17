#! /bin/bash
if [ $# -ne 1 ]
then
    echo "Error: Invalid Argument Count - enter absolute file name or directory for analysis"
    echo "Syntax: $0 absolute_path"
    exit
fi
export PYTHON_PATH="/nfs/seqdb/production/interpro/production/python"
export PYTHON_VERSION="python2.6.4/bin/python"
export SCRIPT_PATH="/nfs/seqdb/production/interpro/development/metagenomics"
export OUTPUT_PATH="/nfs/nobackup/interpro/development/metagenomics/analyses"
$PYTHON_PATH/$PYTHON_VERSION $SCRIPT_PATH/triggerMetagenPipeline.py -v $PYTHON_PATH/$PYTHON_VERSION -o $OUTPUT_PATH -s $SCRIPT_PATH -p $1
