#!/bin/bash

module purge
module load group ce-molsim stack python

path="./fullpostfiles/"
species="Acetone_TRAPPE"
framework="IRMOF1"
tag="postoutput"


echo "Setup all fine, beginning analysis."

python ParallelResultsParser.py -p $path $species $framework $tag > analysis.log
echo "Analysis complete!"
