#!/bin/bash/

module purge
module load group ce-molsim stack python

path="./up/"
species="MeOH"
framework="IRMOF1step0"
postfilename="01.full.postfile"

PYTHONPATH=$PYTHONPATH:~/scratch/pymbar/
export PYTHONPATH

echo "Setup all fine, beginning analysis."

python SerialResultsParser.py -p $path $species $framework $postfilename > analysis.log
echo "Analysis complete!"
