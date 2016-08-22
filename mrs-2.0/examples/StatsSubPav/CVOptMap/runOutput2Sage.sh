#!/bin/bash
FILE="$1"
#cat "$FILE"
# take the output of run.sh and parse it for getting mean, min, max and std in sage using numpy
cat "$FILE" | tail -11 | head -10 | sed 's/\t/,/g' | sed 's/ //g' | tr '\n' '; ' | sed '$s/;$/\n/' | sed -e "s/'/'\\\\''/g;s/\(.*\)/'\1'/" | sed '$i import numpy as np; np.set_printoptions(suppress=True); a=np.matrix(' | sed '$a ); print a.mean(0); print a.min(0); print a.max(0); print a.std(0);\n'
# now paste the output into $sage shell
