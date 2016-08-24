#!/bin/bash
FILE="$1"
#cat "$FILE"
# take the output of run.sh and parse it for getting mean, min, max and std in sage using numpy
cat "$FILE" | tail -10 | head -10 | sed 's/\t/,/g' | sed 's/ //g' | tr '\n' '; ' | sed '$s/;$/\n/' | sed -e "s/'/'\\\\''/g;s/\(.*\)/'\1'/" | sed '$i import numpy as np; np.set_printoptions(suppress=True); a=np.matrix(' | sed '$a ); print "----------------L1 dist and d_KL--------"; print a.mean(0); print a.min(0); print a.max(0); print a.std(0); print "----------------TV dist and sqrt(d_KL/2)--------"; b=a.copy(); b[:,0]=a[:,0]/2.0; b[:,1]=np.sqrt(a[:,1]/2.0); print b.mean(0); print b.min(0); print b.max(0); print b.std(0);\n'
# now paste the output into $sage shell
