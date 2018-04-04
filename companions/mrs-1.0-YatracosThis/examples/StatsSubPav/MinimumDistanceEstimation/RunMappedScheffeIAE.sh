#RunMappedScheffeIAE.sh
#Syntax: d maxLeavesEst n dataSeed critLeaves maxCheck"

rm *.txt

TOTAL=1; #How many simulations

D=1; #dimension
MAXLEAVESEST=100; #maximum #leaves in the estimator
MAXCHECK=1;
STARTLEAVES=5;

#MYDIR=RepMappedGaussian${D}D
#mkdir $MYDIR 
#mkdir ${MYDIR}/n${N}L${CRITLEAVES}

for I in 100
do
	for L in 10
	do 
	    #mkdir ${MYDIR}/n${I}L${L}
		for S in `seq 1 $TOTAL`
		do 
		echo Simulation $TOTAL for n = $I and L = $L
		./MappedScheffeIAE $D $MAXLEAVESEST $I $S $L $MAXCHECK $STARTLEAVES;
		#mv *.txt ${MYDIR}/n${I}L${L}
		done
	done
done
