#RunMappeRosendScheffe.sh
#Syntax: d maxLeavesEst n dataSeed critLeaves maxCheck"

#rm *.txt

D=1; #dimension
MAXLEAVESEST=100; #maximum #leaves in the estimator
MAXCHECK=1;
INTN=1000000;
S=1;
N=100;

for L in 10 20 30 
do
	echo $L $S
	./MappedScheffeKL $D $MAXLEAVESEST $N $S $L $MAXCHECK $INTN
	cat MappedKLIAE$S.txt >> KLresults.txt
	rm MappedKLIAE$S.txt
	S=`expr "$S" + 1`; 	
done

#for I in 10000
#do
#	for L in 10
#	do 
	    #mkdir ${MYDIR}/n${I}L${L}
#		for S in `seq $START $TOTAL`
#		do 
#		echo Simulation $TOTAL for n = $I and L = $L
#		#mv *.txt ${MYDIR}/n${I}L${L}
#		done
#	done
#done

