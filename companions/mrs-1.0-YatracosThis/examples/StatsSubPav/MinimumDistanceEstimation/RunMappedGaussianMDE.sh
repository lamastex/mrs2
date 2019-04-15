#RunMapedGaussianMDE.sh
#Syntax: "....

DATASEED=1
D=1
N=100
MAXLEAVESEST=100
CRITLEAVES=5
MAXCHECK=1

./MappedGaussianMDE $DATASEED $D $N $MAXLEAVESEST $CRITLEAVES $MAXCHECK

# TOTAL=2; #How many simulations

# D=1; #dimension
# MAXLEAVESEST=100; #maximum #leaves in the estimator/MRP
# MAXCHECK=10; #the number to be divided (temporary)
# FLEXI=1; #0 if using maxcheck, 1 if using sequence
# SEQUENCE=( 1 );
# STARTHIST=1;
# ENDHIST=21;
# PADDING=1;
# ITERS=5;
# #echo ${SEQUENCE[@]}

# #N=10000; #sample size
# #CRITLEAVES=10; #the PQ will stop when there are CRITLEAVES 
# #MAXCHECK=2; #this is for checking when to stop splitting based on delta
# #if there are no changes in the delta max after MAXCHECK checks, stop splitting
# #note that this has been switched off     
# #to switch on, go to FinMixScheffe and change boolean stopCrit to true

# #MYDIR=RepMappedGaussian${D}D
# #mkdir $MYDIR 
# #mkdir ${MYDIR}/n${N}L${CRITLEAVES}

# for I in 10000
# do
# 	for L in ${ENDHIST}
# 	do 
# 	    #mkdir ${MYDIR}/n${I}L${L}
# 		for S in `seq 1 $TOTAL`
# 		do 
# 		echo Simulation $TOTAL for n = $I and L = $L
# 		./MappedScheffe $D $MAXLEAVESEST $I $S $L $MAXCHECK $FLEXI $STARTHIST $PADDING $ITERS ${SEQUENCE[@]} #> output.txt
# 		#mv *.txt ${MYDIR}/n${I}L${L}
# 		#append the results for each loop to an overall file #make sure in the correct folder depending on n and L
# 		cat results$S.txt >> n${I}L${L}final.txt
# 		rm results$S.txt
# 		done
# 	done
# done
