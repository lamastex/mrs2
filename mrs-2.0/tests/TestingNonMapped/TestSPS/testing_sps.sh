#!/bin/sh
echo "Testing Hist1SPS.txt same as Hist1SPS_Expected.txt: no differences expected"
diff -s Hist1SPS.txt Hist1SPS_Expected.txt
echo
echo "Testing sps nodes directly: no differences expected"
diff -s Hist1SPS.txt copyConstructorSPSnodeHist1.txt
diff -s Hist1SPS.txt copyAssignmentSPSnodeHist1.txt
diff -s Hist2SPS.txt copyAssignmentSPSnodeHist2.txt
diff -s Hist4SPS.txt copyAssignmentSPSnodeHist4.txt
echo
echo "Testing sps tree node by node against expected output: no differences expected"
diff -s testSPSnodeByNodeOutput.txt testSPSnodeByNodeOutputExpected.txt 
echo
echo "Testing swap: expect no differences"
diff -s Hist1SPS.txt shouldNowBeCopyOfHist1PavingSwap1.txt
diff -s Hist2SPS.txt shouldNowBeCopyOfHist2PavingSwap1.txt
diff -s Hist4SPS.txt shouldNowBeCopyOfHist4PavingSwap2.txt
diff -s Hist2SPS.txt shouldNowBeCopyOfHist2PavingSwap2.txt
echo
echo "Testing swap with swapCheckOutput: expect differences in address of node and dataItrs"
diff -s swapCheckForHist1PavingSwap1.txt  swapCheckForShouldNowBeCopyOfHist1PavingSwap1.txt
diff -s swapCheckForHist2PavingSwap1.txt  swapCheckForShouldNowBeCopyOfHist2PavingSwap1.txt
diff -s swapCheckForHist4PavingSwap2.txt  swapCheckForShouldNowBeCopyOfHist4PavingSwap2.txt
diff -s swapCheckForHist2PavingSwap2.txt  swapCheckForShouldNowBeCopyOfHist2PavingSwap2.txt
echo
echo "Testing unionTreeStructure: differences against orginal pavings should just be that data in union is all 0"
diff -s Hist1SPS.txt unionOfHist1PavingAndEmptyPaving.txt
diff -s Hist2SPS.txt unionOfHist2PavingAndNullPavingPtr.txt
echo
echo "Testing unionTreeStructure: no differences expected"
diff -s unionOfHist1AndHist2Pavings.txt unionOfHist1AndHist2PavingsExpected.txt
diff -s unionOfHist1AndHist2Pavings.txt unionOfHist1CopySwappedAndHist2CopySwappedPavings.txt
echo
echo "Testing clearing data: differences against orginal pavings should just be that data in cleared trees is all 0"
diff -s Hist1SPS.txt clearedCopySPSnodeHist1.txt
diff -s Hist4SPS.txt clearedCopySPSnodeHist4.txt
echo
echo "Testing clearing data differences against union wih empty paving: no differences expected"
diff -s unionOfHist1PavingAndEmptyPaving.txt clearedCopySPSnodeHist1.txt
echo
echo "Testing resetting counts only: no differences expected"
diff -s Hist1SPS.txt countsOnlyResetToTrueCopySPSnodeHist1.txt
diff -s Hist1SPS.txt countsOnlyResetToFalseCopySPSnodeHist1.txt
diff -s Hist4SPS.txt countsOnlyResetToFalseCopySPSnodeHist4.txt
diff -s Hist4SPS.txt countsOnlyResetToTrueCopySPSnodeHist4.txt
echo
echo "End of testing"

