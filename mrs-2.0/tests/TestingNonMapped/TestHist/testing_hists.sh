#!/bin/sh
echo "Testing histogram basics"
echo
echo "Constructors and assignment: expect no differences"
diff -s defaultConstructor.txt  copyOfDefaultConstructor.txt
diff -s defaultConstructor.txt copyAssignmentOfDefaultConstructor.txt
echo
diff -s Hist1basic.txt copyConstructHist1basic.txt
diff -s Hist4basic.txt copyConstructHist4basic.txt
diff -s Hist4basic.txt copyConstructHist4ResetHoldAllStatsTrue.txt
diff -s Hist1basic.txt copyAssignHist1basic.txt
diff -s Hist2basic.txt copyAssignHist2basic.txt
diff -s Hist4basic.txt copyAssignmentHist4ResetHoldAllStatsTrue.txt
echo
echo "Testing clearing hists: expect differences only in counts (zero in cleared copies)"
diff -s Hist1basic.txt clearedCopyHist1.txt
diff -s Hist4basic.txt clearedCopyHist4.txt
echo
echo "Testing reinsertion of same data to cleared hists: expect no differences"
diff -s Hist1basic.txt clearedCopyHist1WithDataReinserted.txt
diff -s Hist4basic.txt clearedCopyHist4WithDataReinserted.txt
echo
echo "Testing resetting holdAllStats: expect no differences"
diff -s Hist1basic.txt holdAllStatsResetToFalseForCopyHist1.txt
diff -s Hist1basic.txt holdAllStatsResetToTrueForCopyHist1.txt
diff -s Hist4basic.txt holdAllStatsResetToTrueForCopyHist4.txt
diff -s Hist4basic.txt holdAllStatsResetToFalseForCopyHist4.txt
echo
echo "Testing swap: expect no differences"
diff -s Hist1basic.txt shouldNowBeCopyOfHist1Swap1.txt
diff -s Hist2basic.txt shouldNowBeCopyOfHist2Swap1.txt
diff -s Hist2basic.txt shouldNowBeCopyOfHist2Swap2.txt
diff -s Hist4basic.txt shouldNowBeCopyOfHist4Swap2.txt
echo
echo "Checking extra details: only differences should be in the addresses of this and data collection"
diff -s copyOfHist1detailsSwap1.txt shouldNowBeCopyOfHist1detailsSwap1.txt
diff -s copyOfHist2detailsSwap1.txt shouldNowBeCopyOfHist2detailsSwap1.txt
diff -s copyOfHist2detailsSwap2.txt shouldNowBeCopyOfHist2detailsSwap2.txt
diff -s copyOfHist4detailsSwap2.txt shouldNowBeCopyOfHist4detailsSwap2.txt
echo
echo
echo "Testing histogram arithmetic"
echo
echo "Addition and no paving or no box histograms: expect no differences"
diff -s Hist1.txt AdditionNoPavingAndHist1.txt
diff -s Hist1.txt AdditionHist1AndNoPaving.txt
diff -s Hist4.txt AdditionNoPavingAndHist4.txt
diff -s Hist4.txt AdditionHist4AndNoPaving.txt
diff -s Hist4.txt AdditionNoPavingHoldAllStatsTrueAndHist4.txt
diff -s Hist4.txt AdditionHist4AndNoPavingHoldAllStatsTrue.txt
diff -s Hist1.txt AdditionHist1AndEmptyBox.txt
diff -s Hist1.txt AdditionEmptyBoxAndHist1.txt
diff -s Hist4.txt AdditionHist4AndEmptyBox.txt
diff -s Hist4.txt AdditionEmptyBoxAndHist4.txt
diff -s Hist4.txt AdditionEmptyBoxAllStatsTrueAndHist4.txt
diff -s Hist4.txt AdditionHist4AndNoPavingHoldAllStatsTrue.txt
echo
echo "PLus equals and no paving or no box histograms: expect no differences"
diff -s Hist1.txt NoPavingPlusEqualsHist1.txt
diff -s Hist1.txt Hist1PlusEqualsHistNoPaving.txt
diff -s Hist4.txt NoPavingPlusEqualsHist4.txt
diff -s Hist4.txt Hist4PlusEqualsNoPaving.txt
diff -s Hist4.txt Hist4PlusEqualsNoPavingHoldAllStatsTrue.txt
diff -s Hist4.txt NoPavingHoldAllStatsTruePlusEqualsHist4.txt
diff -s Hist1.txt EmptyBoxPlusEqualsHist1.txt
diff -s Hist1.txt Hist1PlusEqualsHistEmptyBox.txt
diff -s Hist4.txt EmptyBoxPlusEqualsHist4.txt
diff -s Hist4.txt Hist4PlusEqualsEmptyBox.txt
diff -s Hist4.txt Hist4PlusEqualsEmptyBoxHoldAllStatsTrue.txt
diff -s Hist4.txt EmptyBoxHoldAllStatsTruePlusEqualsHist4.txt
echo
echo "Testing plus equal: expect no differences"
diff -s CopyHist1PlusEqualHist2.txt CopyHist2PlusEqualHist1.txt
diff -s CopyHist1PlusEqualHist2.txt Hist2PlusEqualsHist1ResetHoldAllStatsTrue.txt
diff -s CopyHist1PlusEqualHist2.txt DefaultPlusEqualsCopyHist1PlusEqualHist2.txt
diff -s CopyHist4PlusEqualHist1.txt CopyHist1PlusEqualHist4.txt
diff -s Hist1PlusEqualHist2PlusEqualHist3.txt Hist3PlusEqualHist2PlusEqualHist1.txt
diff -s Hist4PlusEqualHist2PlusEqualHist3.txt Hist4PlusEqualHist2PlusEqualHist3HoldAllStatsTrue.txt
echo
echo "Testing addition: expect no differences"
diff -s CopyHist1PlusEqualHist2.txt AdditionHist1AndHist2.txt
diff -s AdditionHist1AndHist2.txt AdditionHist2AndHist1.txt
diff -s AdditionHist2AndHist4.txt AdditionHist4AndHist2.txt
diff -s AdditionHist2AndHist4.txt AdditionHist2AndHist4ResetHoldAllStatsTrue.txt
diff -s Hist1PlusEqualHist2PlusEqualHist3.txt AdditionHist3AndHist2AndHist1.txt 
diff -s AdditionHist3AndHist2AndHist1.txt AdditionHist1HoldAllStatsFalseAndHist2AndHist3.txt
diff -s Hist4PlusEqualHist2PlusEqualHist3.txt AdditionHist3AndHist2AndHist4.txt 
diff -s AdditionHist3AndHist2AndHist4.txt AdditionHist3AndHist2AndHist4ResetHoldAllStatsTrue.txt
echo

echo
echo "End of testing"
