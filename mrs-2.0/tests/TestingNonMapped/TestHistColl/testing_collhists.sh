#!/bin/sh
echo "Testing collator histograms"
echo "default constructors and operations on default constructed nodes: no differences expected (all should be empty)"
diff -s defaultConstructorCollatorHist.txt copyOfdefaultConstructor.txt
diff -s defaultConstructorCollatorHist.txt assignmentCopyOfdefaultConstructor.txt
echo
echo "adds with defaults: no differences expected (all should be empty)"
diff -s defaultConstructorCollatorHist.txt additionDefaultAndDefault.txt
diff -s defaultConstructorCollatorHist.txt addToCollationDefaultAndDefaultHist.txt
echo
echo "constructors with boxes, no data: no differences expected"
diff -s defaultConstructorWithHistogramWithBoxNoData.txt addToCollationDefaultAndNoDataHist.txt
diff -s zeroValueColl.txt defaultConstructorWithHistogramWithBoxNoData.txt
echo
echo "adds with box but no data: no differences expected"
diff -s zeroValueColl.txt addToCollationDefaultAndNoDataHist.txt
echo
echo "average with box but no data: no differences expected"
diff -s zeroValueColl.txt averageFromZeroValueColl.txt
echo
echo "adding to collations: no differences expected"
diff -s threeHistsCollator.txt threeHistsCollatorFromAddition.txt
diff -s threeHistsCollator.txt threeHistsCollatorFromPlusEqualsAddition.txt
diff -s threeHistsCollator.txt copyAssignmentThreeHists.txt
diff -s threeHistsCollator.txt copyConstructorFromThreeHists.txt
echo
echo "swap"
diff -s threeHistsCollator.txt shouldBeCopyOfThreeHists_afterSwap.txt
diff -s twoHistsCollator.txt shouldBeCopyOfTwoHists_afterSwap.txt
echo
echo "adds of defaults to real hists: no differences expected"
diff -s threeHistsCollator.txt shouldBeCopyOfThreeHists_afterPlusEqualsEmptyCollator.txt
echo
echo "adds of zero value hists to real hists: expect difference of extra 0 value in collations"
diff -s threeHistsCollator.txt threeHistsPlusZeroValueCollator.txt
echo
echo "averaging and normalising: no differences expected"
diff -s averageFromThreeHistsColl.txt normalisationFromThreeHistsColl.txt
diff -s averageFromThreeHistsColl.txt averageThreeHistsPlusNothingCollatedColl.txt
diff -s averageFromThreeHistsColl.txt averageThreeHistsPlusNothingCollatedColl.txt
diff -s collHist1.txt averageOfCollationOfFourOfCollHist1.txt
echo
echo "averaging and normalising: expect differences because average of three + zero value is not same as its normalisation"
diff -s averageThreeHistsPlusZeroValueColl.txt normalisationThreeHistsPlusZeroValueColl.txt
echo
echo "marginalising: no differences expected"
diff -s collHist1.txt makeMarginalCollHist1AllDimsReq.txt
diff -s makeMarginalOnThreeHistsOnDim2.txt collationMarginalsOfThreeHistsOnDim2.txt
diff -s makeMarginalOnAverageOfThreeHistsOnDim2.txt averageOfCollationMarginalsOfThreeHistsOnDim2.txt
diff -s makeMarginalOnAverageOfThreeHistsOnDim2.txt averageOfMakeMarginalThreeHistsOnDim2.txt
echo
echo "marginalising: differences in missing dimensions"
diff -s zeroValueColl.txt zeroValueCollMakeMarginal.txt
diff -s collHist1.txt makeMarginalCollHist1OnDim1.txt
diff -s collHist1.txt makeMarginalCollHist1OnDim2.txt
diff -s threeHistsCollator.txt makeMarginalOnThreeHistsOnDim2.txt
echo
echo "exporting and importing - expect no differences"
diff -s threeHistsExport1.txt  threeHistsImportExport1.txt
diff -s threeHistsCollator.txt threeHistsCollatorImported1.txt
diff -s threeHistsExport2.txt  threeHistsImportExport2.txt
diff -s threeHistsExport1.txt  threeHistsImportExport3.txt
echo
echo "exporting and importing - expect differences because export precision < original output precision"
diff -s threeHistsCollator.txt threeHistsCollatorImported2.txt
echo
echo "End of testing"
echo
