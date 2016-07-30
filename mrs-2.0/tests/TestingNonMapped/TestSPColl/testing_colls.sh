#!/bin/sh
echo "Testing collator nodes directly"
echo "default constructors and operations on default constructed nodes: no differences expected (all should be empty)"
diff -s defaultConstructor.txt copyOfdefaultConstructor.txt
diff -s defaultConstructor.txt assignmentCopyOfdefaultConstructor.txt
diff -s defaultConstructor.txt addPavingDefaultAndNULL.txt
diff -s defaultConstructor.txt addPavingDefaultAndDefault.txt
diff -s defaultConstructor.txt additionDefaultAndDefault.txt
echo
echo "constructors with box and operations on nodes with unsplit box: no differences expected"
diff -s constructorWithBox.txt copyConstructorOfConstructorWithBox.txt
diff -s constructorWithBox.txt copyAssignmentOfConstructorWithBox.txt
echo
echo "constructors and splitting and operations on nodes with split boxes: no differences expected"
diff -s collNode1.txt copyConstructorCollNode1.txt
diff -s collNode1.txt copyAssignmentCollNode1.txt
diff -s collNode2.txt shouldBeSameAsCollNode2.txt
diff -s collNode3.txt shouldBeSameAsCollNode3.txt
diff -s collNode3.txt shouldBeSameAsCollNode3FromNothingCollatedPlusColl3.txt
diff -s collNode3.txt shouldBeSameAsCollNode3FromColl3PlusNothingCollated.txt
echo
echo "addition operations on nodes with split boxes: no differences expected"
diff -s collNode1.txt addPavingsCopyOfColl1AndNULL.txt
diff -s additionColl1AndColl2.txt anotherCopyColl1PlusEqualColl2.txt
echo
echo "make collator nodes via adh collator constructor on SPS node: no differences expected"
diff -s collNode1.txt sameViaSPSForCollNode1.txt
diff -s collNode2.txt sameViaSPSForCollNode2.txt
diff -s collNode3.txt sameViaSPSForCollNode3.txt
diff -s collatorFromSPSWithNoData.txt zeroDataNode.txt
echo
echo "addition of SPS nodes: no differences expected"
diff -s collNode1.txt addPavingDefaultAndColl1ViaSPS.txt
diff -s collNode1.txt addPavingCopyColl1ViaSPSAndNULL.txt
diff -s collNode2.txt additionDefaultAndColl2ViaSPS.txt
diff -s collNode3.txt additionDefaultAndColl3ViaSPS.txt
echo
echo "more addition: no differences expected"
diff -s additionCollNode1CollNode2CollNode3.txt anotherCopyColl1PlusEqualColl2PlusEqualSPSHist3.txt
diff -s additionCollNode1CollNode2CollNode3.txt defaultPlusEqualCollNode1CollNode2CollNode3.txt
diff -s additionCollNode1CollNode2CollNode3.txt copySumOfThreePlusNothingCollatedNode.txtecho
echo
echo "adding node with box but zero data: expect to show differences with 4th zero values collated"
diff -s additionCollNode1CollNode2CollNode3.txt copySumOfThreePlusZeroData.txt
echo
echo "reduce collations: expect no differences"
diff -s reduceCopySumOfThree.txt reduceCopySumOfThreePlusNoData.txt
echo
echo "comparing averages and normalising: expect no differences"
diff -s averageFromSumOfThree.txt normalisationFromSumOfThree.txt
diff -s averageFromSumOfThree.txt averagedCopySumOfThreePlusNothingCollatedNode.txt
diff -s normalisationFromSumOfThree.txt normalisedReduceCopySumOfThree.txt
diff -s normalisationFromSumOfThree.txt normalisedCopySumOfThreePlusZeroData.txt
echo
echo "comparing normalisation of reduced collation and averaging of reduced collation: expect same leaf boxes but differences in values"
diff -s normalisedReduceCopySumOfThree.txt averagedReduceCopySumOfThree.txt
diff -s averagedCopySumOfThreePlusZeroData.txt normalisedCopySumOfThreePlusZeroData.txt
echo
echo "multiplication: expect no differences"
diff -s nothingCollatedNode.txt nothingCollatedNodeMultByTen.txt
diff -s collNode1.txt coll1MultByOne.txt
diff -s reduceSumOfThreeMultPlusAddSumOfThree.txt reduceSumOfThreeMultByTwo.txt
echo
echo "multiplication: expect difference of mult by 2 compared to += self"
diff -s sumOfThreeMultPlusAddSumOfThree.txt sumOfThreeMultByTwo.txt
echo
echo "division: expect no differences"
diff -s nothingCollatedNode.txt nothingCollatedNodeDivByTen.txt
diff -s collNode1.txt coll1DivByOne.txt
diff -s collNode1.txt coll1MultByTenDividedByTen.txt
diff -s additionCollNode1CollNode2CollNode3.txt sumOfThreeMultByTwoDividedByTwo.txt
diff -s averageFromSumOfThree.txt reducedSumOfThreeDividedByThree.txt
echo
echo "averaging on children and making copy of child, dividing by number summarised, and reducing:  expect no differences"
diff -s  coll1RightChildMakeAverage.txt coll1RightChildHandRolledAverage.txt
echo
echo "averaging and normalising on children: expect differences because averaging divides by number collated, normalising by 'area' of child node"
diff -s  coll1RightChildMakeAverage.txt coll1RightChildMakeNormalised.txt
echo
echo "averaging on children compared to copy of child of root normalised: expect no differences"
diff -s  coll1RightChildMakeAverage.txt coll1NormalisedRightChild.txt
echo
echo "collnode1 and collnode1 marginalised on all dimensions: expect no differences"
diff -s collNode1.txt collnode1MarginalisedAllReq.txt
echo 
echo "sumOfThree marginalised on 2 and sum of three marginals on 2: expect no differences"
diff -s sumOfThreeMarginalisedOn2.txt sumOfThreeMarginalsOn2.txt
echo
echo "End of testing"
echo
