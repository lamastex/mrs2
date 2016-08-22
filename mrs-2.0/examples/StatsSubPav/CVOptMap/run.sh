#!/bin/bash
# runs in mathxeon3 (same machine that ran kde and posterior mean MCMC estimates in HarlowThesis2012) for OptMAP RP estimate
./CVOptMAP s 0 0. 5.0 5 2000 2 10 0.0001 1 100 2002 &> log2_2002 &&
./CVOptMAP s 0 0. 5.0 5 10000 2 10 0.0001 1 100 10002 &> log2_10002 &&
./CVOptMAP s 0 0. 5.0 5 50000 2 10 0.0001 1 100 50002 &> log2_50002 &&
./CVOptMAP s 0 0. 5.0 5 100000 2 10 0.0001 1 100 100002 &> log2_100002 &&
./CVOptMAP s 0 0. 5.0 5 2000 3 10 0.0001 1 100 2003 &> log3_2003 &&
./CVOptMAP s 0 0. 5.0 5 10000 3 10 0.0001 1 100 10003 &> log3_10003 &&
./CVOptMAP s 0 0. 5.0 5 50000 3 10 0.0001 1 100 50003 &> log3_50003 &&
./CVOptMAP s 0 0. 5.0 5 100000 3 10 0.0001 1 100 100003 &> log3_100003 &&
./CVOptMAP s 0 0. 5.0 5 2000 4 10 0.0001 1 100 2004 &> log4_2004 &&
./CVOptMAP s 0 0. 5.0 5 10000 4 10 0.0001 1 100 10004 &> log4_10004 &&
./CVOptMAP s 0 0. 5.0 5 50000 4 10 0.0001 1 100 50004 &> log4_50004 &&
./CVOptMAP s 0 0. 5.0 5 100000 4 10 0.0001 1 100 100004 &> log4_100004 &&
./CVOptMAP s 0 0. 5.0 5 2000 5 10 0.0001 1 100 2005 &> log5_2005 &&
./CVOptMAP s 0 0. 5.0 5 10000 5 10 0.0001 1 100 10005 &> log5_10005 &&
./CVOptMAP s 0 0. 5.0 5 50000 5 10 0.0001 1 100 50005 &> log5_50005 &&
./CVOptMAP s 0 0. 5.0 5 100000 5 10 0.0001 1 100 100005 &> log5_100005 &&
touch finished

