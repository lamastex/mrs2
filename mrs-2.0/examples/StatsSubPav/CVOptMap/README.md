# Prior-selected MAP or optimally penalized ML estimate

The simulations are based on independent samples from Gaussian mixture density for comparison with Zhang's MCMC badwidth-selected KDE method.

TODO: For reliable L1-error computations beyond 5 dimensions one needs to simulate from piece-wise constant approximations to target densities.  
The code for generating samples for piecewise constant approx can be found in:

* `mrs-2.0/examples/StatsSubPav/MCMC/MCMCFunctionSimGaussian.cpp` - lines 93 to 178
* The L1 errors calculation: `getIAE()` - line 229 
* The `getIAE` function is defined in `piecewise_constant_function.hpp`

## To run

```%sh
# use appropriate paths!
CXSCDIR=/home/rsa64/all/git/mrs2/companions/cxsc-2-5-4
GSLDIR=/home/rsa64/all/git/mrs2/companions/gsl-2.1
export LD_LIBRARY_PATH=${CXSCDIR}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${GSLDIR}/lib:$LD_LIBRARY_PATH

# for data from a file do
```%sh
./CVOptMAP dataCVOptMAP/datasets/dp.txt 0 0. 5.0 3 1 2
```
# now call the program - see the argc and argv for the input details
./CVOptMAP s 0 0. 5.0 5 2000 2 2 0.0001 1 100 2002 &> logOutput_2_reps.txt
```

# To generate simulated data into a file and get the histogram:

```%sh
$ pwd
$ /git/mrs2/mrs-2.0/examples/StatsSubPav/CVOptMap
```

```%sh
$ pushd ../examples/MooreRejSam/Rosenbrock/
$ ./Rosenbrock 2 10000 100000
$ popd
$ mkdir -p dataCVOptMAP/datasets/
$ cd dataCVOptMAP/datasets/
$ cut -f 2,3 ../../../../MooreRejSam/Rosenbrock/MRS_Rosenbrock.samples > MRS_Rosenbrock.txt 
$ cd ../../
$ ./CVOptMAP dataCVOptMAP/datasets/MRS_Rosenbrock.txt 0 0. 5.0 3 1 2
```


### Summary of the run

```%sh
cat summaryAllReps_d_2_n_2000_r_2_2002.txt

n, dim, minPoints, minVolume, chooseStarts, keep, optimal_temperature, CVScores_opt, lv1outCVScore are : 
2000	2	1	0.0001	100	10.646584	0.061346	 -0.061083
-------------------------------------------------------------------
L1	KL	Timing	Leaves	
  0.340406	  0.268218	6.31921	183
  0.326602	  0.258687	5.80963	196
```

### Here are snippets of the output in `logOutput_2_reps.txt`:

There are two replicated: `rep number 1` and `rep number 2`. The output has been truncated. This program needs more TLC for industrial use. The detailed logs are mainly for mathematical intuition.


```%sh
cat logOutput_2_reps.txt

rep number 1

Generate 2000 random values:
 doing automatic prior selection by CV 
A box is being made for the data.  The box is 
[ -5.034515,  5.262999]  [ -4.790465,  5.166209]  

 checkMaxStep = 10

launching full posterior check from root :
stepNumLeaves 1 and posterior support -9.360291E+003 and posterior -9.360291E+003

Note: posterior search started at leaves 21 never found a higher max than its starting point
Note: posterior search started at leaves 61 never found a higher max than its starting point

...
...// a lot more diagnostic output!!!
...

number of max points is 101
Carved launch points
1	10	18	30	39	44	60	64	80	81	100	110	117	129	131	141	151	170	180	188	196	208	211	223	231	241	251	261	280	290	300	310	320	328	340	349	360	367	379	387	391	401	411	421	431	441	451	461	480	490	500	510	520	530	540	549	553	570	575	590	596	610	619	627	636	649	660	661	671	681	691	701	711	721	731	741	751	761	771	790	800	810	820	830	840	850	860	870	880	890	900	910	920	930	940	950	960	970	980	990	999	
Max posterior points
170	171	173	174	182	127	131	135	149	152	156	166	173	185	187	196	183	204	206	225	217	234	237	244	257	262	285	283	289	299	307	319	327	337	349	358	367	376	386	396	400	410	418	430	438	448	460	470	481	491	501	511	521	531	541	550	554	571	576	591	597	611	620	628	637	650	661	662	672	682	692	702	712	722	732	742	752	762	772	790	800	810	820	830	840	850	860	870	880	890	900	910	920	930	940	950	960	970	980	990	999	
maxPosteriors
-6.550169E+003	-6.551022E+003	-6.552036E+003	-6.552890E+003	-6.538166E+003	-6.533343E+003	-6.540320E+003	-6.539405E+003	-6.506974E+003	-6.505971E+003	-6.512948E+003	-6.518017E+003	-6.519832E+003	-6.507580E+003	-6.509448E+003	-6.495432E+003	-6.494533E+003	-6.502259E+003	-6.510046E+003	-6.516617E+003	-6.524243E+003	-6.534465E+003	-6.538803E+003	-6.543336E+003	-6.550549E+003	-6.549846E+003	-6.559463E+003	-6.557713E+003	-6.561910E+003	-6.570444E+003	-6.577634E+003	-6.587513E+003	-6.594703E+003	-6.598716E+003	-6.607571E+003	-6.609537E+003	-6.613789E+003	-6.619041E+003	-6.622192E+003	-6.626086E+003	-6.630152E+003	-6.641828E+003	-6.648527E+003	-6.658164E+003	-6.660881E+003	-6.668251E+003	-6.670324E+003	-6.677183E+003	-6.681399E+003	-6.689933E+003	-6.698468E+003	-6.707002E+003	-6.715536E+003	-6.724071E+003	-6.732605E+003	-6.740286E+003	-6.743700E+003	-6.761674E+003	-6.766635E+003	-6.780129E+003	-6.788023E+003	-6.801568E+003	-6.802501E+003	-6.807786E+003	-6.815127E+003	-6.827227E+003	-6.837279E+003	-6.837862E+003	-6.846058E+003	-6.855459E+003	-6.865663E+003	-6.873237E+003	-6.885513E+003	-6.896215E+003	-6.907572E+003	-6.918430E+003	-6.923177E+003	-6.937468E+003	-6.946196E+003	-6.960509E+003	-6.969044E+003	-6.977578E+003	-6.986112E+003	-6.994647E+003	-7.003181E+003	-7.011716E+003	-7.020250E+003	-7.028784E+003	-7.037319E+003	-7.045853E+003	-7.054388E+003	-7.062922E+003	-7.071456E+003	-7.079991E+003	-7.088525E+003	-7.097060E+003	-7.105594E+003	-7.112742E+003	-7.118504E+003	-7.129811E+003	-7.138185E+003	

maximium is the one at index 16
maxMaxLaunchPoint = 151
maxMaxPostPoint = 183
maxMaxPosterior = -6.494533E+003

At recreate point with leaves 183
from carving point 151
maxPosterior from earlier pq was -6.494533E+003 and posterior calculated here is -6.492987E+003
(logLik here is -6.211268E+003 and logPrior is -2.817186E+002)
time to get prior-selected adaptive hist = 6.31921

OptMAP estimate with 183 leaves
leave-1-out CV summand =  -0.061083
Getting the L1-distance and KL-distance: 

Number of importance sample points censored from 1000000 to 998979

n, dim, minPoints, minVolume, chooseStarts, keep are : 2000	2	1	0.0001	100	1
optimal temperature, CVScores_opt, lv1outCV, apprx L1 error, KL dist are : 0.646584	0.061346	 -0.061083	  0.340406  0.268218

rep number 2

Generate 2000 random values:
 doing automatic prior selection by CV 
A box is being made for the data.  The box is 
[ -4.778042,  4.928242]  [ -4.559362,  5.037076]  

 checkMaxStep = 10

launching full posterior check from root :
stepNumLeaves 1 and posterior support -9.168330E+003 and posterior -9.168330E+003

Note: posterior search started at leaves 61 never found a higher max than its starting point
Note: posterior search started at leaves 71 never found a higher max than its starting point

...
...// a lot more diagnostic comments!!!
...

number of max points is 101
Carved launch points
1	9	18	30	40	46	60	67	71	81	100	110	120	129	135	141	151	170	180	189	198	205	211	225	231	241	251	261	271	290	300	310	320	330	340	346	360	368	380	390	393	401	418	421	431	441	451	461	471	481	500	510	520	530	540	550	560	570	580	590	598	610	620	629	638	649	659	667	679	681	699	701	714	725	731	741	751	761	771	781	791	801	811	821	840	850	860	870	880	890	900	910	920	930	940	950	960	970	980	990	1000	
Max posterior points
301	287	293	298	183	190	192	203	206	209	210	160	232	179	244	196	194	202	207	216	230	237	243	258	264	274	284	294	325	302	310	348	355	365	376	358	395	403	392	403	406	439	428	431	443	451	463	473	506	491	500	510	520	530	540	550	560	570	580	590	598	610	620	629	638	649	659	667	679	681	699	701	714	725	731	741	751	761	771	781	791	801	811	821	857	855	876	903	885	906	938	935	946	947	959	985	984	1004	1014	995	1011	
maxPosteriors
-6.373378E+003	-6.369317E+003	-6.369121E+003	-6.368222E+003	-6.360388E+003	-6.360621E+003	-6.359055E+003	-6.355789E+003	-6.358134E+003	-6.356010E+003	-6.356437E+003	-6.360283E+003	-6.353059E+003	-6.349710E+003	-6.331885E+003	-6.321026E+003	-6.321341E+003	-6.321826E+003	-6.324501E+003	-6.327159E+003	-6.329402E+003	-6.331483E+003	-6.336394E+003	-6.339389E+003	-6.341820E+003	-6.346370E+003	-6.349111E+003	-6.347647E+003	-6.350334E+003	-6.350304E+003	-6.353458E+003	-6.357127E+003	-6.358955E+003	-6.361908E+003	-6.365850E+003	-6.364763E+003	-6.367303E+003	-6.364813E+003	-6.366514E+003	-6.369749E+003	-6.369419E+003	-6.374401E+003	-6.375898E+003	-6.377946E+003	-6.379392E+003	-6.382636E+003	-6.383955E+003	-6.389966E+003	-6.390338E+003	-6.384925E+003	-6.387164E+003	-6.390117E+003	-6.393071E+003	-6.396024E+003	-6.398978E+003	-6.401931E+003	-6.404884E+003	-6.407838E+003	-6.410791E+003	-6.413051E+003	-6.415414E+003	-6.420344E+003	-6.423298E+003	-6.423876E+003	-6.425841E+003	-6.425808E+003	-6.430174E+003	-6.434813E+003	-6.439220E+003	-6.441197E+003	-6.450558E+003	-6.450608E+003	-6.447739E+003	-6.450473E+003	-6.452288E+003	-6.457045E+003	-6.463061E+003	-6.470580E+003	-6.475208E+003	-6.483684E+003	-6.486925E+003	-6.491976E+003	-6.492841E+003	-6.498325E+003	-6.506627E+003	-6.509420E+003	-6.512125E+003	-6.515241E+003	-6.518280E+003	-6.520549E+003	-6.522982E+003	-6.526289E+003	-6.530231E+003	-6.533207E+003	-6.536006E+003	-6.538390E+003	-6.541637E+003	-6.541359E+003	-6.545939E+003	-6.550767E+003	-6.554684E+003	

maximium is the one at index 15
maxMaxLaunchPoint = 141
maxMaxPostPoint = 196
maxMaxPosterior = -6.321026E+003

At recreate point with leaves 196
from carving point 141
maxPosterior from earlier pq was -6.321026E+003 and posterior calculated here is -6.320037E+003
(logLik here is -6.126818E+003 and logPrior is -1.932197E+002)
time to get prior-selected adaptive hist = 5.80963

OptMAP estimate with 196 leaves
leave-1-out CV summand =  -0.062923
Getting the L1-distance and KL-distance: 

Number of importance sample points censored from 1000000 to 997779

n, dim, minPoints, minVolume, chooseStarts, keep are : 2000	2	1	0.0001	100	1
optimal temperature, CVScores_opt, lv1outCV, apprx L1 error, KL dist are : 1.01165	0.0629522	 -0.062923	  0.326602  0.258687
```

### Replicate runs
For replicate runs see `/outputs/run.sh`.
Also see `/outputs/runOutput2Sage.sh` for converting the L1-error and KL-distance estimate into TV-distance and Pinsker's inequality friendly tranformation of KL-distance, using sageMath (downloadable binaries at [http://www.sagemath.org/download.html](http://www.sagemath.org/download.html)) on commandline as follows:

```%sh
cat outputs/runOutput2Sage.sh

#!/bin/bash
FILE="$1"
#cat "$FILE"
# take the output of run.sh and parse it for getting mean, min, max and std in sage using numpy
cat "$FILE" | tail -10 | head -10 | sed 's/\t/,/g' | sed 's/ //g' | tr '\n' '; ' | sed '$s/;$/\n/' | sed -e "s/'/'\\\\''/g;s/\(.*\)/'\1'/" | sed '$i import numpy as np; np.set_printoptions(suppress=True); a=np.matrix(' | sed '$a ); print "----------------L1 dist and d_KL--------"; print a.mean(0); print a.min(0); print a.max(0); print a.std(0); print "----------------TV dist and sqrt(d_KL/2)--------"; b=a.copy(); b[:,0]=a[:,0]/2.0; b[:,1]=np.sqrt(a[:,1]/2.0); print b.mean(0); print b.min(0); print b.max(0); print b.std(0);\n'
# now paste the output into $sage shell
```
