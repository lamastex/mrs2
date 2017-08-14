# To Use mrs-2.0 docker container 

```%sh
$ docker run -it raazesh/mrs2 /bin/bash
root@bf9520ac1c4c:/# ./mrs2/mrs-2.0/examples/MooreRejSam/Rosenbrock/
```

To test out a routine try:

```%sh
root@bf9520ac1c4c:/# ./mrs2/mrs-2.0/examples/MooreRejSam/Rosenbrock/Rosenbrock
# n_dimensions: 2  n_boxes: 100  n_samples: 50  rng_seed = 0
in FirstBox, before getBoxREInfo. k: 0
0 [ -10.000000000000000,  10.000000000000000] [ -10.000000000000000,  10.000000000000000] RE: [   0.000000000000000,   1.000000000000000] BoxIntegral: [   0.000000000000000, 400.000000000000000]
in FirstBox, after getBoxREInfo 
0 [ -10.000000000000000,  10.000000000000000] [ -10.000000000000000,  10.000000000000000]
Umax:    1.000000000000000
f_scale:    1.000000000000000  1.000000000000000E-200
bottom of updateUmax 
in FirstBox, after updateUmax 
bottom of FirstBox. 
after FirstBox, before Refine 
Umax:    0.623344308959635
Umax:    1.000000000000000
f_scale:    1.000000000000000  1.000000000000000E-200
bottom of updateUmax 
in AdaptPartition after updateUmax2 
in updateIntegral. IL, IU:    0.000000016151201 1.323273005754416E+003
# Adaptive partitioning complete. Boxes: 100  Lower bound on Acceptance Prob.: 5.23875e-09 IL, IU: 1.61511e-08   3.083
#Using log(pi)? 0
#No. of Boxes with proposal mass function <= 1e-16 39
#No. of Boxes with proposal mass function <= 1e-10 40
#No. of Boxes with proposal mass function >= 1e-6 54
#No. of Boxes with proposal mass function >= 1e-3 44
after Refine 
output has been written to MRS_IsIt1or2CoinsRangeDomainSet.txt

before Rej..SampleMany 
n_samples: 50
after Rej..SampleMany 
rs_sample IU, N, Nrs:    3.083004958896712 50 50
RSSampleMany, integral est: 0.297014
RSSampleMany mean: 
   Number of labels or topologies = 1
label: 0  proportion:    1.000000000000000
Labelled Mean:
   0.981980795094569
   1.548461800444056

n interval function calls: 199
n real function calls: 522
# CPU Time (seconds). Partitioning: 0.006193  Sampling: 0.001824  Total: 0.008017
# CPU time (secods) per estimate: 0.00016034
```

# Steps in Making mrs2 docker container
Only do this if you want to rebuild from source the minimal dependencies (GSL and C-XSC for MRS-2.0)
## Step 1 
This Step 1 is already done and can be pulled from docker hub as raazesh/mrs2-gsl-cxsc.
If you want to redo it from source then:
```%sh
$ cd docker_mrs2-gsl-cxsc
$ more README.md
```
## Step 2
This Step 2 is also already done and can be pulled from docker hub as raazesh/mrs2
```%sh
$ cd docker_mrs2

$ make build
Sending build context to Docker daemon  17.41kB
Step 1/5 : FROM raazesh/mrs2-gsl-cxsc
...
...
...
Removing intermediate container aa9aa2d933b1
Successfully built 25df48fabe1f
Successfully tagged mrs2:latest

```

Tag and pus to hub:
```%sh
$ docker tag mrs2 raazesh/mrs2

$ docker rmi mrs2

$ docker images
REPOSITORY              TAG                 IMAGE ID            CREATED             SIZE
raazesh/mrs2            latest              25df48fabe1f        10 minutes ago      1.26GB
raazesh/mrs2-gsl-cxsc   latest              62e46ee1f0f6        24 minutes ago      706MB

$ docker push raazesh/mrs2
The push refers to a repository [docker.io/raazesh/mrs2]
```
