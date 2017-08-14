This Step 1 is already done and can be pulled from docker hub as raazesh/mrs2-gsl-cxsc.
If you want to redo it from source then do the following.

### build mrs2-gsl
```
make build-gsl
```
because cxsc needs interactive build we will call
```
make run-gsl-cxsc 
```
and 


```%sh
$ docker ps
CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS              PORTS               NAMES
e2564d53e479        mrs2-gsl            "/bin/bash"         14 minutes ago      Up 14 minutes                           mrs2-gsl-cxsc
```
### commit mrs2-gsl-cxsc, login and push to dockerhub
Now, you can commit this image **after** it has been built with c-xsc library that requires interaction during its build, as detailed further below on 'how to build mrs2-gsl-cxsc'.

```%sh
$ docker commit mrs2-gsl-cxsc raazesh/mrs2-gsl-cxsc
sha256:089db06cd01bd2ab9e09eaf7a8945945f13a99911397e0a57c814b60bcac68c1

$ docker images
REPOSITORY              TAG                 IMAGE ID            CREATED             SIZE
raazesh/mrs2-gsl-cxsc   latest              089db06cd01b        2 minutes ago       653MB

$ docker login
Login with your Docker ID to push and pull images from Docker Hub. If you don't have a Docker ID, head over to https://hub.docker.com to create one.
Username (raazesh): raazesh
Password: 
Login Succeeded

$ docker push raazesh/mrs2-gsl-cxsc
The push refers to a repository [docker.io/raazesh/mrs2-gsl-cxsc]
0d36b84a7fea: Pushed 
cc73d0632e42: Pushed 
21f825eed08a: Pushed 
f345441be878: Pushed 
ca72b11239e8: Pushed 
e12075c906e8: Pushed 
4d41adcb7936: Pushed 
f215f043863e: Mounted from library/ubuntu 
0c291dc95357: Mounted from library/ubuntu 
a9ee34f9e4e2: Mounted from library/ubuntu 
e15f8eeda399: Mounted from library/ubuntu 
040ba7b9591c: Mounted from library/ubuntu 
latest: digest: sha256:6a5895aa51b2443cc64e81c27206f85428f943c3637ca4749ca7d73cb7ac3a7e size: 2840
```


### how to build mrs2-gsl-cxsc after building mrs2-gsl

```%sh
Successfully built 50d545d6f21d
Successfully tagged mrs2-gsl:latest
raazesh@raazesh-Inspiron-15-7579:~/all/git/mrs2/docker$ make run-gsl-cxsc
docker: Error response from daemon: Conflict. The container name "/mrs2-gsl-cxsc" is already in use by container "77b09cead07a0d06ba280bc21c41f715837c6a153d70f52f6671ce39409c1662". You have to remove (or rename) that container to be able to reuse that name.
See 'docker run --help'.
Makefile:12: recipe for target 'run-gsl-cxsc' failed
make: *** [run-gsl-cxsc] Error 125
raazesh@raazesh-Inspiron-15-7579:~/all/git/mrs2/docker$ docker rm mrs2-gsl-cxsc > /dev/null || true
raazesh@raazesh-Inspiron-15-7579:~/all/git/mrs2/docker$ make run-gsl-cxsc
root@e2564d53e479:/# cd mrs2/companions/cxsc-2-5-4/
root@e2564d53e479:/mrs2/companions/cxsc-2-5-4# echo $CXSCDIR
/mrs2/companions/cxsc-2-5-4
root@e2564d53e479:/mrs2/companions/cxsc-2-5-4# ./install_cxsc

##
##  CXSC is a C++ library for eXtended Scientific Computing
##
##  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
##                          Universitaet Karlsruhe, Germany
##            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
##                          Universitaet Wuppertal, Germany
##
##  This library is free software; you can redistribute it and/or
##  modify it under the terms of the GNU Library General Public
##  License as published by the Free Software Foundation; either
##  version 2 of the License, or (at your option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
##  Library General Public License for more details.
##
##  You should have received a copy of the GNU Library General Public
##  License along with this library; if not, write to the Free
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##

Type 'yes' to accept this license offer.
Type 'no' to decline this license offer.

Do you accept the terms of either license? [no] yes

---------------------------------------------------------------------------------

checking host system type....x86_64
checking build system type...Linux

Configuring for a x86_64-linux-gnu host.

Section: Compiler
=================
Possible values:    gnu    for the GNU C/C++ Compiler
                    intel  for the INTEL C/C++ Compiler
Which C++ compiler? [gnu] 
Would you like to select a special GNU C/C++ Compiler-Version? [no] 
Would you like to use compiler optimization options? [yes] 
Which optimization level (e.g. -O1, -O2, ...)? [-O3] 

Using GNU V.5.4.0 Compiler
Would you like to generate 32 or 64 Bit Code? [64] 


Section: Rounding Operations
============================
Possible values: asm    for hardware support for IEEE 754 arithmetic implemented in Assembler
                 hard   for hardware support for IEEE 754 arithmetic
                 soft   for software emulations for directed rounded floating-point operations
Which method? [asm] 

cxscconf.h succesfully created


Section: Prefix
===============
(un-)installation prefix
e.g. /usr/local/cxsc or local home directory
Prefix [/root/cxsc] /mrs2/companions/cxsc-2-5-4

CXSC will be installed into /mrs2/companions/cxsc-2-5-4


Section: Makefile
=================
Creating "Makefile" in /mrs2/companions/cxsc-2-5-4 ... Done


Section: Makefile.example
=========================
Creating "Makefile.example" in /mrs2/companions/cxsc-2-5-4/examples ... Done


Section: Library
================
Do you want to create a dynamic or static library? [dynamic]/static/both both

Create a dynamic and static library.


Section: Automatic Code Generation
==================================
Proceed with automatic code generation? [yes]/no 


Section: Build CXSC Library
===========================
make[1]: Entering directory '/mrs2/companions/cxsc-2-5-4/src'
make[2]: Entering directory '/mrs2/companions/cxsc-2-5-4/src/rts'
gcc -I. -I../.. -I../rts -I../asm -mfpmath=sse -msse2 -fPIC -O3 -fno-strict-aliasing -c a_abs_.c
gcc -I. -I../.. -I../rts -I../asm -mfpmath=sse -msse2 -fPIC -O3 -fno-strict-aliasing -c a_add_.c

...
...
...

***********************************************************
*                      FINAL RESULT                       *
***********************************************************

Number of equalities :128
Number of subsets    :0
Number of supersets  :16
Other relations      :0

*********************************************************
*    All results are equal to the expected solutions.   *
*  It seems the libraries (C-XSC and CToolbox) are OK!  *
*********************************************************
make[1]: Leaving directory '/mrs2/companions/cxsc-2-5-4/CToolbox'

Section: CXSC Results
=====================
Library successfully created
Examples successfully created
Library successfully installed


Section: CToolbox Results
=========================
Examples successfully created
Library and programs successfully installed and tested
```
