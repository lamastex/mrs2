# To Use mrs-2.0 docker container 
(See next heding for how to build this docker container from source.)

Due to changes in gcc compiler collection and the frozen c-xsc project, we are using docker container with the compatible compiler collection for c-xsc to do the compilation and development of mrs2 for now.

```%sh
# run the docker image, from hub if needed while mounting current dir in docker image's /git directory
$ docker run -d -it -v "`pwd`":/git lamastex/mrs2 
$ docker ps # list the docker processes runing to get the container-ID = 1b486f581c1e THIS will be different for you! 
CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS              PORTS               NAMES
1b486f581c1e        lamastex/mrs2       "/bin/bash"         34 minutes ago      Up 34 minutes                           loving_booth
```

Now, you can execute bash in the runnign container (use YOUR conteiner-ID output from above command: `docker ps` for `lamastex/mrs2` image):

```
$  docker exec -it 1b486f581c1e bash
```

**Tip**: See https://stackoverflow.com/questions/39794509/how-to-open-multiple-terminals-in-docker to be able to use multiple terminals into the docker container.

Now, go into the git cloned repository and start working.

You can follow the same steps of bootstrap, configure and make in `cd git/mrs2/companions/mrs-1.0-YatracosThis/`.

```%sh
root@1b486f581c1e:/git/mrs2/mrs-2.0# history
    1  cd git/
    2  ls
    3  git clone  https://github.com/lamastex/mrs2.git
    4  cd mrs2/mrs-2.0/
    5  ./bootstrap 
    6  ./custom_config.sh 
    7  make
    8  ./examples/MooreRejSam/Rosenbrock/Rosenbrock
```


# Steps in Making mrs2 docker container
Only do this if you want to rebuild from source the minimal dependencies (GSL and C-XSC for MRS-2.0). 
Currently the lamastex/mrs2 image is built from raazesh/mrs2 that was successfully built and tested when c-xsc2 was compatible with gcc version from 2017 (see Step 1 and 2 below).

## Step 0
This is to build, tag and push lamastex/mrs2 to docker hub.

```
$ cd docker_mrs2
$ make build
$ docker tag mrs2 lamastex/mrs2
$ docker push lamastex/mrs2
$ docker run -it lamastex/mrs2 bash
```

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

