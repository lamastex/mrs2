#Define the environment
FROM ubuntu

# Define environment variables
ENV SHARE /mrs2
ENV MRSDIR $SHARE/mrs-2.0
ENV CXSCDIR $SHARE/companions/cxsc-2-5-4
ENV GSLDIR $SHARE/companions/gsl-2.1
#export LD_LIBRARY_PATH=$CXSCDIR/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$GSLDIR/lib:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH $CXSCDIR/lib:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH $GSLDIR/lib:$LD_LIBRARY_PATH

RUN apt-get -y update
RUN apt-get install -y wget zip make automake gawk gcc g++ vim

# Move over JDK and Scala
#ADD . /
RUN wget https://github.com/lamastex/mrs2/archive/master.zip 
RUN unzip master.zip && mv mrs2-master mrs2

#and build from mrs-2/README.md

# cleaning any such companion library directories
#RUN rm -rf /mrs2/companions/gsl-2.1 


# building from source
RUN cd /mrs2/companions/tarballs/ && tar zxf gsl-2.1.tar.gz && mv gsl-2.1 ../ && cd ../gsl-2.1/ && ./configure --prefix=$GSLDIR && make && make install


#CMD ["/bin/bash"]
# for cxsc library we need interactive build: https://stackoverflow.com/questions/40854482/interactive-command-in-dockerfile
RUN cd /mrs2/companions/tarballs/ && tar zxf cxsc-2-5-4.tar.gz && mv cxsc-2-5-4 ../ && cd ../cxsc-2-5-4/ && echo "because cxsc needs interactive build we will call\n 'make run-gsl-cxsc' and 'docker commit mrs2-gsl-cxsc raazesh/mrs2-gsl-cxsc' (see ./README.md for other steps to push to docker hub)"

