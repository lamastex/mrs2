FROM raazesh/mrs2
#Define the environment build again if needed from docker_mrs2-gsl-cxsc
#2.1 building from source - is broken as the latest gcc compiler is incompatible with cxsc 
#FROM lamastex/mrs2-gsl-cxsc

# Define environment variables
#ENV SHARE /mrs2
#ENV MRS_HOME $SHARE/mrs-2
#ENV CXSCDIR $SHARE/companions/cxsc-2-5-4
#ENV GSLDIR $SHARE/companions/gsl-2.1
#ENV LD_LIBRARY_PATH $CXSCDIR/lib:$LD_LIBRARY_PATH
#ENV LD_LIBRARY_PATH $GSLDIR/lib:$LD_LIBRARY_PATH

RUN apt-get -y update
RUN apt-get install -y vim git

#RUN wget https://github.com/lamastex/mrs2/archive/master.zip 
#RUN unzip master.zip && mv mrs2-master mrs2

# we remove the older mrs-2.o in raazesh/mrs2 and simply gitclone the latest and 
# compile there but using the companions needed for gsl and cxsc in raazesh/mrs2 image
RUN rm -rf /mrs2/mrs-2.0 && mkdir -p /git &&  echo "DO THIS NOW:" && echo "cd /git" && echo "git clone https://github.com/lamastex/mrs2.git" && echo "you can compile from docker at this https cloned git repos by typing" && echo "cd /git/mrs2/mrs-2.0/" && echo "./bootstrap" && echo "./custom_config.sh" && echo "make" && echo "cd /git/mrs2/companions/mrs-1.0-YatracosThis/" && echo "./bootstrap" && echo "./custom_config.sh" && echo "make"

#2.2 building from source - is broken as the latest gcc compiler is incompatible with cxsc 
#RUN cd /mrs2/mrs-2.0/ && ./bootstrap && ./custom_config.sh && make && make install
#CMD ["/bin/bash"]

