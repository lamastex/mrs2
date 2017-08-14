#Define the environment build again if needed from docker_mrs2-gsl-cxsc
FROM raazesh/mrs2-gsl-cxsc

# Define environment variables
#ENV SHARE /mrs2
#ENV MRS_HOME $SHARE/mrs-2
#ENV CXSCDIR $SHARE/companions/cxsc-2-5-4
#ENV GSLDIR $SHARE/companions/gsl-2.1
#ENV LD_LIBRARY_PATH $CXSCDIR/lib:$LD_LIBRARY_PATH
#ENV LD_LIBRARY_PATH $GSLDIR/lib:$LD_LIBRARY_PATH

RUN apt-get -y update
RUN apt-get install -y vim

#RUN wget https://github.com/lamastex/mrs2/archive/master.zip 
#RUN unzip master.zip && mv mrs2-master mrs2


# building from source
RUN cd /mrs2/mrs-2.0/ && ./bootstrap && ./custom_config.sh && make && make install

#CMD ["/bin/bash"]

