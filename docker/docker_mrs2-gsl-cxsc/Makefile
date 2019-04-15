build-gsl:
	@docker build -t lamastex/mrs2-gsl .

run-gsl:
	@docker run \
	  -it \
	  --name mrs2-gsl mrs2-gsl

# because cxsc needs interactive build we will call 
# 'make run-gsl-cxsc' and 'docker commit mrs2-gsl-cxsc lamastex/mrs2:latest' 
run-gsl-cxsc:
	@docker run \
	  -it \
	  --name mrs2-gsl-cxsc lamastex/mrs2-gsl /bin/bash
### to mount a volume from host file system
	  #-v /home/raazesh/all/git/mrs2:/mrs2 \

clean:
	@docker rm mrs2-gsl > /dev/null || true

stop-all:
	@docker stop mrs2-gsl
	@docker stop mrs2-gsl-cxsc

start:
	@docker start mrs2
####################################################################################
## others as needed
####################################################################################
#mount:
#	sudo mount -t iso9660 -o loop,ro -v $(shell pwd)/dist/centos.iso dist/centos
#
#tty:
#	@docker exec -it mrs2 /bin/bash
#import:
#	@docker exec -it mrs2 mrs2 import --name=centos7 --arch=x86_64 --path=/mnt
#
#all: build mount run import
#
#vbox:
#	vboxmanage hostonlyif ipconfig vboxnet0 --ip 192.168.56.1
####################################################################################
