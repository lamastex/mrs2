build:
	@docker build -t mrs2 .

run:
	@docker run \
	  -it \
	  --name mrs2 mrs2 /bin/bash
### to mount a volume from host file system
	  #-v /home/raazesh/all/git/mrs2:/mrs2 \

clean:
	@docker rm mrs2 > /dev/null || true

stop:
	@docker stop mrs2

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
