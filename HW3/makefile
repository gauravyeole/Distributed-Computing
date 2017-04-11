#Gaurav Anil Yeole
#EEL6935 - Distributed Computing,Spring 2017
#UFID: 54473949
#Homework - 3A

CC = gcc
CFLAGS = -g -DRPC_SVC_FG
RPCGEN_FLAG = -C

all: client server



client: client.o clnt.o 
	$(CC) -o client client.o clnt.o -lnsl

server: server.o svc.o  
	$(CC) -o server server.o svc.o -lnsl



server.o: server.c timeServer.h
	$(CC) $(CFLAGS) -c server.c

client.o: client.c timeServer.h
	$(CC) $(CFLAGS) -c client.c


svc.o: svc.c timeServer.h
	$(CC) $(CFLAGS) -c svc.c

clnt.o: clnt.c timeServer.h
	$(CC) $(CFLAGS) -c clnt.c

clean:
	rm -f client client.o server server.o clnt.o svc.o 