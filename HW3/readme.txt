// Gaurav Anil Yeole
// EEL 6935: Distributed Computing, Spring 2017
// Homework 3A

Operating System: Ubuntu 16.04LTS
Language Written: C
Compiler Version: gcc version 5.4.0 20160609 (Ubuntu 5.4.0-6ubuntu1~16.04.4)

Software Needed: Kindly install rpcbind first before running the codes, command to install rpcbind:
				$ sudo apt-get install rpcbind
				(Note: kindly make sure rpcbind is running before running the code. $sudo rpcbind)\
				
This homework contains following programs:

svc.c - This is program for RPC Server stub. This program creates UDP service. This program ensures that server is always running. Server stub receieves request from client
		and calls the service routine. 

server.c - This Program contains the routine for the service provided by server. In this case, it contains fuction timeserver_1_svc which returns time to the client

clnt.c - This is program for RPC client stub. This calls rpc function client_call, to send datapacket to the server. 

client.c - This program contains the code for procedure of client. 

To run the server enter following command:
	$./server
	
To run client enter following command:
	$./client localhost

