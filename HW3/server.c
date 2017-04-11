//Gaurav Anil Yeole
//UFID: 54473949
//EEL6935 - Distributed Computing, Spring 2017
//Homework - 3A

#include "timeServer.h"
#include<time.h>

#define SPLIT_S_ADDR_INTO_BYTES( \
    s_addr) \
    ((s_addr) >> 24) & 0xFF, \
    ((s_addr) >> 16) & 0xFF, \
    ((s_addr) >>  8) & 0xFF, \
    ((s_addr)      ) & 0xFF


long *
timeserver_1_svc(void *argp, struct svc_req *rqstp)
{
	static long  result;

	static time_t rawtime;
	
	
	time(&rawtime);
	
	printf( "IP Address of Client is: %hu.%hu.%hu.%hu\t",  SPLIT_S_ADDR_INTO_BYTES(ntohl(rqstp->rq_xprt->xp_raddr.sin_addr.s_addr)));
	printf("Port is %d\n",ntohs(rqstp->rq_xprt->xp_raddr.sin_port));
	printf("-------------------------------------------------------\n");
	
	return &rawtime;
}
