//Gaurav Anil Yeole
//UFID: 54473949
//EEL6935 - Distributed Computing, Spring 2017
//Homework - 3A

#include <memory.h> 
#include "timeServer.h"


static struct timeval TIMEOUT = { 25, 0 };

long *
timeserver_1(void *argp, CLIENT *clnt)
{
	static long clnt_res;

	memset((char *)&clnt_res, 0, sizeof(clnt_res));
	if (clnt_call (clnt, TIMESERVER,
		(xdrproc_t) xdr_void, (caddr_t) argp,
		(xdrproc_t) xdr_long, (caddr_t) &clnt_res,
		TIMEOUT) != RPC_SUCCESS) {
		return (NULL);
	}
	return (&clnt_res);
}
