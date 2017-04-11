//Gaurav Anil Yeole
//UFID: 54473949
//EEL6935 - Distributed Computing, Spring 2017
//Homework - 3A

#ifndef _TIMESERVER_H_RPCGEN
#define _TIMESERVER_H_RPCGEN

#include <rpc/rpc.h>


#ifdef __cplusplus
extern "C" {
#endif


#define TIMESERVER_PROG 0x23451111
#define TIMESERVER_VER 1

#if defined(__STDC__) || defined(__cplusplus)
#define TIMESERVER 1
extern  long * timeserver_1(void *, CLIENT *);
extern  long * timeserver_1_svc(void *, struct svc_req *);
extern int timeserver_prog_1_freeresult (SVCXPRT *, xdrproc_t, caddr_t);

#else /* K&R C */
#define TIMESERVER 1
extern  long * timeserver_1();
extern  long * timeserver_1_svc();
extern int timeserver_prog_1_freeresult ();
#endif /* K&R C */

#ifdef __cplusplus
}
#endif

#endif 
