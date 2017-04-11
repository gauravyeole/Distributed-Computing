//Gaurav Anil Yeole
//UFID: 54473949
//EEL6935 - Distributed Computing, Spring 2017
//Homework - 3A

#include "timeServer.h"
#include<time.h>

void
timeserver_prog_1(char *host)
{
	CLIENT *clnt;
	long  *result_1;
	char *timeserver_1_arg;
	struct tm* tm_info;
	char buffer[26];

#ifndef	DEBUG
	clnt = clnt_create (host, TIMESERVER_PROG, TIMESERVER_VER, "udp");
	if (clnt == NULL) {
		clnt_pcreateerror (host);
		exit (1);
	}
#endif	/* DEBUG */

	result_1 = timeserver_1((void*)&timeserver_1_arg, clnt);
	if (result_1 == (long *) NULL) {
		clnt_perror (clnt, "call failed");
	}
	else{
		tm_info = localtime(result_1);
		
		strftime(buffer, 26, "%m/%d/%Y:%H:%M:%S", tm_info);
		
		if(tm_info->tm_hour <= 12){
			for (int i=0;i<19;i++){
			printf("%c",buffer[i]);
			}
			printf(" AM\n");
		}
		else{
			
			tm_info->tm_hour = tm_info->tm_hour - 12;
			strftime(buffer, 26, "%m/%d/%Y:%H:%M:%S", tm_info);
			for (int i=0;i<19;i++){
			printf("%c",buffer[i]);
			}
			printf(" PM\n");
		}
    	
	}
#ifndef	DEBUG
	clnt_destroy (clnt);
#endif	 /* DEBUG */
}


int
main (int argc, char *argv[])
{
	char *host;

	if (argc < 2) {
		printf ("usage: %s server_host\n", argv[0]);
		exit (1);
	}
	host = argv[1];
	timeserver_prog_1 (host);
exit (0);
}
