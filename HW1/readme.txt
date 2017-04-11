Read Me

File Containing Main Method: TMAN.cpp

Programming Language used: C++
Operating System: Ubuntu 16.04 LTS
Compiler Information: 	Thread model: posix
						gcc version 5.4.0 20160609 (Ubuntu 5.4.0-6ubuntu1~16.04.4) 
						
Library(API) used for plotting: Gnuplot	(Gnuplot_i.hpp header file of library is included in the folder.)

To run the program type following commands to the terminal:

	$ make
	$ ./TMAN N k D n r1,r2,r3,...rn 		(for Dynamic Ring Topology)
	$ ./TMAN N k B							(for Binary Tree Topology)
	$ ./TMAN N k C							(for Crescent Moon Topology)
	
	where, N = total Number of Peers 
		   k = number of neighbors maintained by each node/peer(partial view of peer)
		   n = number of radii
		   r1,r2,..rn = values of each radius
