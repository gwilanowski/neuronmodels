C++ implementation of a simple motoneuron model presented in:
Traub, R. D.: Motorneurons of different geometry and the size principle.
Biological Cybernetics 25(3), 163-176 (1977)

The program reproduces Fig. 6(top) from this paper (file "traub.png"). However, the injected current is 35 nA, 
instead of 30 nA used in the paper. The program creates a two column "traub.dat" file, ready to be plotted. 

Program uses boost::odeint library for numerical integration

Comments: gwilanowski@ibib.waw.pl
