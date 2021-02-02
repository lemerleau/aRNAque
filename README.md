# RNAEVOL
Inverse RNA Folding using an Evolutionary Algorithm
## To compile, please use the following command and make sure ViennaRNA is installed
    
    $ g++ -std=c++11 -I../H -MD -MP -c -o rnaevol.o rnaevol.cpp
    $ g++ -std=c++11 -o rnaevol rnaevol.o -lstdc++ -lRNA -lm -ldl
    
## To run the program, please use.

    $ ./rnaevol -T 1000 --target="((....)).((....)).((....)).((....))" -n 100 -mm "o"

## To test the c-code in [test/](test/), please use:
    
    $ gcc -I../H -MD -MP -c -o test.o test.c
    $ gcc -o test test.o -lstdc++ -lRNA -lm -ldl
