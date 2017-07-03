#############################################################################################################
Annealing-Algorithm in VLSI Application
Introduction:
There is one scene in the VLSI design that the designer is required to minimize the area cost for the chip. 
Annealing-Algorithm can be applied for this purpose. The algorithm randomly permutes those modules of a 
chip and creates a trend that a smaller area cost can be accepted easier. Thus, the algorithm can finally
get a minimized floorplan.
I have improved the program to make it possible to generate the plot within the overall
constraints which is not considered in the original algorithm.
The program will keep checking the overall ratio after every permutation and will multiply the area cost 
to a parameter when the ratio rans out of final constrains. That means current permutation has a greater 
possibility to be denied. This scheme will finally ensure that the program can generate the final plot 
within the overall constrains.
#############################################################################################################
C++ 
Author: KANGLI CHU@SYRACUSE UNIVERSITY
Reference: EasyBMP lib

This program gets input file(.txt) of modules info and overall constraints 
and output 3 type of plots: 
1. The plot of final floorplan, 
2. The plot of initial floorplan,
3. The plot of temperature vs area cost.
And also, the program will show the final normalized polish expression in the terminal.

The EasyBMP lib is applied here to general the result plot only.
