# interval_hybrid
Programs for a paper under development on using a hybrid approach to solve interval FEA problems

The baseline optimization programs for the x-braced truss structures (papova1,5,10,20, and 100) are
hybrid1optps.m   optimization by particle swarm  all angles the same
hybrid1optga.m   optimization by genetic algorithm   all angles the same

The above programs do not need the intlab library.

hybrid1ps_fp.m   hybrid particle swarm - fixed point all angles the same for x braced truss structures
hybrid2ps_sen.m  hybrid particle swarm - sensitivity all angles the same for x braced truss structures

hybrid2plot.m   Since the problem is only a function of one angle, it can just be plotted as done in this code  other values are by sensitivity


