# interval_hybrid
Programs for a paper under development on using a hybrid approach to solve interval FEA problems

The baseline optimization programs for the x-braced truss structures (papova1,5,10,20, and 100) are
hybrid1optfmin.m   optimization by gradient   all angles the same
hybrid1optps.m   optimization by particle swarm  all angles the same
hybrid1optga.m   optimization by genetic algorithm   all angles the same

The above programs do not need the intlab library.   The fixed point methods require INTLAB which can be started using the start.m routine.   The rest of the library is available from Prof Rump

hybrid1XX_fp.m   hybrid XX optimization  - fixed point all angles the same for x braced truss structures
hybrid1ps_sen.m  hybrid XX optimation  - sensitivity all angles the same for x braced truss structures
where XX is fmin, ps or ga.

hybrid2plot.m   Since the problem is only a function of one angle, it can just be plotted as done in this code  other values are by sensitivity

hybrid3xx. files are used for the second problem set.   The angle of each load is an independent interval variable.  
