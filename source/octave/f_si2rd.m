#!/usr/bin/env octave
# =====================================================================
# NAME
# 
#    f_si2rd
#
# DESCRIPTION
#
#    f_si2rd is a real, vector-valued function that computes the system 
#    of equations that describe the time-dependent evolution of a 
#    SIIRD compartmentalized epidemiological model. 
#
# INPUT ARGUMENTS
#
#    x is a 5-dimensional, real-valued vector that holds the total 
#    number of susceptible (S), unconfirmed infectious (I_c), 
#    confirmed infectious (I_c), recovered (R), and dead (D) 
#    individuals at some time (t). 
#
#       x(1) = S(t)
#       x(2) = I_u(t)
#       x(3) = I_c(t)
#       x(3) = R(t)
#       x(4) = D(t)
#
#    t is time.
#
#    beta is the rate at which suseptible (S) individuals in the 
#    population become infected (I) with the diease.
#
#    nu is the rate at which 
#
#    gamma is the rate at which infected (I) individuals in the
#    population recover (R) from the disease.
#
#    mu is the rate at which infected (I) individuals in the population 
#    die (D) from the disease.
#
#    n is the total number of indiviudals in the population.
#
# OUTPUT
#
# AUTHOR
#
#    Marty Kandes, Ph.D.
#    Computational & Data Science Research Specialist
#    High-Performance Computing User Services Group
#    San Diego Supercomputer Center
#    University of California, San Diego
# 
# LAST UPDATED
#
#    Monday, March 30th, 2020
# ---------------------------------------------------------------------
function xdot = f_si2rd(t, x, p)
  xdot(1) = -p(1)*((x(2)+x(3))/(x(1)+x(2)+x(3)+x(4)+x(5)))*x(1);
  xdot(2) = p(1)*((x(2)+x(3))/(x(1)+x(2)+x(3)+x(4)+x(5)))*x(1) - p(2)*(x(2)/(x(1)+x(2)+x(3)+x(4)+x(5)))*x(2) - p(3)*x(2);
  xdot(3) = p(2)*(x(2)/(x(1)+x(2)+x(3)+x(4)+x(5)))*x(2) - p(4)*x(3) - p(5)*x(3);
  xdot(4) = p(3)*x(2) + p(4)*x(3);
  xdot(5) = p(5)*x(3);
endfunction
# =====================================================================