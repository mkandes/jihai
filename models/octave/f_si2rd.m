%!/usr/bin/env octave
% =====================================================================
% NAME
% 
%    f_si2rd
%
% DESCRIPTION
%
%    f_si2rd is a real, vector-valued function that computes the system 
%    of equations that describe the time-dependent evolution of the 
%    SI2RD compartmentalized epidemiological model, a modifed version 
%    of the SIRD model that attempts to account for the rate at which
%    testing for the diesase is occuring in a population of unconfirmed
%    infections individuals.
%
% INPUT ARGUMENTS
%
%    t is time. 
%
%    x is a 5-dimensional, real-valued vector that holds the total 
%    number of susceptible (S), unconfirmed infectious (I_u), 
%    confirmed infectious (I_c), recovered (R), and dead (D) 
%    individuals at some time (t). 
%
%       x(1) = S(t)
%       x(2) = I_u(t)
%       x(3) = I_c(t)
%       x(4) = R(t)
%       x(5) = D(t)
%
%    p = [beta nu gamma_u gamma_c mu]
%
%    beta is the rate at which suseptible (S) individuals in the 
%    population become infected (I_u) with the diease. i.e., the 
%    transmission rate of the diease. 
%
%    nu is the rate at which testing for the diesase in the population
%    is detecting unconfirmed infectious (I_u) individuals, converting
%    these uncormfirmed cases to comfirmed ones (I_c). 
%
%    gamma_u and gamma_c are the rates at which infected (I_u, I_c)
%    individuals in the population recover (R) from the disease.
%
%    mu is the rate at which infected (I_c) individuals in the
%    population die (D) from the disease. i.e., the mortality rate. 
%
% OUTPUT
%
%    xdot is a 5-dimensional, real-valued vector that holds the rate
%    of change in the number of susceptible (S), unconfirmed 
%    infectious (I_u), confirmed infectious (I_c), recovered (R), and
%    dead (D) individuals between each time step. 
%
% AUTHOR
%
%    Marty Kandes, Ph.D.
%    Computational & Data Science Research Specialist
%    High-Performance Computing User Services Group
%    San Diego Supercomputer Center
%    University of California, San Diego
% 
% LAST UPDATED
%
%    Saturday, April 4th, 2020
% ---------------------------------------------------------------------
function xdot = f_si2rd(t, x, p)
  xdot(1) = -p(1)*((x(2)+x(3))/(x(1)+x(2)+x(3)+x(4)+x(5)))*x(1);
  xdot(2) = p(1)*((x(2)+x(3))/(x(1)+x(2)+x(3)+x(4)+x(5)))*x(1) - p(2)*(x(2)/(x(1)+x(2)+x(3)+x(4)+x(5)))*x(2) - p(3)*x(2);
  xdot(3) = p(2)*(x(2)/(x(1)+x(2)+x(3)+x(4)+x(5)))*x(2) - p(4)*x(3) - p(5)*x(3);
  xdot(4) = p(3)*x(2) + p(4)*x(3);
  xdot(5) = p(5)*x(3);
endfunction
% =====================================================================