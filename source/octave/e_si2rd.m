#!/usr/bin/env octave
# =====================================================================
# NAME
# 
#    e_si2rd
#
# DESCRIPTION
#
#    e_si2rd is a real, scalar function that computes ...
#
# INPUT ARGUMENTS
#
#    x is a 5-dimensional, real-valued vector
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
#    Sunday, March 29th, 2020
# ---------------------------------------------------------------------
function error = e_si2rd(t, x0, x_data, p)
  g_si2rd = @(x,t) f_si2rd(t, x, p);
  lsode_options ( 'absolute tolerance', 0.01 );
  lsode_options ( 'relative tolerance', 0.01 );
  x = lsode(g_si2rd, x0, t);
  error = sum((x(:,3)-x_data(:,2)).^2+(x(:,5)-x_data(:,4)).^2);
endfunction
# =====================================================================