%!/usr/bin/env octave
% =====================================================================
% NAME
% 
%    covid-19-g20-si2rd-projection
%
% DESCRIPTION
%
%    covid-19-g20-si2rd-projection computes a projection for the 
%    progression of covid-19 in the G20 countries (ex-China) using 
%    the SI2RD model. 
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
%    Sunday, April 5th, 2020
% ---------------------------------------------------------------------

% Clean up workspace. 
close all;
clear all;

% Enter the number of days in the SIRD timeseries dataset.
t_data = linspace(1,75,75);

% Read in the SIRD data from a csv file. The format of the csv file is
% assumed to be as follows:
%
%    Row 1: Date 
%    Row 2: Number of individuals who have been infected (I) with the disease. 
%    Row 3: Number of individuals who have recovred (R) from the the disease
%    Row 4: Number of individuals who have died (D) from the disease. 
%
% Note: The date row is not read in.
x_data = csvread('../../data/covid-19/jhu/g20/covid-19-g20-infections-recoveries-deaths-jhu-20200405-865c933.csv',1,0);

% Transpose data from row-major to column-major format.
x_data = x_data';

% Compute the active number of infectious (I) individuals in the population.
x_data(:,1) = x_data(:,1) - x_data(:,2) - x_data(:,3);

% Add a new leading column to the data matrix.
x_data = [zeros(size(x_data, 1), 1) x_data];

% Compute the number of individuals who are still suseptible (S) to the
% disease and save the data into the new leading column. 
x_data(:,1) = 3172902207 - x_data(:,2) - x_data(:,3) - x_data(:,4);

# Remove initial days of epidemic, if desired.
%t_c=42;
%t_data(1:t_c)=[];
%x_data(1:t_c,:)=[];

% Define initial conditions for the population, where x0 = [s0; iU; iC, r0; d0].
iU = 1048576;
x0 = [3172902207-iU; iU; 7; 0; 0]; % t_c=0

% Estimate the value for each of the key parameters that will dictate
% the expected progression of the disease in the population, including 
% the transmission rate (beta), the testing rate (nu), the recovery rates of 
% both unconfirmed and confirmed infections(gamma_u and gamma_c, respectively), 
% and the mortality rate (mu), where p = [beta nu gamma_u gamma_c mu].
p = [1/3.2 0.0001 1/10 1/16 0.005];

% Fit key parameters to the data. 
p = fminsearch(@(p) e_si2rd(t_data, x0, x_data, p), p);
err = e_si2rd(t_data, x0, x_data, p);

% Run the model out into the future for a projection. 
t = linspace(1,300,300);
g_si2rd = @(x,t) f_si2rd(t, x, p);
x = lsode(g_si2rd, x0, t);

% Plot the current data and the model projection. 
plot(t_data, x_data(:,2), 'or', 'markerfacecolor', 'r', 'markersize', 7, ...
     t_data, x_data(:,4), 'ok', 'markerfacecolor', 'k', 'markersize', 7, ...
     t, x(:,1), '-b', 'linewidth', 2, ...
     t, x(:,2), '-.r', 'linewidth', 2, ...
     t, x(:,3), '-r', 'linewidth', 2, ...
     t, x(:,4), '-g', 'linewidth', 2, ...
     t, x(:,5), '-k', 'linewidth', 2)
title('G20 COVID-19 SI2RD Epidemiological Model');
legend("Confirmed Infectious (data)", ...
       "Confirmed Dead (data)", ...
       "Susceptible (model)", ...
       "Unconfirmed Infectious (model)", ...
       "Confirmed Infectious (model)", ...
       "Recovered (model)", ...
       "Dead (model)", ...
       "location", ...
       "northeast");
ylabel('Number of Susceptible, Infectious, Recovered, and Dead');
xlabel('Time (days)');
axis ([1 300 0 1500000]);
text(225, 900000, strcat('beta=', mat2str(p(1),3)));
text(225, 825000, strcat('nu=', mat2str(p(2),3)));
text(225, 750000, strcat('gamma_u=', mat2str(p(3),3)));
text(225, 675000, strcat('gamma_c=', mat2str(p(4),3)));
text(225, 600000, strcat('mu=', mat2str(p(5),3)));
text(225, 525000, strcat('err=', num2str(err,3)));
print -dpng covid-19-g20-si2rd-projection-jhu-20200405-865c933.png


% =====================================================================