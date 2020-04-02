%!/usr/bin/env octave
% =====================================================================
% NAME
% 
%    covid-19-g20-si2rd-projection
%
% DESCRIPTION
%
%    covid-19-g20-si2rd-projection computes a projection for covid-19 
%
% INPUT ARGUMENTS
%
% OUTPUT
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
%    Wednesday, April 1st, 2020
% ---------------------------------------------------------------------

% Clean up workspace. 
close all;
clear all;

% Enter the number of days in the SIRD timeseries dataset.
t_data = linspace(1,71,71);

% Read in the SIRD data from a csv file. The format of the csv file is
% assumed to be as follows:
%
%    Row 1: Date 
%    Row 2: Number of individuals who have been infected (I) with the disease. 
%    Row 3: Number of individuals who have recovred (R) from the the disease
%    Row 4: Number of individuals who have died (D) from the disease. 
%
% Note: The date row is not read in.
x_data = csvread('../../data/covid-19/jhu/covid-19-g20-infections-recoveries-deaths-jhu-20200401-33640a5.csv',1,0);

% Transpose data from row-major to column-major format.
x_data = x_data';

% Compute the active number of infectious (I) individuals in the population.
x_data(:,1) = x_data(:,1) - x_data(:,2) - x_data(:,3);

% Add a new leading column to the data matrix.
x_data = [zeros(size(x_data, 1), 1) x_data];

% Compute the number of individuals who are still suseptible (S) to the
% disease and save the data into the new leading column. 
% 44938712+25664300+210147125+37894799+67022000+83149300+1352642280+267670543+60317546+126150000+126577691+146745098+34218169+58775022+51709098+83154997+67886004+328239523
x_data(:,1) = 3172902207 - x_data(:,2) - x_data(:,3) - x_data(:,4);

# Remove initial days of epidemic, if desired.
%t_c=42;
%t_data(1:t_c)=[];
%x_data(1:t_c,:)=[];

% Define initial conditions for the population, where x0 = [s0; iU; iC, r0; d0].
iU = 120000;
x0 = [3172902207-iU; iU; 7; 0; 0]; % t_c=0

% Estimate the value for each of the key parameters that will dictate
% the expected progression of the disease in the population, including 
% the infectivity rate (beta), the recovery rate (gamma), and the death
% rate (mu), where p = [beta nu gamma_u gamma_c mu].
p = [1/3.4 0.00009 1.00*(1/9) 0.99*(1/12) 0.01*(1/3)];

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
axis ([1 300 0 1000000]);
text(225, 600000, strcat('beta=' , mat2str(p(1),3)));
text(225, 550000, strcat('nu=', mat2str(p(2),3)));
text(225, 500000, strcat('gamma_u=', mat2str(p(3),3)));
text(225, 450000, strcat('gamma_c=', mat2str(p(4),3)));
text(225, 400000, strcat('mu=', mat2str(p(5),3)));
text(225, 350000, strcat('err=', num2str(err,3)));
print -dpng covid-19-g20-si2rd-projection-jhu-20200330-ce9872c.png


% =====================================================================