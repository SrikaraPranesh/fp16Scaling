% This script generates all the tables in the paper 
% "Squeezing a Matrix Into Half Precision,
% with an Application to Solving Linear Systems -- N.J Higham,
% S. Pranesh, and M. Zounon"


% On the laptop which was used for numerical experiments
% the execution of the code takes around 90 minutes.

clear all; close all
addpath('matrices_for_testing')
addpath('GMRES_IR')
addpath('DiagScaling')
addpath('Scaling_Algorithms')
addpath('Misc_Scripts')
addpath('main_scripts')
addpath('all_suitable_SuitSparse_Matrices')


% Main_SimpleScaling;
Main_DiagonalScaling;
condition_number;
Underflow_Percentage


movefile('*.txt','results')


