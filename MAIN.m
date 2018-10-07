% This script generates all the tables in the paper 
% "Squeezing a Matrix Into Half Precision,
% with an Application to Solving Linear Systems -- N.J Higham,
% S. Pranesh, and M. Zounon"


% On the mac book of S.Pranesh the code would usually
% take around 45 minutes for complete execution.

clear all; close all
addpath('matrices_for_testing')
addpath('GMRES_IR')
addpath('DiagScaling')
addpath('Scaling_Algorithms')
addpath('Misc_Scripts')
addpath('main_scripts')


Main_SimpleScaling;
Main_DiagonalScaling;
Main_Rank1Scaling;
condition_number;
DiagFailExample;


movefile('*.txt','results')


