% This script generates all the tables in the paper 
% "Squeezing a Matrix Into Half Precision,
% with an Application to Solving Linear Systems -- N.J Higham,
% S. Pranesh, and M. Zounon"


clear all; close all

addpath('matrices_for_testing')
addpath('GMRES_IR')
addpath('DiagScaling')
addpath('Scaling_Algorithms')
addpath('Misc_Scripts')
addpath('main_scripts')


Main_SimpleScaling;
keyboard
Main_DiagonalScaling;
condition_number;
Underflow_Percentage


movefile('*.txt','results')
movefile('*.mat','results')

