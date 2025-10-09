%% Readme
%
% Decentralized Small Gain and Phase Stability Conditions for Power Systems: Limitations and Improvement
% Validation Scripts
%
% Run this script to import the system specified in the excel file
% "UserData.xlsm". It contains network specification and device definition.
% Note: To avoid errors, the file name must be exactly "UserData". When
% testing with different system, rename the file as UserData.
%
% This script generates the linearize model for each device and the network
% admittance, which will be used in the other scripts.
%
% After have run this script, run dec_conditions_inf_bus when dealing with
% the infite bus test cases, or run dec_conditions_multi_bus for a generic
% multi-bus system (IEEE14 system is consider here and in the paper as test
% case).
%
% Diego Cifelli 19/09/2025

%% Clear matlab
clear all; clc; close all; 

%% User data
UserDataName = 'UserData'; 

%% Change the current folder of matlab
cd(fileparts(mfilename('fullpath')));

%% Set user data type
% If user data is in excel format, please set 1. If it is in json format,
% please set 0.
UserDataType = 1;

%% Run toolbox
SimplusGT.Toolbox.Main();  