%% Script to calculate and plot faults and project faults to depth, gridded according to the specified grid size
% This version also plots the shear zones beneath the brittle portion.
% Written by Zoe Mildon, 2016

% ASSUMPTIONS:
%   - the slip vector is preserved down dip
%   - the trace at the surface continues to depth
%   - the dip of the faults are consistent with depth (ie not listric
%       geometry)

% INPUTS:
% requires five individual text files with the following information: 
%       -list of fault names
%       -fault names with relevant dip
%       -fault names with relevant rake
%       -fault names with relevant projection direction
%       -for faults shorter than the depth of the seismogenic zone, a file
%       with the names and lengths (usually rounded to nearest kilometer)

% Requires the kml files of each fault to be saved individually in a single directory

% OUTPUTS:
% Writes a .inr file which can be used in Coulomb. HOWEVER BEFORE using
% in Coulomb, the #fixed value needs to be changed.

clear
format short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%      DATA TO CHANGE  %%%%%%%%%

grid_size=1000;      % Units in metres

filename='Montereale'      % Name of the outputfile that will be created

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% SETTING THE BRITTLE/DUCTILE DEPTHS %%
depth=-15000; % Depth of the seismogenic zone in metres
depth_ductile=-24000; % Depth of the ductile zone

%%%%%%%     DEFINING THE LIMITS OF THE GRID
UTM_zone=33; % Specify the UTM zone and letter
UTM_letter='T';

min_x=285; 
max_x=475;
min_y=4570;
max_y=4810; 
COUL_GRID_SIZE=10;

%%%%%%%     LOCATION OF DATA FILES
% Format of data files must be "fault_name  data"
DIPS_FILE='Data/all_faults_dips_1dp.txt'; 
%constant_dip=35;
%shear_zone_dip=65; % use this option for making all shear zones have the
%same dip, irrespective of the dip of the brittle fault
RAKES_FILE='Data/fault_rakes.txt';
PROJECTION_DIRECTION_FILE='Data/fault_projection_directions.txt';
SHORT_FAULT_LENGTHS_FILE='Data/short_fault_lengths.txt';
FAULT_DATA_PATH='Data/Field_data/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    DO NOT CHANGE ANY OTHER FILES FOR NORMAL OPERATION   %%%%%%%%%%%
model_3D_shear_zones


