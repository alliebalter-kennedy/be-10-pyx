% This script performs a parameter optimzation for the He-3/Be-10
% production ratio using He-3 and Be-10 measurements on Ferrar Dolerites
% from Antarctica

clear all
close all

%% Constants

constants.P3_SLHL = 124.03;             % Production rate of He-3 in
                                        % pyroxene at sea level, high
                                        % latitude; [atoms g^-1 yr^-1];
                                        % calculated using the v3 online
                                        % calculator

constants.cronusPAccepted = 5.02E+09;   % Accepted CRONUS-P concentration; 
                                        % Blard et al. 2015

constants.l10 = 4.9975E-07;             % Be-10 decay rate; [yr^-1];
                                        % Nishizumii 2007? Check reference.

constants.rho = 2.94;                   % Density of Ferrar Dolerite; 
                                        % [g cm^-3]

constants.N3_nonCosmogenic = 6e6;       % Concentration of non-cosmogenic 
                                        % He-3 in Ferrar Dolerite; [atoms
                                        % g^-1]; Balter-Kennedy et al.,
                                        % 2020 and references therein

constants.L = 150;                      % spallation attenuation length for 
                                        % sample thickness corrections.
                                        % Should be appropriate for
                                        % Antarctica.

%% Add relevant paths and load data

addpath("data/")
addpath(['/Users/alexandrabalter/MATLAB/projects/bedrock_core/Common/' ...
    'v3_20190620'])
addpath('/Users/alexandrabalter/MATLAB/projects/bedrock_core/Common/depth_profiles')

filename = 'SULA_data.txt';    % file where sample, He-3, and 
                                        % Be-10 data are stored.

data = readtable(filename);             % load data

%% Unpack data into usable form

% put all sample information into structural array for use below.

% put Sample IDs in array
samples.ID = table2cell(data(:, 'Sample_ID'));   % Sample IDs

% put sample data in array
samples.thickness = table2array(data(:, 'thickness'));   % Sample thickness

% put He-3 data in array
samples.N3 = table2array(data(:, 'N3_LDEO'));         % He-3 concentrations
samples.dN3 = table2array(data(:, 'dN3_LDEO'));       % and uncertainty,
                                                      % measured at LDEO;
                                                      % [atoms g^-1]
                                                      
% CRONUS-P (CPX-2) He-3 concentration measured at LDEO; [atoms ^g-1]
samples.cronusPMeasured = zeros(5, 1);
samples.cronusPMeasured(:) = samples.N3(strcmp(samples.ID, 'CPX-2'));

% put Be-10 data in array
samples.N10 = table2array(data(:, 'N10'));      % Be-10 concentrations
samples.dN10 = table2array(data(:, 'dN10'));    % and uncertainty; prepared
                                                % at LDEO, measured at
                                                % LLNL-CAMS; [atoms g^-1]

% put location data in array
samples.lat = table2array(data(:, 'lat'));      % Latitude; [DD]
samples.long = table2array(data(:, 'long'));    % Longitude; [DD]
samples.elv = table2array(data(:, 'elv'));      % Elevation; [m]
samples.atm = table2cell(data(:, 'atm'));       % Amosphere model; 'ant' is
                                                % Antarctic; 'std' is
                                                % standard
%% correct data for thickness and shielding

% calculate thickness correction factor. 'thickness' function is from 
% Balco et al. (2008)

correctionFactor.thickness = [1./thickness(samples.thickness, ...
    constants.L, constants.rho)]'; 

% calculate CRONUS-P correction factor for standardizing He-3
% concentrations

correctionFactor.cronusP = constants.cronusPAccepted./samples.cronusPMeasured; 

% shielding correction
correctionFactor.shielding = table2array(data(:, 'shielding'));   % Shielding factor

% correct He-3 concentrations
samples.N3_corrected = samples.N3 .* correctionFactor.shielding .* ...
    correctionFactor.thickness .* correctionFactor.cronusP;

% correct Be-10 concentrations
samples.N10_corrected = samples.N10 .* correctionFactor.shielding .* ...
    correctionFactor.thickness;

%% find He-3 production rate at each sample location

production.constants = bedrock_constants();

% Atmospheric pressure at site
production.siteP = ERA40atm(samples.lat,samples.long ...
    ,samples.elv); % site air pressure

% Define production rate info
production.SFsp = stone2000(samples.lat,production.siteP,1); % scaling factor

% Build a data structure with production rate information
production.p.P3sp = production.constants.P3p_St .* production.SFsp; % He-3 spallation production rate at surface in pyroxene

%% Optimize for exposure time

% opts = optimset('fminsearch');
% opts = optimset(opts,'display','iter');

x0 = [35 5 5 5 5 5]; 

[optX_t, fmin_t] = fminsearch(@(X) fit_t(X, samples, constants, production), x0);

disp('He3/Be-10 production ratio optimized using constant exposure equations')
disp(['He-3/Be-10 production ratio is ' sprintf('%0.2f',optX_t(1))])
disp(['Be-10 production rate is ' sprintf('%0.2f',constants.P3_SLHL./optX_t(1))])
disp(['Reduced chi-squared is ' sprintf('%0.2f',fmin_t./(2*length(samples.ID)-length(x0)))]);

%% Optimize for steady erosion

x0 = [35 5 5 5 5 5];

[optX_ee, fmin_ee] = fminsearch(@(X) fit_ee(X, samples, constants, production), x0);

disp('---')
disp('He3/Be-10 production ratio optimized using steady erosion equations')
disp(['He-3/Be-10 production ratio is ' sprintf('%0.2f',optX_ee(1))])
disp(['Be-10 production rate is ' sprintf('%0.2f',constants.P3_SLHL./optX_ee(1))])
disp(['Reduced chi-squared is ' sprintf('%0.2f',fmin_ee./(2*length(samples.ID)-length(x0)))]);
%% Optimization using constant exposure
function [miss] = fit_t(X, samples, constants, production)

l10 = constants.l10;
Lsp = constants.L;
p3 = production.p.P3sp;

% X is [R3/10 t1...tN], with ti in Myr
t = X(2:end).*1e6;

% calculate final He-3 to Be-10 concentrations; X(1) = 3/10 P. ratio
predN10 = p3./(X(1) .* l10) .* (1-exp(-l10.*t));
predN3 = p3 .* t;

miss10 = (predN10' - samples.N10_corrected)./(samples.dN10); 
miss3 = (predN3' - samples.N3_corrected)./(samples.dN3);
miss = sum(miss10.^2 + miss3.^2);

end

%% Optimization using steady erosion
function [miss] = fit_ee(X, samples, constants, production)

% X is [R3/10 e1...eN], with ei in cm Myr

l10 = constants.l10;
Lsp = constants.L;
p3 = production.p.P3sp;

ee =  X(2:end).*1e-6.*2.94; % g cm-2 yr-1

% calculate ratio of final He-3 to Be-10 concentrations; X(1) = 3/10 P. ratio
predN10 = p3./X(1) ./ (l10 + ee./Lsp);
predN3 = p3 .* Lsp ./ ee;

miss10 = (predN10' - samples.N10_corrected)./(samples.dN10); 
miss3 = (predN3' - samples.N3_corrected)./(samples.dN3);
miss = sum(miss10.^2 + miss3.^2);

end