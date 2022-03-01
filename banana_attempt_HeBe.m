% Banana plots would have lines of constant exposure and lines of constant
% erosion. 

clear all
close all

t = [0:1e4:50e6];
erosion = [0:1e-6:50e-5];

HeBeRatio = 34;
l10 = 4.99e-7;
rho = 2.94;
Lsp = 150;

constants.mask = [1:6 8:9];             % use to exlude samples if necessary

constants.cronusPAccepted = 5.02E+09;   % Accepted CRONUS-P concentration; 
                                        % Blard et al. 2015
constants.N3_nonCosmogenic = 6e6;       % Concentration of non-cosmogenic 
                                        % He-3 in Ferrar Dolerite; [atoms
                                        % g^-1]; Balter-Kennedy et al.,
                                        % 2020 and references therein
%% Load LABCO data
% Location data

loc.latitude = -77.54976;       % [DD]
loc.longitude = 160.9578;       % [DD]
loc.elevation = 990.2;          % [m]

loc.atm = 'ant';                % Atmosphere model flag for v3
loc.shielding = 1;              % shielding
loc.yearCollected = 2010;       % year core was collected


%% Add relevant paths

addpath("data/")
addpath(['/Users/alexandrabalter/MATLAB/projects/bedrock_core/Common/' ...
    'v3_20190620'])
addpath(['/Users/alexandrabalter/MATLAB/projects/bedrock_core/Common/' ...
    'depth_profiles'])

%% Load data

% Use LDEO He-3 data only. Still deciding whether BGC data should be used
% in paper, not presented here.

filename = 'LABCO_data.txt';             % file where sample, He-3, and 
                                        % Be-10 data are stored.

data = readtable(filename);             % load data

%% Unpack data

% put Sample IDs in array
samples.ID = table2cell(data(strcmp(data.lab, 'LDEO'), 'sample_ID'));
samples.ID = samples.ID(constants.mask);

% put He data in arrays
samples.N3 = table2array(data(strcmp(data.lab, 'LDEO'), 'N3_LDEO'));
samples.N3 = samples.N3(constants.mask);

samples.dN3 = table2array(data(strcmp(data.lab, 'LDEO'), 'dN3_LDEO'));
samples.dN3 = samples.dN3(constants.mask);

samples.N4 = table2array(data(strcmp(data.lab, 'LDEO'), 'N4_LDEO'));
samples.N4 = samples.N4(constants.mask);

samples.dN4 = table2array(data(strcmp(data.lab, 'LDEO'), 'dN4_LDEO'));
samples.dN4 = samples.dN4(constants.mask);

% get info for CRONUS-P correction
samples.cronusPMeasured = table2array(data(strcmp(data.lab, 'LDEO'), 'CPX_LDEO'));
samples.cronusPMeasured = samples.cronusPMeasured(constants.mask);

% put Be-10 data in array
samples.N10 = table2array(data(strcmp(data.lab, 'LDEO'), 'N10_LDEO'));
samples.N10 = samples.N10(constants.mask);
samples.dN10 = table2array(data(strcmp(data.lab, 'LDEO'), 'dN10_LDEO'));
samples.dN10 = samples.dN10(constants.mask);

% put sample data in array
samples.avgDepth = table2array(data(strcmp(data.lab, 'LDEO'), 'depth'));
samples.avgDepth = samples.avgDepth(constants.mask);

samples.td = table2array(data(strcmp(data.lab, 'LDEO'), 'top_depth'));
samples.td = samples.td(constants.mask);
samples.bd = table2array(data(strcmp(data.lab, 'LDEO'), 'bottom_depth'));
samples.bd = samples.bd(constants.mask);
%% Standardize He-3

% calculate CRONUS-P correction factor for standardizing He-3
% concentrations

correctionFactor.cronusP = constants.cronusPAccepted./samples.cronusPMeasured; 

samples.N3_standardized = correctionFactor.cronusP .* samples.N3;

samples.N3_cosmogenic = samples.N3_standardized - constants.N3_nonCosmogenic;


%% Be-10 production profile set up

productionRates.constants = bedrock_constants();

productionRates.constants.P10q_St = (1./HeBeRatio).*productionRates.constants.P3p_St;

% Atmospheric pressure at site
productionRates.siteP = ERA40atm(loc.latitude,loc.longitude ...
    ,loc.elevation); % site air pressure

% Build and load muon profile
% build_muon_profile.m builds a production rate profile defined on a grid
% for efficient integration later after Balco, 2017. 
productionRates.m = build_muon_profile_w14c(productionRates.siteP, ...
    productionRates.constants,0);

% Define production rate info
productionRates.SFsp = stone2000(loc.latitude,productionRates.siteP,1); % scaling factor

productionRates.p.P3sp = productionRates.constants.P3p_St .* productionRates.SFsp; 
productionRates.p.P10sp = (productionRates.p.P3sp ./ HeBeRatio) - productionRates.m.P10mu(1); 

productionRates.p.Lsp = Lsp;

productionRates.p3Depth = productionRates.p.P3sp.*exp(-((samples.avgDepth.*rho)./Lsp));
productionRates.p10Depth = PofZ(samples.avgDepth.*rho, productionRates.m, productionRates.p, 10);

%% Normalize Be and He concentrations to PR

banana.HeNorm = samples.N3_cosmogenic ./ productionRates.p3Depth;
banana.BeNorm = samples.N10 ./ productionRates.p10Depth;

%% plot

banana_background(HeBeRatio, Lsp, l10, rho, erosion, t);
hold on
plot(banana.HeNorm, banana.BeNorm./banana.HeNorm, 'bo', 'MarkerFaceColor', 'b')
set(gca, 'xscale', 'linear', 'xlim', [1e6 1.5e6], 'FontSize', 14)