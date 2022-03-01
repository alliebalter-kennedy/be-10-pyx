% if doing erosion, need to include that in the decay correction? included
% in the theoreitcal He-3 and Be-10 depth profiles, but Be-10 conc are only
% This script performs decay/erosion corrections on He-3 and Be-10 data
% from the LABCO bedrock core, collected from Ferrar Dolerite at the
% Labyrinth, Dry Valleys, Antarctica by John Stone and colleagues.

% Inputs are sample information and He-3 and Be-10 concentrations stored in
% the .txt file 'LABCO_data.txt'. Outputs include figures and tables for the
% paper Balter-Kennedy et al., (in prep), working title, "10Be analysis in
% pyroxene - a method for routine chemical extraction"._data.txt'.

% The surface exposure age of the site is calculated using v3 of the
% Online Exposure Age Calculator described by Balco et al., (2008). Code
% for sending sample information to the calculator and retrieving exposure
% ages was modified from a script distributed at an ICE-D:GREENLAND
% workshop at the University at Buffalo, 2021. Orginal script by Greg
% Balco.

% The theoretical depth profiles are calculated using code published
% alongside Schaefer et al. (2016), original scripts by Greg Balco. 

% BGC data not used in this script. Can see "LABCO_analysis_erosion_old.m"
% file for ways to incoportate BGC data if needed. All this does
% differently is use the surface concentration from the BGC measurements,
% rather than from the data fit. 

clear all
close all


%% Constants

constants.HeBeRatio = 34;            % calculated from SULA samples,
                                        % Balter-Kennedy, in prep.

constants.P3_SLHL = 124.03;             % Production rate of He-3 in
                                        % pyroxene at sea level, high
                                        % latitude; [atoms g^-1 yr^-1];
                                        % calculated using the v3 online
                                        % calculator
loc.erosionRate = 0;              % Subaerial erosion 
                                        % rate; [cm yr^-1]. Use 0 for
                                        % Figure 1. 44.1 is steady state
                
constants.HeMuons = 'Larsen';   % How to calculate percentage of He-3 
                                        % production from muons. 
                                        % 'Larsen erosion' is from Larsen et al. 2021, with erosion; 
                                        % 'Larsen' is from Larsen 2021 w/out erosion, 
                                        % 'Be' is the same % as Be Production in qtz.

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

depthProfile.depth_cm = [0:0.1:200]';    % Depth profile for calculations;
                                        % [cm]

%% Location data

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

% put He data in arrays
samples.N3 = table2array(data(strcmp(data.lab, 'LDEO'), 'N3_LDEO'));

samples.dN3 = table2array(data(strcmp(data.lab, 'LDEO'), 'dN3_LDEO'));

samples.N4 = table2array(data(strcmp(data.lab, 'LDEO'), 'N4_LDEO'));

samples.dN4 = table2array(data(strcmp(data.lab, 'LDEO'), 'dN4_LDEO'));

% get info for CRONUS-P correction
samples.cronusPMeasured = table2array(data(strcmp(data.lab, 'LDEO'), 'CPX_LDEO'));

% put Be-10 data in array
samples.N10 = table2array(data(strcmp(data.lab, 'LDEO'), 'N10_LDEO'));
samples.dN10 = table2array(data(strcmp(data.lab, 'LDEO'), 'dN10_LDEO'));

% put sample data in array
samples.avgDepth = table2array(data(strcmp(data.lab, 'LDEO'), 'depth'));

samples.td = table2array(data(strcmp(data.lab, 'LDEO'), 'top_depth'));
samples.bd = table2array(data(strcmp(data.lab, 'LDEO'), 'bottom_depth'));

%% Standardize He-3

% calculate CRONUS-P correction factor for standardizing He-3
% concentrations

correctionFactor.cronusP = constants.cronusPAccepted./samples.cronusPMeasured; 

samples.N3_standardized = correctionFactor.cronusP .* samples.N3;

samples.N3_cosmogenic = samples.N3_standardized - constants.N3_nonCosmogenic;

samples.pctNonCosmogenic = (constants.N3_nonCosmogenic./samples.N3).*100;

%% calcualte surface age

% get surface concentration from fit to LDEO He-3 data
exponentialFitHe.coefficients = polyfit(samples.avgDepth, ...
    log(samples.N3_cosmogenic), 1);

exponentialFitHe.L3 = -1./exponentialFitHe.coefficients(1);
exponentialFitHe.surfaceConcentration = exp(exponentialFitHe.coefficients(2));

calculator.formatInput = ['LABCO_SURFACE_LDEO %.2f %.2f %.0f %s 0 %.2f %.0f %.10f %.0f; ' ...
    'LABCO_SURFACE_LDEO He-3 pyroxene %i %i NONE 0;'];
calculator.inputText = sprintf(calculator.formatInput, loc.latitude, loc.longitude, ...
    loc.elevation, loc.atm, constants.rho, loc.shielding, loc.erosionRate, ...
    loc.yearCollected, exponentialFitHe.surfaceConcentration, ...
    [0.025*exponentialFitHe.surfaceConcentration]); % typical measurement uncertainty


%% Send sample info to cosmo calculator
calculator.url = "https://hess.ess.washington.edu/cgi-bin/matweb";
calculator.xmlResult = webread(calculator.url,'mlmfile','age_input_v3', ...
    'reportType','XML', 'resultType','long','plotFlag','no','text_block', ...
    calculator.inputText);

%Load the parser and parse string from data returned from cosmo calculator
import matlab.io.xml.dom.*
calculator.xDoc = parseString(Parser, calculator.xmlResult);

%% Get the content of the stated elements to return as number

calculator.stAgeList = getElementsByTagName(calculator.xDoc,'t3pyroxene_St');
    for i = 1:length(exponentialFitHe.surfaceConcentration)
        samples.StAge(i, 1) = str2num(getTextContent(node(calculator.stAgeList,i)));
    end

calculator.stIntList = getElementsByTagName(calculator.xDoc,'delt3pyroxene_int_St');
    for i = 1:length(exponentialFitHe.surfaceConcentration)
        samples.StInternalUncertainty(i) = str2num(getTextContent(node(calculator.stIntList,i)));
    end

%% Calculate correction factors 

% calculate decay term for Be-10; includes radioactive decay + erosion
% (l + ep/L)
constants.decayTerm = constants.l10 + ((loc.erosionRate .* constants.rho)...
    ./ constants.L);

% calculate decay correction factor (should be > 1). 
correctionFactor.decay = (samples.StAge .* (constants.decayTerm)) ./ ...
    (1 - exp(-(constants.decayTerm).*samples.StAge));

if sum(correctionFactor.decay < 1) ~= 0
    print('WARNING: at least one decay correction factor >1')
end

% correct He-3 for erosion if erosion is > 0 cm yr-1

constants.erosionTerm = (loc.erosionRate .* constants.rho)...
    ./ constants.L;

correctionFactor.erosion = (samples.StAge .* (constants.erosionTerm)) ./ ...
    (1 - exp(-(constants.erosionTerm).*samples.StAge));

%% decay-correct Be-10 concentrations (includes erosion)

% perform decay correction on measured Be-10 concentrations
samples.N10_decayCorrected = samples.N10 .* correctionFactor.decay;

%% Erosion-correct He-3 concentrations if erosion > 0

if loc.erosionRate == 0
    samples.N3_erosionCorrected = samples.N3_cosmogenic;
else 

    samples.N3_erosionCorrected = samples.N3_cosmogenic .* ...
        correctionFactor.erosion;
end

%% Calculate ratios

% Ratio from measurements - He-3 is standardized to CRONUS-P
samples.measuredHeBeRatio = samples.N3_standardized ./ samples.N10;
samples.dMeasuredHeBeRatio = sqrt((samples.dN3./samples.N3).^2 + ...
    (samples.dN10./samples.N10).^2) .* samples.measuredHeBeRatio;


% Decay-corrected ratios
samples.decayCorrectedHeBeRatio = samples.N3_erosionCorrected ./...
    samples.N10_decayCorrected;
samples.dDecayCorrectedHeBeRatio = sqrt((samples.dN3./samples.N3).^2 + ...
    (samples.dN10./samples.N10).^2) .* samples.decayCorrectedHeBeRatio;

%% Plot LDEO He-3 data


% calculate He-3 concentrations at sample depths using expoential fits to
% calculate mistfit.
exponentialFitHe.N3_cosmoAtSampleDepths = polyval(exponentialFitHe.coefficients, ...
    samples.avgDepth);

% calculate mismatch between measurements and fit

exponentialFitHe.misfit = sum(((samples.N3_cosmogenic - ...
    exp(exponentialFitHe.N3_cosmoAtSampleDepths))./samples.dN3).^2);

exponentialFitHe.reducedChiSq = exponentialFitHe.misfit ...
    ./(length(samples.ID)-1);

% line for plotting
exponentialFitHe.line = polyval(exponentialFitHe.coefficients, ...
    depthProfile.depth_cm);

if loc.erosionRate > 0
    exponentialFitHe.coefficients_erosionCorrected = polyfit(samples.avgDepth, ...
    log(samples.N3_erosionCorrected), 1);

    exponentialFitHe.N3_erosionCorrectedAtSampleDepths = polyval( ...
        exponentialFitHe.coefficients_erosionCorrected, samples.avgDepth);

    exponentialFitHe.N3_erosionCorrectedAtSampleDepths = polyval( ...
        exponentialFitHe.coefficients_erosionCorrected, samples.avgDepth);
    
    exponentialFitHe.line_erosionCorrected = polyval( ...
        exponentialFitHe.coefficients_erosionCorrected, depthProfile.depth_cm);
end

figure(1)
subplot(1, 2, 1)
hold on
errorbar(samples.N3./1e8, samples.avgDepth, samples.dN3./1e8, 'horizontal', 'sk', 'MarkerFaceColor', 'k')
errorbar(samples.N3_cosmogenic./1e8, samples.avgDepth, samples.dN3./1e8, ...
    'horizontal', 'sb', 'MarkerFaceColor', 'b')
plot(exp(exponentialFitHe.line)./1e8, depthProfile.depth_cm, 'b');
plot(exp(exponentialFitHe.line(1))./1e8, depthProfile.depth_cm(1), ...
    'bp', 'MarkerFaceColor', 'b', 'MarkerSize', 10)
% if loc.erosionRate > 0
%     errorbar(samples.N3_erosionCorrected./1e8, samples.avgDepth, samples.dN3./1e8, ...
%     'horizontal', 'sg', 'MarkerFaceColor', 'g')
%     plot(exp(exponentialFitHe.line_erosionCorrected)./1e8, depthProfile.depth_cm, 'g');
% end
xlabel(['^{3}He Concentration (10^{8} atoms/g)'])
ylabel('Depth (cm)')
legend('Measured ^{3}He', 'Standardized cosmogenic ^{3}He', ...
    'Exponential Fit', 'Modeled surface concentration', 'Location', 'SE')
set(gca, 'YDir', 'reverse', 'FontSize', 14, 'ylim', [0 220], 'xscale', 'log')
grid on

%% fit to raw be-10 data

exponentialFitBe.coefficientsMeasured = polyfit(samples.avgDepth, ...
    log(samples.N10), 1);

% exponentialFitBe.L10 = -1./exponentialFitBe.coefficientsMeasured(1);
exponentialFitBe.surfaceConcentrationMeasured = exp(exponentialFitBe.coefficientsMeasured(2));

exponentialFitBe.N10_MeasuredAtSampleDepths = polyval( ...
    exponentialFitBe.coefficientsMeasured, samples.avgDepth);

% plot lines
exponentialFitBe.lineMeasured = polyval(exponentialFitBe.coefficientsMeasured, depthProfile.depth_cm);

%% fit to decay-corrected Be-10 data

exponentialFitBe.coefficientsDecayCorrected = polyfit(samples.avgDepth, ...
    log(samples.N10_decayCorrected), 1);

exponentialFitBe.L10 = -1./exponentialFitBe.coefficientsDecayCorrected(1);
exponentialFitBe.surfaceConcentration = exp(exponentialFitBe.coefficientsDecayCorrected(2));

exponentialFitBe.N10_decayCorrectedAtSampleDepths = polyval( ...
    exponentialFitBe.coefficientsDecayCorrected, samples.avgDepth);

% calculate reduced x^2 as a measure of fit

exponentialFitBe.misfit = sum(((samples.N10_decayCorrected - exp( ...
    exponentialFitBe.N10_decayCorrectedAtSampleDepths))./samples.dN10).^2, 1);

exponentialFitBe.reducedChiSq = exponentialFitBe.misfit./(length(samples.ID)-1);

% plot lines
exponentialFitBe.line = polyval(exponentialFitBe.coefficientsDecayCorrected, depthProfile.depth_cm);

subplot(1, 2, 2)
hold on
errorbar(samples.N10./1e6, samples.avgDepth, samples.dN10./1e6, ...
    'horizontal', 'sk', 'MarkerFaceColor', 'k')
% errorbar(samples.N10_decayCorrected./1e6, samples.avgDepth, ...
%     samples.dN10./1e6, 'horizontal', 'sr', 'MarkerFaceColor', 'r')
plot(exp(exponentialFitBe.lineMeasured)./1e6, depthProfile.depth_cm, 'r');
plot(exponentialFitBe.surfaceConcentrationMeasured./1e6, 0, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 10)
xlabel(['^{10}Be Concentration (10^{6} atoms/g)'])
ylabel('Depth (cm)')
legend('Measured ^{10}Be', 'Exponential Fit', 'Modeled Surface Concenration', 'Location', 'SE')
set(gca, 'YDir', 'reverse', 'FontSize', 14, 'ylim', [0 220], 'xscale', 'log')
grid on

%% Calculate ratio of fits to decay/erosion corrected values.

if loc.erosionRate == 0
    exponentialFitBe.correctedHeBeRatio = exp(exponentialFitHe.line) ./ ...
        exp(exponentialFitBe.line);
else
    exponentialFitBe.correctedHeBeRatio = exp(exponentialFitHe.line_erosionCorrected) ./ ...
        exp(exponentialFitBe.line);
end

%% Theoretical Be-10 and He-3 concentrations with depth

% Calculate theoretical Be-10 and He-3 concentrations with depth given age and erosion rate

% define depth profile
depthProfile.depth_gcm2 = depthProfile.depth_cm.*constants.rho; % [g cm^-2];

%% Be-10 production profile set up

depthProfile.constants = bedrock_constants();

% Atmospheric pressure at site
depthProfile.siteP = ERA40atm(loc.latitude,loc.longitude ...
    ,loc.elevation); % site air pressure

% Build and load muon profile
% build_muon_profile.m builds a production rate profile defined on a grid
% for efficient integration later after Balco, 2017. 
depthProfile.m = build_muon_profile_w14c(depthProfile.siteP, ...
    depthProfile.constants,0);

% Define production rate info
depthProfile.SFsp = stone2000(loc.latitude,depthProfile.siteP,1); % scaling factor

% Build a data structure with production rate information
depthProfile.p.P3tot = depthProfile.constants.P3p_St .* depthProfile.SFsp; % the calculated PR is a total PR because v3 calculator does not remove muons for 3He
depthProfile.p.P10sp = depthProfile.constants.P10q_St .* depthProfile.SFsp; % Be-10 spallation production rate at surface IN QUARTZ
depthProfile.p.P26sp = depthProfile.p.P10sp.*depthProfile.constants.R2610q; % Al-26 spallation production rate at surface
depthProfile.p.P14sp = depthProfile.constants.P14q_St.*depthProfile.SFsp; % C-14 spallation production rate at surface

% Attenuation
depthProfile.p.Lsp = constants.L; % g/cm2.

%% He-3 depth profile

% muon production rate (erosion included) from Larsen et al., 2021
depthProfile.muPctLarsenErosion = 0.23./depthProfile.constants.P3p_St; 

% muon production rate (no erosion) from Larsen et al., 2021
depthProfile.muPctLarsen = 0.45./depthProfile.constants.P3p_St; 

% based on theory in quartz, should be similar to Be-10
depthProfile.muPctBe = depthProfile.m.P10mu(1)./depthProfile.p.P10sp; 

if strcmp(constants.HeMuons, 'Larsen')
    depthProfile.muPctHe = depthProfile.muPctLarsen;
elseif strcmp(constants.HeMuons, 'Larsen erosion')
    depthProfile.muPctHe = depthProfile.muPctLarsenErosion; 
elseif strcmp(constants.HeMuons, 'Be')
    depthProfile.muPctHe = depthProfile.muPctBe;
else
    print('WARNING: indicate percent He-3 production by muons')
end

depthProfile.p.P3mu = depthProfile.muPctHe.* ... 
    depthProfile.p.P3tot; % calculate production of He-3 by muons 
depthProfile.p.P3sp = depthProfile.p.P3tot - ...
    depthProfile.p.P3mu; % calculate production of He-3 by spallation

depthProfile.p.P3Lmu = 8780; % muon attenuation length for He-3 in g/cm2, from Larsen 2021

%% Be-10 production profile calc

% re-write P10sp to take on ratio from surface samples. Since we consider
% ratio from surface samples to represent total Be-10 and He-3 production,
% we must separate out muon and spallation production. First find total P10
% at SLHL based on ratio and P3.
depthProfile.p.P10tot = (depthProfile.constants.P3p_St./constants.HeBeRatio); 

% Now subtract percentage of production by muons and re-write P10sp to take 
% on ratio from surface samples. Based on percentage of Be-10 in quartz
% produced by muons.
depthProfile.p.P10sp = (depthProfile.p.P10tot - (depthProfile.muPctBe.* ...
    depthProfile.p.P10tot)).*depthProfile.SFsp;

%% Calculate concentrations based on He-3 age

% calculate theoretical concentrations at depths

for i = 1:length(samples.ID)
    model.pfun10 = @(zz, t) PofZ_he((zz + loc.erosionRate.*constants.rho.*t), ...
        depthProfile.m,depthProfile.p,10) .* exp(-depthProfile.constants.l10.*t);
    model.pfun3 = @(zz, t) PofZ_he((zz + loc.erosionRate.*constants.rho.*t), ...
        depthProfile.m,depthProfile.p,3);

    model.N10_sampleDepths(i, 1) = integral2(model.pfun10, ...
        samples.td(i).*constants.rho, samples.bd(i).*constants.rho, ...
        0, samples.StAge,'RelTol',1e-3,'AbsTol',1e-3) ...
        ./(samples.bd(i).*constants.rho - samples.td(i).*constants.rho);
    model.N3_sampleDepths(i, 1) = integral2(model.pfun3, ...
        samples.td(i).*constants.rho, samples.bd(i).*constants.rho, ...
        0, samples.StAge,'RelTol',1e-3,'AbsTol',1e-3) ...
        ./(samples.bd(i).*constants.rho - samples.td(i).*constants.rho);
end

% calculate theoretical depth profiles

for i = 1:length(depthProfile.depth_gcm2)
    model.pfun10 = @(t) PofZ_he((depthProfile.depth_gcm2(i) + loc.erosionRate.* ...
        constants.rho.*t),depthProfile.m,depthProfile.p,10) ...
        .* exp(-depthProfile.constants.l10.*t);
    model.pfun3 = @(t) PofZ_he((depthProfile.depth_gcm2(i) + loc.erosionRate.* ...
        constants.rho.*t),depthProfile.m,depthProfile.p,3);

    model.N10_line(i, 1) = integral(model.pfun10, 0, samples.StAge, ...
        'RelTol',1e-3,'AbsTol',1e-3);
    model.N3_line(i, 1) = integral(model.pfun3, 0, samples.StAge, ...
        'RelTol',1e-3,'AbsTol',1e-3);
end
        
model.HeBeRatio = model.N3_line./model.N10_line;

%% Plot measured and theoretical ratio

figure(2)
hold on
errorbar(samples.measuredHeBeRatio, samples.avgDepth, samples.dMeasuredHeBeRatio, ...
    'horizontal', 'ks', 'MarkerFaceColor', 'k')
% plot(R_fits, profile_cm, 'k', 'Linewidth', 2) % should plot raw ratio
% fits
plot(model.HeBeRatio, depthProfile.depth_cm, 'b', 'LineWidth', 2)
set(gca, 'XLim', [20 80], 'YLim', [0 200], 'YDir', 'reverse', 'FontSize', 14)
% legend(['^{3}He/^{10}Be ratio' 10 'from measurements'], 'Linear fit to measurements', 'Theoretical ^{3}He/^{10}Be ratio', 'Location', 'NW')
ylabel('Depth (cm)')
xlabel('^{3}He/^{10}Be Ratio')
grid on

figure(3)
subplot(1, 2, 1)
hold on
errorbar(samples.N10, samples.avgDepth, samples.dN10, 'horizontal', 'ks', 'MarkerFaceColor', 'k');
plot(model.N10_sampleDepths, samples.avgDepth, 'rs')
plot(model.N10_line, depthProfile.depth_cm, 'r')
set(gca, 'YDir', 'reverse', 'FontSize', 14, 'xscale', 'log', 'ylim', [0 220])
xlabel('^{10}Be Concentration')
ylabel('Depth (cm)')
legend('Measured ^{10}Be', 'Modeled ^{10}Be', 'Modeled ^{10}Be profile', 'Location', 'Southeast')
grid on 

figure(3)
subplot(1, 2, 2)
hold on
errorbar(samples.N3_cosmogenic, samples.avgDepth, samples.dN3, 'horizontal', 'ks', 'MarkerFaceColor', 'k');
plot(model.N3_sampleDepths, samples.avgDepth, 'bs')
plot(model.N3_line, depthProfile.depth_cm, 'b')
set(gca, 'YDir', 'reverse', 'FontSize', 14, 'xscale', 'log', 'ylim', [0 220])
xlabel('^{3}He Concentration')
ylabel('Depth (cm)')
legend('Measured ^{3}He', 'Modeled ^{3}He', 'Modeled ^{3}He profile', 'Location', 'Southeast')
grid on


%% output table

% out = table(samples_LDEO.sample_IDs, samples_LDEO.depth, samples_LDEO.N3, samples_LDEO.dN3, samples_LDEO.N4, samples_LDEO.dN4, samples_LDEO.cronusp_meas, samples_LDEO.N3_stand ,samples_LDEO.N3_cosmo, samples_LDEO.N10, samples_LDEO.dN10, samples_LDEO.N3_cosmo./samples_LDEO.N10, samples_LDEO.N10_decaycorr,  R3_10_corr); 
% out.Properties.VariableNames = {'Sample_ID', 'Sample Depth (cm)', 'Measured [3He](atoms g-1)', '[3He] Uncertainty (atoms g-1)', '[4He] (atoms g-1)', '[4He] Uncertinaty (atoms g-1)', '[3He] of CRONUS-P (atoms g-1)', 'Standardized [3He] (atoms g-1)', 'Cosmogenic [3He] (atoms g-1)', 'Measured [10Be] (atoms g-1)', '[10Be] Uncertainty (atoms g-1)', 'Raw He-3/Be-10 Ratio','Decay-corrected [Be-10 (atoms/g)', 'He-3/Be-10 Ratio'};
% 
% addpath("tables/")
% filename = 'tables/LABCO_tables.xlsx';
% writetable(out,filename,'Sheet',1)
% 
