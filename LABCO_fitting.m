% this script finds the best fitting combination of erosion rate and age
% for the Labyrinth Ferrar Dolerite using Be-10 and He-3 measurements in 9 
% pyroxene samples. 

% The theoretical depth profiles are calculated using code published
% alongside Schaefer et al. (2016), original scripts by Greg Balco. 

clear all
close all

%% Constants

constants.mask = [1:6 8:9];             % use to exlude samples if necessary

constants.HeBeRatio = 39;            % calculated from SULA samples,
                                        % Balter-Kennedy, in prep.

constants.P3_SLHL = 124.03;             % Production rate of He-3 in
                                        % pyroxene at sea level, high
                                        % latitude; [atoms g^-1 yr^-1];
                                        % calculated using the v3 online
                                        % calculator
loc.erosionRate = [0:1e-6:0.6e-4];     % Subaerial erosion 
                                        % rates to test; [cm yr^-1]. Use 0 
                                        % for Figure 1. 44.1 is steady state

loc.ages = [5e5:1e5:8e6];               % ages to test            
                

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
depthProfile.p.P3sp = depthProfile.constants.P3p_St .* depthProfile.SFsp; % He-3 spallation production rate at surface in pyroxene
depthProfile.p.P10sp = depthProfile.constants.P10q_St .* depthProfile.SFsp; % Be-10 spallation production rate at surface IN QUARTZ
depthProfile.p.P26sp = depthProfile.p.P10sp.*depthProfile.constants.R2610q; % Al-26 spallation production rate at surface
depthProfile.p.P14sp = depthProfile.constants.P14q_St.*depthProfile.SFsp; % C-14 spallation production rate at surface

% Attenuation
depthProfile.p.Lsp = constants.L; % g/cm2.

% percentage of Be-10 production by muons in quartz
depthProfile.muPctBe = depthProfile.m.P10mu(1)./depthProfile.p.P10sp; 

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

for a = 1:length(loc.erosionRate)
    for b = 1:length(loc.ages)

        for i = 1:length(samples.ID)
            model.pfun10 = @(zz, t) PofZ_he_nomu((zz + loc.erosionRate(a).*constants.rho.*t), ...
                depthProfile.m,depthProfile.p,10) .* exp(-depthProfile.constants.l10.*t);
            model.pfun3 = @(zz, t) PofZ_he_nomu((zz + loc.erosionRate(a).*constants.rho.*t), ...
                depthProfile.m,depthProfile.p,3);
        
            model.N10_sampleDepths(i, a, b) = integral2(model.pfun10, ...
                samples.td(i).*constants.rho, samples.bd(i).*constants.rho, ...
                0, loc.ages(b),'RelTol',1e-3,'AbsTol',1e-3) ...
                ./(samples.bd(i).*constants.rho - samples.td(i).*constants.rho);
            model.N3_sampleDepths(i, a, b) = integral2(model.pfun3, ...
                samples.td(i).*constants.rho, samples.bd(i).*constants.rho, ...
                0, loc.ages(b),'RelTol',1e-3,'AbsTol',1e-3) ...
                ./(samples.bd(i).*constants.rho - samples.td(i).*constants.rho);
        end

    end
end

%% calculate reduced chi-squared

model.mismatch = squeeze(sum(((samples.N10 - model.N10_sampleDepths)./...
    samples.dN10).^2, 1)) + squeeze(sum(((samples.N3_cosmogenic - ...
    model.N3_sampleDepths)./samples.dN3).^2, 1));

model.reducedChi2 = model.mismatch./((2.*length(samples.ID))-2);

model.bestFit = min(min(model.reducedChi2));

[model.minErosion, model.minAge] = find(model.reducedChi2 == model.bestFit);

model.bestErosion = loc.erosionRate(model.minErosion);
model.bestAge = loc.ages(model.minAge);

%% Calculate steady state erosion rate

model.steadyErosion = (depthProfile.p.P3sp ./4.7679e8).*constants.L; %He-3 surface conc. from exponential fit to data
model.steadyErosion = (model.steadyErosion ./ constants.rho);
%% Calculate depth profiles for 
for i = 1:length(depthProfile.depth_gcm2)
    model.pfun10 = @(t) PofZ_he_nomu((depthProfile.depth_gcm2(i) + loc.erosionRate(model.minErosion).* ...
        constants.rho.*t),depthProfile.m,depthProfile.p,10) ...
        .* exp(-depthProfile.constants.l10.*t);
    model.pfun3 = @(t) PofZ_he_nomu((depthProfile.depth_gcm2(i) + loc.erosionRate(model.minErosion).* ...
        constants.rho.*t),depthProfile.m,depthProfile.p,3);

    model.N10_line(i, 1) = integral(model.pfun10, 0, loc.ages(model.minAge), ...
        'RelTol',1e-3,'AbsTol',1e-3);
    model.N3_line(i, 1) = integral(model.pfun3, 0, loc.ages(model.minAge), ...
        'RelTol',1e-3,'AbsTol',1e-3);
 end

%% Modeled and measured ratios

samples.HeBeRatio = samples.N3_cosmogenic./samples.N10;
samples.HeBeRatioUncertainty = sqrt((samples.dN3./samples.N3).^2 + ...
    (samples.dN10./samples.N10).^2) .* samples.HeBeRatio;

model.BestFitHeBeRatio = squeeze(model.N3_sampleDepths(:, model.minErosion, model.minAge)) ./ ...
    squeeze(model.N10_sampleDepths(:, model.minErosion, model.minAge));

%% Plot colormap

[map.X, map.Y] = meshgrid(loc.ages, loc.erosionRate);

figure(1)
subplot(1, 4, 1)
pcolor(loc.erosionRate, loc.ages, model.reducedChi2');
hold on
plot([model.steadyErosion model.steadyErosion], [0 10e6], 'k--', 'LineWidth', 1.5)
colormap(parula);
c = colorbar;
caxis([0 30]);
c.Label.String = 'Red-chi^{2}';
c.Label.FontSize = 14;
set(gca, 'xlim', [0 6e-5], 'ylim', [1e6 8e6])
xlabel('Erosion rate (cm/yr)');
ylabel('Age (years)')

shading flat


%% Plot modeled vs. measured concentrations 

figure(1)
subplot(1, 4, 2)
hold on
errorbar(samples.N10, samples.avgDepth, samples.dN10, 'horizontal', 'ks', 'MarkerFaceColor', 'k');
plot(squeeze(model.N10_sampleDepths(:, model.minErosion, model.minAge)), samples.avgDepth, 'rs', 'MarkerFaceColor', 'r')
plot(model.N10_line, depthProfile.depth_cm, 'r')
set(gca, 'YDir', 'reverse', 'FontSize', 14, 'xscale', 'log', 'ylim', [0 220])
xlabel('^{10}Be Concentration')
ylabel('Depth (cm)')
legend('Measured ^{10}Be', 'Modeled ^{10}Be', 'Modeled ^{10}Be profile', 'Location', 'Southeast')
grid on 

figure(1)
subplot(1, 4, 3)
hold on
errorbar(samples.N3_cosmogenic, samples.avgDepth, samples.dN3, 'horizontal', 'ks', 'MarkerFaceColor', 'k');
plot(squeeze(model.N3_sampleDepths(:, model.minErosion, model.minAge)), samples.avgDepth, 'bs', 'MarkerFaceColor', 'b')
plot(model.N3_line, depthProfile.depth_cm, 'b')
set(gca, 'YDir', 'reverse', 'FontSize', 14, 'xscale', 'log', 'ylim', [0 220])
xlabel('^{3}He Concentration')
ylabel('Depth (cm)')
legend('Measured ^{3}He', 'Modeled ^{3}He', 'Modeled ^{3}He profile', 'Location', 'Southeast')
grid on

%% Plot Ratios

figure(1)
subplot(1, 4, 4)
hold on
errorbar(samples.HeBeRatio, samples.avgDepth, samples.HeBeRatioUncertainty, ...
    'horizontal', 'ks', 'MarkerFaceColor', 'k')
plot(model.BestFitHeBeRatio, samples.avgDepth, 'bs', 'LineWidth', 2, 'MarkerFaceColor', 'b')
set(gca, 'XLim', [20 80], 'YLim', [0 200], 'YDir', 'reverse', 'FontSize', 14)
% legend(['^{3}He/^{10}Be ratio' 10 'from measurements'], 'Linear fit to measurements', 'Theoretical ^{3}He/^{10}Be ratio', 'Location', 'NW')
ylabel('Depth (cm)')
xlabel('^{3}He/^{10}Be Ratio')
grid on