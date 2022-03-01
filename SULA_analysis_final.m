%% Purpose

% This script calculates the He-3/Be-10 production ratio in pyroxene using
% samples from Ferrar Dolerite bedrock and boulders collected from several
% locations in the Dry Valleys.

% Inputs are sample information and He-3 and Be-10 concentrations stored in
% the .txt file 'SULA_data.txt'. Outputs include figures and tables for the
% paper Balter-Kennedy et al., (in prep), working title, "10Be analysis in
% pyroxene - a method for routine chemical extraction".

% The surface exposure age of each sample is calculated using v3 of the
% Online Exposure Age Calculator described by Balco et al., (2008). Code
% for sending sample information to the calculator and retrieving exposure
% ages was modified from a script distributed at an ICE-D:GREENLAND
% workshop at the University at Buffalo, 2021. Orginal script by Greg
% Balco.

clear all
close all

%% Constants

constants.P3_SLHL = 124.03;             % Production rate of He-3 in
                                        % pyroxene at sea level, high
                                        % latitude; [atoms g^-1 yr^-1];
                                        % calculated using the v3 online
                                        % calculator

constants.erosionRate = 0e-6;           % Steady-state subaerial erosion 
                                        % rate; [cm yr^-1]. In current
                                        % version, if this is a single
                                        % value, or a vector with length
                                        % equal to the number of samples,
                                        % this will override the input
                                        % datafile. If the length of this
                                        % vector is >1 and doesn't match
                                        % the number of samples, the
                                        % erosion rates from the datafile
                                        % will be used.

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

filename = 'SULA_data.txt';    % file where sample, He-3, and 
                                        % Be-10 data are stored.

data = readtable(filename);             % load data

%% Unpack data into usable form

% put all sample information into structural array for use below.

% put Sample IDs in array
samples.ID = table2cell(data(:, 'Sample_ID'));   % Sample IDs

% put sample data in array
samples.thickness = table2array(data(:, 'thickness'));   % Sample thickness

samples.density = table2array(data(:, 'density'));       % Sample density

samples.shielding = table2array(data(:, 'shielding'));   % Shielding factor

samples.year = table2array(data(:, 'yrCollected'));      % Year of sample 
                                                         % collection

% choose correct erosion rate
% if length(constants.erosionRate) == 1
%     samples.erosionRate(1:length(samples.ID), 1) = constants.erosionRate;
% elseif length(constants.erosionRate) == length(samples.ID)
%     samples.erosionRate = constants.erosionRate;
% else
%     samples.erosionRate = table2array(data(:, 'erosion'));
% end

samples.erosionRate = table2array(data(:, 'erosion'));

% put He-3 data in array
samples.N3 = table2array(data(:, 'N3_LDEO'));         % He-3 concentrations
samples.dN3 = table2array(data(:, 'dN3_LDEO'));       % and uncertainty,
                                                      % measured at LDEO;
                                                      % [atoms g^-1]
                                                      

samples.N3_ETH = table2array(data(:, 'N3_Zurich'));   % He-3 concentrations   
samples.dN3_ETH = table2array(data(:, 'dN3_Zurich')); % and uncertainty,
                                                      % measured at ETH
                                                      % [atoms g^-1]

                                                      
samples.N3_GFZ = table2array(data(:, 'N3_Potsdam'));  % He-3 concentrations   
samples.dN3_GFZ =table2array(data(:, 'dN3_Potsdam')); % and uncertainty,
                                                      % measured at
                                                      % GFZ [atoms g^-1]
                                                      

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

%% Calculate Surface Exposure Ages

% Code in this section modified from script distributed at an ICE-D:GREENLAND
% workshop at the University at Buffalo, 2021. Orginal script by Greg
% Balco. 

% define input text, v3 format for online calculator
calculator.formatInput = ['%s %.2f %.2f %.0f %s %.2f %.2f %.0f %.10f %.0f; ' ...
    '%s He-3 pyroxene %i %i CRONUS-P %i; '];

for i = 1:length(samples.ID)
    calculator.inputText{i} = sprintf(calculator.formatInput, samples.ID{i}, ...
        samples.lat(i), samples.long(i), samples.elv(i), samples.atm{i}, ...
        samples.thickness(i), samples.density(i), samples.shielding(i), ...
        samples.erosionRate(i), samples.year(i), samples.ID{i}, ...
        samples.N3(i), samples.dN3(i), samples.cronusPMeasured(i));
end

calculator.inputText = strjoin(calculator.inputText);

% Send sample info to cosmo calculator
calculator.url = "https://hess.ess.washington.edu/cgi-bin/matweb";
calculator.xmlResult = webread(calculator.url,'mlmfile','age_input_v3', ...
    'reportType','XML', 'resultType','long','plotFlag','no','text_block', ...
    calculator.inputText);

% Load the parser and parse string from data returned from cosmo calculator
import matlab.io.xml.dom.*
calculator.xDoc = parseString(Parser,calculator.xmlResult);

% Get St ages as numbers

calculator.sampleNameList = getElementsByTagName(calculator.xDoc,'sample_name');

calculator.stAgeList = getElementsByTagName(calculator.xDoc,'t3pyroxene_St');
    for i = 1:length(samples.ID)
        samples.StAge(i, 1) = str2num(getTextContent(node(calculator.stAgeList,i)));
    end

calculator.stIntList = getElementsByTagName(calculator.xDoc,'delt3pyroxene_int_St');
    for i = 1:length(samples.ID)
        samples.StInternalUncertainty(i) = str2num(getTextContent(node(calculator.stIntList,i)));
    end

                                                
%% Calculate correction factors

% calculate thickness correction factor. 'thickness' function is from 
% Balco et al. (2008)

correctionFactor.thickness = [1./thickness(samples.thickness, ...
    constants.L, samples.density)]'; 

% calculate CRONUS-P correction factor for standardizing He-3
% concentrations

correctionFactor.cronusP = constants.cronusPAccepted./samples.cronusPMeasured; 

% calculate decay term for Be-10; includes radioactive decay + erosion
% (l + ep/L)
constants.decayTerm = constants.l10 + ((samples.erosionRate .* constants.rho)...
    ./ constants.L);

% calculate decay correction factor (should be > 1). 
correctionFactor.decay = (samples.StAge .* (constants.decayTerm)) ./ ...
    (1 - exp(-(constants.decayTerm).*samples.StAge));

if sum(correctionFactor.decay < 1) ~= 0
    print('WARNING: at least one decay correction factor <1')
end

%% Standardize He-3 for use in production rate calibration

% standardize He-3 concentrations using accepted value for CRONUS-P
samples.N3_standardized = samples.N3 .* correctionFactor.cronusP;

% correct He-3 for erosion if erosion is > 0 cm yr-1

constants.erosionTerm = (samples.erosionRate .* constants.rho)...
    ./ constants.L;

correctionFactor.erosion = (samples.StAge .* (constants.erosionTerm)) ./ ...
    (1 - exp(-(constants.erosionTerm).*samples.StAge));

for i = 1:length(samples.ID)
    if samples.erosionRate == 0

        % correct standardized He-3 concentrations for sample thickness (i.e.,
        % obtain surface concentration)
        samples.N3_standardizedSurface = samples.N3_standardized ...
            .* correctionFactor.thickness;
    else 

        samples.N3_standardizedSurface = samples.N3_standardized .* ...
            correctionFactor.erosion .* correctionFactor.thickness;
    end
end

%% Decay correct Be-10

% perform decay correction on measured Be-10 concentrations
samples.N10_decayCorrected = samples.N10 .* correctionFactor.decay;

% correct decay-corrected Be-10 concentrations for sample thickness (i.e.,
% obtain surface concentration)

samples.N10_decayCorrectedSurface = samples.N10_decayCorrected... 
    .* correctionFactor.thickness;

%% Calculate ratios

% Ratio from measurements - not corrected for thickness, but He-3 is
% standardized to CRONUS-P
samples.measuredHeBeRatio = samples.N3_standardized ./ samples.N10;

% Decay-corrected ratios - not explicitly corrected for thickness, but
% ratio would be the same.
samples.decayCorrectedHeBeRatio = samples.N3_standardizedSurface ./...
    samples.N10_decayCorrected;

%% Perform linear fit for production rate calibration

% use concentrations corrected for thickness. 
% perform linear fit
[linearFit.coefficients, linearFit.stats] = polyfit(samples.N3_standardizedSurface./1e7, ...
    samples.N10_decayCorrectedSurface./1e7, 1);

% assign values for plotting
linearFit.xValuesForPlot = [0:100:60e2];

% calculate confidence intervals
[linearFit.yfit, linearFit.dy] = polyconf(linearFit.coefficients, ...
    linearFit.xValuesForPlot, linearFit.stats, 'predopt', 'curve');

% calculate uncertainy in slope and y intercept
linearFit.coefficientUncertainty = polyparci(linearFit.coefficients, ...
    linearFit.stats);
linearFit.slopeUncertainty = abs(linearFit.coefficientUncertainty(1, 1) - ...
    mean(linearFit.coefficientUncertainty(:, 1)));
linearFit.yInterceptUncertainty = abs(linearFit.coefficientUncertainty(1, 2) - ...
    mean(linearFit.coefficientUncertainty(:, 2)));

% calculate R-squared and coefficient of variance (CV)
linearFit.Rsq = 1 - (linearFit.stats.normr/norm(samples.N10_decayCorrectedSurface ...
    ./1e7 - mean(samples.N10_decayCorrectedSurface./1e7)))^2;
linearFit.slopeCV = linearFit.slopeUncertainty ./ linearFit.coefficients(1) .*100;

% plot
figure(1)
hold on

% plot fit line
plot(linearFit.xValuesForPlot, linearFit.yfit, 'r', 'LineWidth', 1)

% plot 95% confidence intervals
plot(linearFit.xValuesForPlot, linearFit.yfit - linearFit.dy, 'k:')
plot(linearFit.xValuesForPlot, linearFit.yfit + linearFit.dy, 'k:')

% plot data
errorbar(samples.N3_standardizedSurface./1e7, samples.N10_decayCorrectedSurface./1e7, ...
    samples.dN10./1e7, samples.dN10./1e7, samples.dN3./1e7, samples.dN3./1e7, ...
    'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k')

% add sample ID labels
text(samples.N3_standardizedSurface./1e7, (samples.N10_decayCorrectedSurface./1e7)+0.3, ...
    samples.ID,'VerticalAlignment','bottom','HorizontalAlignment','right', ...
    'FontSize', 10)

% add R-squared value
text(250, 3, ['R^{2} = ' sprintf('%0.2f', linearFit.Rsq)], 'FontSize', 12);

% add equation
txt = ['N_{10, decay-corr.} =  (' sprintf('%0.3f', linearFit.coefficients(1))...
    ' ' char(177) ' ' sprintf('%0.3f', linearFit.slopeUncertainty) ') * N_{3} + (' ...
    sprintf('%0.1f', linearFit.coefficients(2)) ' ' char(177) ' ' ...
    sprintf('%0.1f', linearFit.yInterceptUncertainty) ')'];

text(250, 4.5, txt,'FontSize', 12);

% set axes and legend
set(gca, 'FontSize', 16, 'ylim', [-1 22])
legend('Linear Fit', '','', 'Measured ^{10}Be (decay-corrected) and ^{3}He concentrations', ...
    'Location', 'Southeast', 'Fontsize', 10)
ylabel('Decay-corrected ^{10}Be concentrations (10^{7} atoms/g)', 'FontSize', 14);
xlabel('^{3}He concentrations (10^{7} atoms/g)', 'FontSize', 14);
title('SULA Samples')


%% Monte Carlo fit

% To take advantage of the low measurement uncertainty, we use a Monte
% Carlo simulaton to quantify overall uncertainty in the He-3/Be-10 ratio,
% and therefore the Be-10 production rate.


% number of simulations
monteCarlo.numSim = 10000;

for a = 1:monteCarlo.numSim
        % create empty vectors for 
        monteCarlo.N3 = zeros(length(samples.ID), 1);   
        monteCarlo.N10 = zeros(length(samples.ID), 1);

    for i = 1:length(samples.ID)
        % fill vectors with random concentration using 1-sig measurement
        % uncertainty
        monteCarlo.N3(i) = normrnd(samples.N3_standardizedSurface(i), samples.dN3(i));
        monteCarlo.N10(i) = normrnd(samples.N10_decayCorrectedSurface(i), samples.dN10(i));
    end
    
    % linear fit for each simulation
    monteCarlo.coefficients(a, :) = polyfit(monteCarlo.N3, monteCarlo.N10, 1);

end

% make histogram of slopes
figure(2)
hist(monteCarlo.coefficients(:, 1))

% get mean and 95% confidence intervals on slope and y-intercept

monteCarlo.slopeMean = mean(monteCarlo.coefficients(:, 1));   % mean of slopes across all sims
monteCarlo.slope95CI = 2*std(monteCarlo.coefficients(:, 1));  % 2 sigma of slopes across all sims
monteCarlo.slope_CV = (monteCarlo.slope95CI./monteCarlo.slopeMean)*100;  % coefficient of variance

monteCarlo.yIntMean = mean(monteCarlo.coefficients(:, 2));    % mean of y intercepts 
monteCarlo.yInt95CI = 2*std(monteCarlo.coefficients(:, 2));   % 2 sigma of y intercepts
monteCarlo.yIntCV = (monteCarlo.yInt95CI./monteCarlo.yIntMean)*100;   % coefficient of variance

%% Plot Monte Carlo results

figure(3)

for a = 1:1000
    h = plot(linearFit.xValuesForPlot, (monteCarlo.coefficients(a, 1).*linearFit.xValuesForPlot)+(monteCarlo.coefficients(a, 2)./1e7), 'Color', [0 0 0 0.007]);
    hold on
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

plot(linearFit.xValuesForPlot, (monteCarlo.slopeMean.*linearFit.xValuesForPlot)+(monteCarlo.yIntMean./1e7), 'r', 'LineWidth', 1)
hold on

errorbar(samples.N3_standardizedSurface./1e7, samples.N10_decayCorrectedSurface./1e7, samples.dN10./1e7, samples.dN10./1e7, samples.dN3./1e7, samples.dN3./1e7, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k')
text(samples.N3_standardizedSurface./1e7, (samples.N10_decayCorrectedSurface./1e7)+0.5, samples.ID,'VerticalAlignment','bottom','HorizontalAlignment','right', 'FontSize', 10)

set(gca, 'FontSize', 16, 'ylim', [-1 21])
legend('Linear Fit', 'Measured ^{10}Be (decay-corrected) and ^{3}He concentrations', 'Location', 'Southeast', 'Fontsize', 10)
ylabel('Decay-corrected ^{10}Be concentrations (10^{7} atoms/g)', 'FontSize', 14);
xlabel('^{3}He concentrations (10^{7} atoms/g)', 'FontSize', 14);
title('SULA Samples')

% equation from MC fit
txt = ['N_{10, decay-corr.} =  (' sprintf('%0.3f', monteCarlo.slopeMean) ' ' char(177) ' ' sprintf('%0.3f', monteCarlo.slope95CI) ') * N_{3} + (' sprintf('%0.1f', monteCarlo.yIntMean./1e7) ' ' char(177) ' ' sprintf('%0.1f', monteCarlo.yInt95CI./1e7) ')'];

hold on
text(200, 3, txt,'FontSize', 12);
hold on

%% Calculate Be-10 production rate

monteCarlo.HeBeRatio = 1./monteCarlo.slopeMean;

samples.P10 = (1./samples.decayCorrectedHeBeRatio) .* constants.P3_SLHL;

monteCarlo.P10 = monteCarlo.slopeMean .* constants.P3_SLHL;



