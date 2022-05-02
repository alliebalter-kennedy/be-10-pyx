function out = LABCO_objective_zeroerosion(X, data, p, dFlag)

% This is the objective function for fitting exercise 1, which is assuming
% zero erosion and a floating exposure time. 
% out = BCO_objective_zeroerosion(X,d,p,dFlag)
% 
% where X has the following params:
%   X(1) exposure age (Myr)
%   X(2) He-3/Be-10 ratio (non-dimensional)
%   X(3) Lsp (g/cm2) spallogenic e-folding length
%   X(4) k_neg_10 (fraction) negative muon capture cross section
%   X(5) k_neg_3 (fraction) negative muon capture cross section
%   X(6) depthMissing (cm) missing depth below 18 cm
%   X(7) nonCosmoHe3 (atoms g-1) amount of non-cosmogenic He-3
%
% samples is data structure with core data
% p is additional struct with things like surface production rate, etc. 
% dflag is 0 for objective fctn; 1 for diagnostic struct
%
% Modified from script by Greg Balco - Berkeley Geochronology Center - 2018
% Allie Balter-Kennedy - Lamont-Doherty Earth Observatory - 2022
% Not licensed for reuse or distribution

%% Get X in correct dimensions

t = X(1) .* 1e6; % to Ma
HeBeRatio = X(2); % No conversion
Lsp = X(3); % No conversion
k_neg_10 = X(4); % No conversion
k_neg_3 = X(5); % No conversion
depthMissing = X(6); % cm
nonCosmoHe3 = X(7); % Atoms g-1


%% constants

constants.l10 = 4.9975E-07;             % Be-10 decay rate; [yr^-1];
                                        % Nishizumii 2007? Check reference.

constants.rho = 2.94;                   % Density of Ferrar Dolerite; 
                                        % [g cm^-3]

%% Define production

zcm = data.avgDepth;
zcm(zcm>=18) = zcm(zcm>=18) + depthMissing;

zgcm2 = zcm .* constants.rho; % depth, g cm-2

% Spallation

result.P10sp = p.P3sp./HeBeRatio; % ratio is being fit

result.N3sp = p.P3sp .* exp(-zgcm2./Lsp) .* t; % spallogenic He-3 produced at sample depths
result.N10sp = (result.P10sp .* exp(-zgcm2(~isnan(data.N10))./Lsp))./ constants.l10 .* (1-exp(-constants.l10 .*t)); % spallogenic Be-10 produced at sample depths

% Negative muon capture

result.N10neg = (data.stub_mneg10.*k_neg_10)./constants.l10 .* (1-exp(-constants.l10 .*t)); % Be-10 produced by negative muon capture at sample depths
result.N3neg = data.stub_mneg3.*k_neg_3.*t; % He-3 produced by negative muon capture at sample depths

% Fast muons

result.N10fast = data.mfast10./ constants.l10 .* (1-exp(-constants.l10 .*t)); % Be-10 produced by fast muons at sample depths
result.N3fast = data.mfast3 .* t; % He-3 produced by fast muons at sample depths

% Total
result.N10tot = result.N10sp + result.N10neg + result.N10fast;
result.N3tot = result.N3sp + result.N3neg + result.N3fast + nonCosmoHe3;

% Miss

result.miss10 = (result.N10tot - data.N10(~isnan(data.N10)))./(data.dN10(~isnan(data.N10))); 
result.miss3 = (result.N3tot - data.N3_standardized)./(data.dN3);
result.x2 = sum(result.miss10.^2) + sum(result.miss3.^2);

result.depth = zcm;

if dFlag == 0
    out = result.x2;
elseif dFlag == 1
    out = result;
end