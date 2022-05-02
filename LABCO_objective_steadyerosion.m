function out = LABCO_objective_steadyerosion(X, data, p, dFlag)

% This is the objective function for fitting exercise 2, which is assuming
% steady erosion (rate floats) and a fixed exposure time. 
% out = LABCO_objective_zeroerosion(X,d,p,dFlag)
% 
% where X has the following params:
%   X(1) erosion rate (cm Myr-1)
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

%% constants

constants.l10 = 4.9975E-07;             % Be-10 decay rate; [yr^-1];
                                        % Nishizumii 2007? Check reference.

constants.rho = 2.94;                   % Density of Ferrar Dolerite; 
                                        % [g cm^-3]

%% Get X in correct dimensions

erosionRate = X(1) .* constants.rho .* 1e-6; % g/cm2/yr
HeBeRatio = X(2); % No conversion
Lsp = X(3); % No conversion
k_neg_10 = X(4); % No conversion
k_neg_3 = X(5); % No conversion
depthMissing = X(6); % cm
nonCosmoHe3 = X(7); % Atoms g-1

%% Define production

zcm = data.avgDepth;
zcm(zcm>=18) = zcm(zcm>=18) + depthMissing;

zgcm2 = zcm .* constants.rho; % depth, g cm-2

% Spallation

result.P10sp = p.P3sp./HeBeRatio; % ratio is being fit

result.N3sp = exp(-zgcm2./Lsp).*(p.P3sp.*Lsp./erosionRate).*(1-exp(-p.texp.*(erosionRate./Lsp))); % spallogenic He-3 produced at sample depths
result.N10sp = ((result.P10sp .* exp(-zgcm2(~isnan(data.N10))./Lsp))./(constants.l10 + erosionRate./Lsp)).*(1-exp(-(constants.l10 + erosionRate./Lsp).*p.texp)); % spallogenic Be-10 produced at sample depths

% muons

% helium-3
for a = 1:length(data.N3_standardized)
    result.N3fast(a) = integral(@(t) interp1(log(p.predz),p.m3stub_fast,log(zgcm2(a)+(erosionRate.*t))),0,p.texp,'reltol',1e-3,'abstol',1e-1); % He-3 produced by fast muons at sample depths
    result.N3neg(a) = integral(@(t) k_neg_3.*interp1(log(p.predz),p.m3stub_neg,log(zgcm2(a)+(erosionRate.*t))),0,p.texp,'reltol',1e-3,'abstol',1e-1); % He-3 produced by negative muon capture at sample depths
end

% beryllium-10
for a = 1:length(data.N10(~isnan(data.N10)))
    result.N10fast(a) = integral(@(t) interp1(log(p.predz),p.m10stub_fast,log(zgcm2(a)+erosionRate.*t)),0,p.texp,'reltol',1e-3,'abstol',1e-1) .* exp(-constants.l10.*p.texp); % Be-10 produced by fast muons at sample depths
    result.N10neg(a) = integral(@(t) k_neg_10.*interp1(log(p.predz),p.m10stub_neg,log(zgcm2(a)+erosionRate.*t)),0,p.texp,'reltol',1e-3,'abstol',1e-1) .* exp(-constants.l10.*p.texp); % Be-10 produced by negative muon capture at sample depths
end


% Total
result.N10tot = result.N10sp + result.N10neg' + result.N10fast';
result.N3tot = result.N3sp + result.N3neg' + result.N3fast' + nonCosmoHe3;

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