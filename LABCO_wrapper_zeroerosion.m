% This script performs the fitting exercise for the zero erosion case. 

% Allie Balter-Kennedy - Lamont-Doherty Earth Observatory - March 2022
% Utilizes scripts written by Greg Balco, with permission

clear all
close all

%% Define constants

%% constants

constants.P3_SLHL = 124.03;             % Production rate of He-3 in
                                        % pyroxene at sea level, high
                                        % latitude; [atoms g^-1 yr^-1];
                                        % calculated using the v3 online
                                        % calculator

constants.l10 = 4.9975E-07;             % Be-10 decay rate; [yr^-1];
                                        % Nishizumii 2007? Check reference.

constants.rho = 2.94;                   % Density of Ferrar Dolerite; 
                                        % [g cm^-3]

constants.Lsp = 140;                    % spallation attenuation length for 
                                        % sample thickness corrections.
                                        % Should be appropriate for
                                        % Antarctica.

constants.cronusPAccepted = 5.02E+09;   % Accepted CRONUS-P concentration; 
                                        % Blard et al. 2015


%% Load data 

[loc data] = create_LABCO_data(constants);

% Be-10 mask
mask10 = ~isnan(data.N10);

%% surface production

h = antatm(loc.elevation); % site pressure
SFsp = stone2000(loc.latitude,h,1); % scaling factor

p.P3sp = constants.P3_SLHL.*SFsp;

%% Generate muon fluxes and stopping rates at all sample depths

mc10.Natoms = 1.5684e22; % for O in average Ferrar Dolerite pyroxenes
mc10.sigma0 = 0.280e-30; % ubarns; Balco 2017
mc10.k_neg = 1; % Dummy

% mc3.Natoms = 2.61E+22; % for total atoms in augite (average atomic weight ~23)
mc3.Natoms = 2.7373e+22; % for total atoms standard basalt (average atomic weight ~22)
mc3.sigma0 = 5.70e-30; % ubarns; from Balco fit to Larsen data
mc3.k_neg = 1; % Dummy

% muon fluxes & stopping rate for He-3
for a = 1:length(data.avgDepth)
    m3 = P_mu_total_alpha1(data.avgDepth(a).*constants.rho,h,mc3,'yes');    
    data.mfast3(a) = m3.P_fast; % unscaled PRs at different depths
    data.stub_mneg3(a) = m3.P_neg;
end;

% muon fluxes & stopping rate for Be-10
for a = 1:length(data.avgDepth(mask10))
    m10 = P_mu_total_alpha1(data.avgDepth(a).*constants.rho,h,mc10,'yes');
    data.mfast10(a) = m10.P_fast; % unscaled PRs at different depths
    data.stub_mneg10(a) = m10.P_neg; % same as above, but in correct dimension for N10
end

data.mfast10 = data.mfast10'; 
data.mfast3 = data.mfast3'; 

data.stub_mneg10 = data.stub_mneg10';
data.stub_mneg3 = data.stub_mneg3';

%% Initial guess

% Initial guess 

t0 = 1; % Ma
HeBeRatio0 = 36.9; % from SuLa optimization 
Lsp0 = 140; % g/cm2
k_neg_10_0 = 0.00191 .* 0.704 .*0.1828; % fstar * fC * fD; fC and fD from Heisinger 2002b for O in SiO2. Will fit whole k_neg for pyroxene. 
k_neg_3_0 = 0.0045; % guess from Nesterenok
depthMissing0 = 0; 
nonCosmoHe3_0 = 3.3e6; % from Balco blog post

x0 = [t0 HeBeRatio0 Lsp0 k_neg_10_0 k_neg_3_0 depthMissing0 nonCosmoHe3_0];
xmin = [0 0 0 0 0 -Inf nonCosmoHe3_0];
xmax = [Inf Inf Inf Inf Inf 20 nonCosmoHe3_0];

% fit age only (assume 3 cm missing below 18 cm)
% xmin = [0 HeBeRatio0 Lsp0 k_neg_10_0 0.0045 -3 3.3e6];
% xmax = [Inf HeBeRatio0 Lsp0 k_neg_10_0 0.0045 -3 3.3e6];
%% Try optimizing

opts = optimset('fmincon');
opts = optimset(opts,'display','iter','tolfun',1e-7);

[optx,fval] = fmincon(@(x) LABCO_objective_zeroerosion(x,data,p,0),x0,[],[],[],[],xmin,xmax,[],opts);

opt.t = optx(1).*1e6;
opt.HeBeRatio = optx(2);
opt.Lsp = optx(3);
opt.k_neg_10 = optx(4);
opt.k_neg_3 = optx(5);
opt.depthMissing = optx(6);
opt.nonCosmoHe3 = optx(7);

% Report results
disp(['texp = ' sprintf('%0.2f',optx(1)) ' Ma; He/Be Ratio = ' sprintf('%0.1f',optx(2))]);
disp(['Lsp = ' sprintf('%0.2f',optx(3)) ' g cm^{-2}; k_neg_10 = ' sprintf('%0.5f',optx(4)) '; k_neg_3 = ' sprintf('%0.5f',optx(5))]);
disp(['Depth missing = ' sprintf('%0.2f',optx(6)) ' cm;']);

%% Plot results

% Get diagnostics
optd = LABCO_objective_zeroerosion(optx,data,p,1);

pz = [0:275].*constants.rho;

% He-3
for a = 1:length(pz);
    m3 = P_mu_total_alpha1(pz(a),h,mc3,'yes');    
    pstub_mfast3(a) = m3.P_fast;
    pstub_mneg3(a) = m3.P_neg;
end;

% Be-10
for a = 1:length(pz);
    m10 = P_mu_total_alpha1(pz(a),h,mc10,'yes');    
    pstub_mfast10(a) = m10.P_fast;
    pstub_mneg10(a) = m10.P_neg;
end;

% calculate predicted depth profiles
pred.N10sp = ((optd.P10sp.*exp(-pz./opt.Lsp))./constants.l10) .* (1-exp(-constants.l10 .* opt.t));
pred.N10fast = ((pstub_mfast10.*mc10.sigma0)./constants.l10) .* (1-exp(-constants.l10.*opt.t));
pred.N10neg = ((pstub_mneg10.*opt.k_neg_10)./constants.l10) .* (1-exp(-constants.l10.*opt.t));
pred.N10mu = pred.N10fast + pred.N10neg;
pred.N10 = pred.N10sp + pred.N10fast + pred.N10neg;

pred.N3sp = opt.t.*p.P3sp.*exp(-pz./opt.Lsp);
pred.N3fast = opt.t.*pstub_mfast3.*mc3.sigma0;
pred.N3neg = opt.t.*pstub_mneg3.*opt.k_neg_3;
pred.N3mu = pred.N3fast + pred.N3neg;
pred.N3 = pred.N3sp + pred.N3fast + pred.N3neg;


%% Plots

% set up
ldeo = find(strcmp(data.lab, 'LDEO'));
bgc = find(strcmp(data.lab, 'BGC'));

lsidex = 0.1;
lsidew = 0.53;
rsidex = 0.68;
rsidew = 0.28;
ht1 = 0.4;
bot1 = 0.08;
ht2 = 0.4;
bot2 = bot1 + ht1 + 0.1; 

ax1 = axes('pos',[lsidex bot1 lsidew ht1]);
ax11 = axes('pos',[rsidex bot1 rsidew ht1]);

ax2 = axes('pos',[lsidex bot2 lsidew ht2]);
ax21 = axes('pos',[rsidex bot2 rsidew ht2]);

% new sample depths 

zcm = data.avgDepth;
zcm(data.avgDepth>=18) = zcm(data.avgDepth>=18) + opt.depthMissing;

%% Plot concentrations

% beryllium
figure(1)
axes(ax2)
hold on
plot(pred.N10mu, pz./constants.rho, '--', 'Color', [0.5 0.5 0.5])
plot(pred.N10sp, pz./constants.rho, 'k--')
plot(pred.N10, pz./constants.rho, 'k')
errorbar(data.N10(~isnan(data.N10)), zcm(~isnan(data.N10)), data.dN10(~isnan(data.N10)), 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
% plot(optd.N10tot, zcm(~isnan(data.N10)), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 6)

set(gca, 'YDir', 'reverse', 'xscale', 'log', 'FontSize', 14, 'ylim', [0 pz(end)./constants.rho])
xlabel('[^{10}Be]', 'FontSize', 16)
ylabel('Depth (cm)', 'FontSize', 16)
grid on
box on

txt = {['exposure age = ' sprintf('%0.2f',optx(1)) ' Ma'], ...
    ['He/Be Ratio = ' sprintf('%0.1f',optx(2))], ['Lsp = ' sprintf('%0.2f',optx(3)) ' g cm^{-2}'] ...
    ['kneg_{10} = ' sprintf('%0.5f',optx(4))], ['kneg_3 = ' sprintf('%0.5f',optx(5))], ...
    ['Depth missing = ' sprintf('%0.2f',optx(6)) ' cm'], ...
    ['\chi^{2} = ' sprintf('%0.1f', fval)]};

text(1.05e6, 175, txt)

axes(ax21)
hold on
plot([0 0], [0 pz(end)./constants.rho], 'k-')
errorbar(optd.miss10, zcm(ldeo), ones(length(ldeo), 1), 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)

set(gca, 'YDir', 'reverse', 'FontSize', 14, 'ylim', [0 pz(end)./constants.rho])
xlabel('Normalized residual', 'FontSize', 16)
yticklabels([])
grid on
box on

% helium
axes(ax1)
hold on
errorbar(data.N3_standardized(ldeo), zcm(ldeo), data.dN3(ldeo), 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(data.N3_standardized(bgc), zcm(bgc), data.dN3(bgc), 'horizontal', 'k^', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
plot(pred.N3mu, pz./constants.rho, '--', 'Color', [0.5 0.5 0.5])
plot(pred.N3sp, pz./constants.rho, 'k--')
plot(pred.N3, pz./constants.rho, 'k')

% plot(optd.N3tot, zcm, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 6)

set(gca, 'YDir', 'reverse', 'xscale', 'log', 'FontSize', 14, 'ylim', [0 pz(end)./constants.rho])
xlabel('[^{3}He]', 'FontSize', 16)
ylabel('Depth (cm)', 'FontSize', 16)
grid on
box on
legend('LDEO', 'BGC', 'Muons', 'Spallation', 'Total', 'Location', 'Southeast')

axes(ax11)
hold on
plot([0 0], [0 pz(end)./constants.rho], 'k-')
errorbar(optd.miss3(ldeo), zcm(ldeo), ones(length(ldeo), 1), 'horizontal', 'ks', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)
errorbar(optd.miss3(bgc), zcm(bgc), ones(length(bgc), 1), 'horizontal', 'k^', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 8)

set(gca, 'YDir', 'reverse', 'FontSize', 14, 'ylim', [0 pz(end)./constants.rho])
xlabel('Normalized residual', 'FontSize', 16)
yticklabels([])
grid on
box on
