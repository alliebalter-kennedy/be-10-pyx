% This script performs the fitting exercise for the steady erosion case. 

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

p.texp = 14.5e6;


%% Load data 

[loc data] = create_LABCO_data(constants);

% Be-10 mask
mask10 = ~isnan(data.N10);

%% surface production

p.h = antatm(loc.elevation); % site pressure
SFsp = stone2000(loc.latitude,p.h,1); % scaling factor

p.P3sp = constants.P3_SLHL.*SFsp;

%% Generate muon fluxes and stopping rates at all sample depths

p.mc10.Natoms = 1.5684e22; % for O in average Ferrar Dolerite pyroxenes
p.mc10.sigma0 = 0.280e-30; % ubarns; Balco 2017
p.mc10.k_neg = 1; % Dummy

% mc3.Natoms = 2.61E+22; % for total atoms in augite (average atomic weight ~23)
p.mc3.Natoms = 2.7373e+22; % for total atoms standard basalt (average atomic weight ~22)
p.mc3.sigma0 = 5.70e-30; % ubarns; from Balco fit to Larsen data
p.mc3.k_neg = 1; % Dummy

p.predz = logspace(0,4,100); %[g cm^-2]

% muon fluxes & stopping rate for He-3 and Be-10
for a = 1:length(p.predz)
    % helium
    this_mu3 = P_mu_total_alpha1(p.predz(a),p.h,p.mc3,'yes');
    p.m3stub_fast(a) = this_mu3.P_fast;
    p.m3stub_neg(a) = this_mu3.P_neg;
end

for a = 1:length(p.predz)
    % beryllium 
    this_mu10 = P_mu_total_alpha1(p.predz(a),p.h,p.mc10,'yes');
    p.m10stub_fast(a) = this_mu10.P_fast;
    p.m10stub_neg(a) = this_mu10.P_neg;
end;


%% Initial guess

% Initial guess 

erosionRate0 = 5; % cm/Myr
HeBeRatio0 = 36.9; % from SuLa optimization 
Lsp0 = 140; % g/cm2
k_neg_10_0 = 0.00191 .* 0.704 .*0.1828; % fstar * fC * fD; fC and fD from Heisinger 2002b for O in SiO2. Will fit whole k_neg for pyroxene. 
k_neg_3_0 = 0.0045; % guess from Nesterenok
depthMissing0 = 0; 
nonCosmoHe3_0 = 3.3e6; % from Balco blog post

x0 = [erosionRate0 HeBeRatio0 Lsp0 k_neg_10_0 k_neg_3_0 depthMissing0 nonCosmoHe3_0];
xmin = [0 0 0 0 0 -20 nonCosmoHe3_0];
xmax = [Inf Inf Inf Inf Inf 20 nonCosmoHe3_0];

% fit age only (assume 3 cm missing below 18 cm)
% xmin = [0 HeBeRatio0 Lsp0 k_neg_10_0 0.0045 -3 3.3e6];
% xmax = [Inf HeBeRatio0 Lsp0 k_neg_10_0 0.0045 -3 3.3e6];
%% Try optimizing

opts = optimset('fmincon');
opts = optimset(opts,'display','iter','tolfun',1e-7);

[optx,fval] = fmincon(@(x) LABCO_objective_steadyerosion(x,data,p,0),x0,[],[],[],[],xmin,xmax,[],opts);

opt.erosionRate = optx(1);
opt.HeBeRatio = optx(2);
opt.Lsp = optx(3);
opt.k_neg_10 = optx(4);
opt.k_neg_3 = optx(5);
opt.depthMissing = optx(6);
opt.nonCosmoHe3 = optx(7);

% Report results
disp(['erosion rate = ' sprintf('%0.2f',optx(1)) ' cm Myr^{-1}; He/Be Ratio = ' sprintf('%0.1f',optx(2))]);
disp(['Lsp = ' sprintf('%0.2f',optx(3)) ' g cm^{-2}; k_neg_10 = ' sprintf('%0.5f',optx(4)) '; k_neg_3 = ' sprintf('%0.5f',optx(5))]);
disp(['Depth missing = ' sprintf('%0.2f',optx(6)) ' cm;']);

%% Plot results

% Get diagnostics
optd = LABCO_objective_steadyerosion(optx,data,p,1);

pz = [0:5:275].*constants.rho;

mc10.Natoms = 1.5684e22; % for O in average Ferrar Dolerite pyroxenes
mc10.sigma0 = 0.280e-30; % ubarns; Balco 2017
mc10.k_neg = opt.k_neg_10; % Dummy

% mc3.Natoms = 2.61E+22; % for total atoms in augite (average atomic weight ~23)
mc3.Natoms = 2.7373e+22; % for total atoms standard basalt (average atomic weight ~22)
mc3.sigma0 = 5.70e-30; % ubarns; from Balco fit to Larsen data
mc3.k_neg = opt.k_neg_3; % Dummy


% calculate predicted depth profiles
    %spallation
        pred.N10sp = (optd.P10sp.*exp(-pz./opt.Lsp))./(constants.l10 + ((opt.erosionRate .* constants.rho .* 1e-6)./opt.Lsp)) .* (1-exp(-(constants.l10 + ((opt.erosionRate .* constants.rho .* 1e-6)./opt.Lsp)) .* p.texp));
        pred.N3sp = (p.P3sp.*exp(-pz./opt.Lsp)./((opt.erosionRate .* constants.rho .* 1e-6)./opt.Lsp)) .* (1-exp(-p.texp .* ((opt.erosionRate .* constants.rho .* 1e-6)./opt.Lsp)));

    %muons
        for a = 1:length(pz);
            pred.N10mu(a) = integral(@(t) P_mu_total_alpha1((pz(a)+(opt.erosionRate .* constants.rho .* 1e-6)*t),p.h,mc10),0,p.texp,'reltol',1e-3) .* exp(-constants.l10.*p.texp);
            pred.N3mu(a) = integral(@(t) P_mu_total_alpha1((pz(a)+(opt.erosionRate .* constants.rho .* 1e-6)*t),p.h,mc3),0,p.texp,'reltol',1e-3);
        end;

pred.N10 = pred.N10sp + pred.N10mu;


pred.N3tot = pred.N3sp + pred.N3mu + opt.nonCosmoHe3;
pred.N3cosmo = pred.N3sp + pred.N3mu;


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


txt = {['erosion rate = ' sprintf('%0.2f',optx(1)) ' cm Myr^{-1}'], ...
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
plot(pred.N3cosmo, pz./constants.rho, 'k')
plot(pred.N3tot, pz./constants.rho, 'k:')
% plot(optd.N3tot, zcm, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 6)

set(gca, 'YDir', 'reverse', 'xscale', 'log', 'FontSize', 14, 'ylim', [0 pz(end)./constants.rho])
xlabel('[^{3}He]', 'FontSize', 16)
ylabel('Depth (cm)', 'FontSize', 16)
grid on
box on
legend('LDEO', 'BGC', 'Muons', 'Spallation', 'Total Cosmogenic He-3', 'Total He-3', 'Location', 'Southeast')

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
