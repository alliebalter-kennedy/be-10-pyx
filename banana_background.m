% This script creates the background for a Be-10/He-3 banana plot. 

function [bananaPlot] = banana_background(HeBeRatio, Lsp, l10, rho, erosion, constantAge);

% [figHandle productionRates.constants.P10q_St]
%% Set up production rates
    addpath(['/Users/alexandrabalter/MATLAB/projects/bedrock_core/Common/' ...
        'depth_profiles'])
    
    latitude = 85.5;
    longitude = 0;
    elevation = 0;
    
    productionRates.constants = bedrock_constants();
    
    % Atmospheric pressure at site
    productionRates.siteP = ERA40atm(latitude,longitude ...
        ,elevation); % site air pressure
    
    % Build and load muon profile
    % build_muon_profile.m builds a production rate profile defined on a grid
    % for efficient integration later after Balco, 2017. 
    productionRates.m = build_muon_profile_w14c(productionRates.siteP, ...
        productionRates.constants,0);
    
    % Define production rate info
    productionRates.SFsp = 1; % scaling factor
    productionRates.p.Lsp = Lsp;
    
    productionRates.p.P3sp = productionRates.constants.P3p_St .* productionRates.SFsp; 
    productionRates.p.P10sp = (productionRates.p.P3sp ./ HeBeRatio) - productionRates.m.P10mu(1); 

    %% line of constant exposure
    
    for i = 1:length(constantAge)
        model.pfun10 = @(t) PofZ_he_nomu(0.*t, productionRates.m,...
            productionRates.p, 10).* exp(-l10.*t);
        model.pfun3 = @(t) PofZ_he_nomu(0.*t, productionRates.m, productionRates.p, 3);
    
        N10_constantExposure(i, 1) = integral(model.pfun10, 0, constantAge(i))./(productionRates.p.P10sp + productionRates.m.P10mu(1));
        N3_constantExposure(i, 1) = integral(model.pfun3, 0, constantAge(i))./productionRates.p.P3sp;
    end

    ratioConstantExposure = N10_constantExposure./N3_constantExposure;

    %% line of steady erosion

    for i = 1:length(erosion)
        model.pfun10 = @(t) PofZ_he_nomu((erosion(i).*rho.*t), ...
            productionRates.m,productionRates.p, 10).* exp(-l10.*t);
        model.pfun3 = @(t) PofZ_he_nomu((erosion(i).*rho.*t), ...
            productionRates.m,productionRates.p,3);
        
        N10_steadyErosion(i, 1) = integral(model.pfun10, 0, 100e6)./(productionRates.p.P10sp + productionRates.m.P10mu(1));
        N3_steadyErosion(i, 1) = integral(model.pfun3, 0, 100e6)./productionRates.p.P3sp;
    end

    ratioSteadyErosion = (N10_steadyErosion./N3_steadyErosion);

    %% lines of constant age and exposure

    constantAge = [0:0.5e6:5e6];
    constantErosion = [5e-6:5e-6:5.5e-5];
    
    for a = 1:length(constantErosion)
        for b = 1:length(constantAge)
            model.pfun10 = @(t) PofZ_he_nomu((constantErosion(a).*rho.*t), ...
                productionRates.m,productionRates.p, 10).* exp(-l10.*t);
            model.pfun3 = @(t) PofZ_he_nomu((constantErosion(a).*rho.*t), ...
                productionRates.m,productionRates.p,3);

            N10_ageErosion = integral(model.pfun10, 0, constantAge(b))./(productionRates.p.P10sp + productionRates.m.P10mu(1));
            N3_ageErosion(a, b) = integral(model.pfun3, 0, constantAge(b))./productionRates.p.P3sp;;
            ratioAgeErosion(a, b) = N10_ageErosion ./ N3_ageErosion(a, b);

        end
    end
     

    %% plot
    
    bananaPlot = figure;
    hold on

    plot(N3_constantExposure, ratioConstantExposure, 'k', 'LineWidth', 1.5)
    plot(N3_steadyErosion, ratioSteadyErosion, 'r', 'LineWidth', 1.5)
    % plot lines of constant erosion
    for i = 1:length(constantErosion)
        plot(N3_ageErosion(i, :), ratioAgeErosion(i, :), 'b:')
        hold on
    end
    % plot lines of constant age
    for i = 1:length(constantAge)
        plot(N3_ageErosion(:, i), ratioAgeErosion(:, i), 'r:')
        hold on
    end

    xlabel('He-3 Concentration (at g^{-1} yr^{-1})', 'FontSize', 14)
    ylabel('Be-10/He-3', 'FontSize', 14)
    grid on
    box on

end 