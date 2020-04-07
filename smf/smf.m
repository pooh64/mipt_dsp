function [alfa, D, beta2, gamma] = smf(lambda);
% световод

%% потери (затухание)
alfa_dB = 0.2; % дБ/км
alfa = 0.1*alfa_dB/log10(exp(1)); % 1/км

%% дисперсия
S0 = 0.09; % пс2/(нм*км)
lambda0 = 1312; % нм
D = 0.25*S0*(lambda-lambda0^4/lambda^3); % пс/(нм*км)
beta2 = -D*lambda^2/(6*pi*10^5);

%% нелинейность показателя преломления
n2_Aeff = 2.79; % 1/(Вт)
gamma = 2*pi*n2_Aeff/(0.01*lambda); % 1/(Вт*км)
endfunction
