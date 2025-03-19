%cd /Users/shschan/NMR/F19_NMR/emma_19F_trmd_NTD_68TAG_080125/4
clear all;  % Removes all variables and functions from memory
clc;        % Clears the command window
close all; 
load data.txt

gradpt = 6;

diff_raw = reshape(data, [length(data(:,1))./gradpt length(data(1,:)).*gradpt]);

%%
diff = [diff_raw(:,1) diff_raw(:,gradpt*2+1:gradpt*3)];

figure('position',[300 300 800 300])
subplot(1,2,1)
hold on
for i=1:gradpt
plot(diff(:,1),diff(:,i+1),'color',[0 0 1].*i/gradpt)
end
xlim([-63 -62])
    set(gca,'ytick',[])
    set(gca,'xdir','reverse');
    set(gca,'tickdir','out')
    set(gca,'box','off')
    xlabel('^1H chemical shift (ppm)')
    

%%

[a lowerCS] = min(abs(diff(:,1)- -62.65));
[a upperCS] = min(abs(diff(:,1)- -62.49));

[a noiselowerCS] = min(abs(diff(:,1)- -61));
[a noiseupperCS] = min(abs(diff(:,1)- -60.84));

for i=1:gradpt,
    integral(i,1) = sum(diff(upperCS:lowerCS, i+1));
    int_noise(i,1) = sum(diff(noiseupperCS:noiselowerCS, i+1));
end

noise = std(int_noise');

for i=1:gradpt,
    integral_err(i,1) = noise;
end

csvwrite('diff_data.csv', diff);


%%

grad_low = 0.05;
grad_high = 0.95; 

% Gmax = 0.558; % G/cmA for 700
% Gmax = 0.677  %for 800
Gmax = 0.67 %for 500

grad_change = (grad_high-grad_low)/(gradpt-1);

for i=1:gradpt,
    G_percent(i,1) = grad_low+grad_change*(i-1);
    G(i,1) = G_percent(i,1)*Gmax;
    xValue(i,1) =  G(i,1)^2;
    yValue(i,1) = log(integral(i,1)./max(integral));
    yValue_err(i,1) = integral_err(i)/integral(i,1)./max(integral);
end

% figure
% hold on
% errorbar(xValue, yValue,yValue_err,'o-')

syms x
ft = fittype('m*x + c', 'independent', 'x', 'dependent', 'y')
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

opts.Startpoint=[1 0.3]
% Fit model to data.
[fitresult, gof] = fit( xValue, yValue, ft, opts );
fitvalues = coeffvalues(fitresult);         %results of fit


% figure
subplot(1,2,2)
hold on
errorbar(xValue, yValue,yValue_err,'o')
f1 = fitvalues(2)*x +fitvalues(1);
fplot(f1,[0 0.5])
set(gca,'tickdir','out')
xlabel('G^2 (T^2 m^{-2})')
ylabel('ln(I/I_0)')
% xlim([0 0.15])

data_to_save = [real(xValue), real(yValue), abs(real(yValue_err))];
csvwrite('errorbar_data.csv', data_to_save);




%%
gamma = 2.68E+08;
delta = 0.004;
DELTA = 0.1;
sigma = 0.9;
err = confint(fitresult);
prefactor = gamma^2 * delta^2  * sigma^2 * (DELTA-delta/3);

D = -fitvalues(2)/prefactor
D_err = (fitvalues(2)-err(1,2))/prefactor


%%
T = 298; %K
n = 0.89E-3 ; %Pa s


% T = 283; %K
% n = 1.3076E-3 ; %Pa/s

kb = 1.3806488e-23 ; %J K-1

rH = kb*T/(6*pi*n*D) %in metres, Stokes Einstein eq 

rH_err = rH*(D_err/D)



MW = rH*1E9*10;  %kDa for folded

MW_err = rH_err*1E9*10; 


%%

 
figure
hold on
for i=1:gradpt,
    plot3(diff(:,1),[zeros(2048,1)+G(i)], diff(:,i+1),'color',[0 0 1].*i/gradpt)
end
xlim([5 12])
view([50 10])
set(gca,'ztick',[])
xlabel('^1H chemical shift (ppm)')
ylabel('G (T m^{-1})')
set(gca,'xdir','reverse');