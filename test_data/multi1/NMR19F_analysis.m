%%
close all
clear all

%% 1
load data.txt;   % convert test.ft1 file to data.txt using the following command:
                %            pipe2txt.tcl -index ppm test.ft1>data.txt
                % you'll need to have this script in the same directory as
                % the data.txt file

NS = 36 * 2048;   % number of experiments summed times the NS for each experiment
conc = 14;       % conc in uM

normalise = NS * conc;

%% 2
figure
plot(data(:,1),data(:,2));
hold on

baseline_raw = data;                            %select region with signal
baseline_raw(5500:9600,:)=[];
plot(baseline_raw(:,1),baseline_raw(:,2),'r.'); %check signal is not selected

%% 3

%baseline using 4th order polynomial

% Fit: 'untitled fit 1'.

xData = baseline_raw(:,1);
yData = baseline_raw(:,2);

% syms x y

% Set up fittype and options.e

ft = fittype( 'poly4' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

[fitresult, gof] = fit( xData, yData, ft, opts );


fitvalues = coeffvalues(fitresult);         %results of fit

baseline_fit = data;
data_corr = data;
for i=1:length(baseline_fit(:,1)),
    baseline_fit(i,2) = fitvalues(1)*baseline_fit(i,1)^4+ fitvalues(2)*baseline_fit(i,1)^3 + fitvalues(3)*baseline_fit(i,1)^2 + fitvalues(4)*baseline_fit(i,1) +fitvalues(5) ;
    data_corr(i,2) = data(i,2) - baseline_fit(i,2);
end

figure
subplot(3,1,1)
plot(data(:,1),data(:,2),'k-')
set(gca,'xdir','reverse');
    xlim([-70 -50])
    set(gca,'ytick',[])
    set(gca,'tickdir','out')
    set(gca,'box','off')
    hold off
subplot(3,1,2)
% plot(baseline_raw(:,1),baseline_raw(:,2),'k-')
plot(data(:,1),data(:,2),'k-')
hold on
plot(baseline_fit(:,1),baseline_fit(:,2),'r-')
set(gca,'xdir','reverse');
    xlim([-70 -50])
    set(gca,'ytick',[])
    set(gca,'tickdir','out')
    set(gca,'box','off')
    hold off
subplot(3,1,3)
hold on
plot(data_corr(:,1),data_corr(:,2),'k-')
set(gca,'xdir','reverse');
% syms x
y= @(x) x*0;
plot(data_corr(:,1),data_corr(:,2),'k-')
fplot(y)
set(gca,'xdir','reverse');
    xlim([-70 -50])
    set(gca,'ytick',[])
    set(gca,'tickdir','out')
    set(gca,'box','off')
    hold off

data_extract = [data_corr(:,1) data_corr(:,2)./normalise];

% adjust phasing of spectrum and repeat if needed
%%

% Fitting (1 state with residuals) 
data = data_extract; 

clear xData
clear yData
    
for i=1:length(data_extract(:,1)),
    xData(i,1) = data(i,1)+62.2; %shift to set centre between 2 peaks as 0 Hz
end

yData = data(:,2)./max(data(:,2));


syms x y
x1 = -4;
x2 = +4;


% Set up fittype and options. 
ft = fittype( 'h1/(1+((CS1-x)/LW1)^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

    opts.StartPoint = [0.8 0.1 1]; %for shorter length

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
fitvalues = coeffvalues(fitresult);         %results of fit

% Plot fit with data.
figure
subplot(2,1,1)
hold on
h = plot(xData-62.2, yData.*max(data(:,2)));
f1=fitvalues(3).*max(data(:,2))/(1+((fitvalues(1)-62.2-x)/fitvalues(2))^2);
fplot(f1,'b','linewidth',0.5)
set(gca,'xdir','reverse');

xlim([-63.5 -59.5]) %ras F32

% Label axes
xlabel x
grid on

%residuals plot
for i=1:length(yData),
    fitted(i,1) = fitvalues(3)/(1+((fitvalues(1)-xData(i))/fitvalues(2))^2); 
    residual(i,1) = yData(i,1)-fitted(i,1);
end
subplot(2,1,2)
hold on
h = plot(xData-62.2, residual.*max(data(:,2)),'color',[0.1 0.1 0.1].*5);
title('residuals')
set(gca,'xdir','reverse');
xlim([-63.5 -59.5]) %ras F32
    set(gca,'tickdir','out')
    set(gca,'box','off')
% Label axes
xlabel x
ylabel y
grid on


clear fiterror

confidence = confint(fitresult);
[m n] = size(confidence);
for i=1:n,
    fiterror(1,i) = abs(confidence(2,i)-confidence(1,i))/2;       %error of fit values
end

fit_data = fitvalues';
fit_data_err =  fiterror';

CS = [fit_data(1)-62.2 fit_data_err(1)]
LW = abs([fit_data(2) fit_data_err(2)].*470.611*2)
h = [fit_data(3) fit_data_err(3)].*max(data(:,2));

for i=1:1,
    integral(i,1) = LW(i)*h(i);
end

populations_frac = integral./sum(integral)


results = [CS(:,1) LW(:,1) populations_frac]

RMSE = ((sum(residual.^2))/length(residual)^0.5)
%%

p0_1 = [fit_data(3) (fit_data(1)).*470.611*2*pi/1000 abs(fit_data(2)*470.611*2*pi)]
%%

cd /Users/shschan/NMR/RAS_NMR/analysis_lineshape
%error analysis by bootstrapping
xData_bs = xData.*470.611*2*pi; %convert ppm to s-1
yData_bs = yData;  

clear pfit pfitErr
[pfit, pfitErr] = fit_mc_1peak(xData_bs, yData_bs);

%%

fit_bs = reshape(pfit, [3 1])';
fit_bs_err = reshape(pfitErr, [3 1])';

%convert to units
h_bs = [fit_bs(:,1).*max(data(:,2)) fit_bs_err(:,1).*max(data(:,2))]
CS_bs = [fit_bs(:,2)*1000./(470.611*2*pi)-62.2 fit_bs_err(:,2)*1000./(470.611*2*pi)]
LW_bs = [fit_bs(:,3)./pi+10 fit_bs_err(:,3)./pi]

for i=1,
    integral_bs(i,1) = LW_bs(i,1)*h(i,1);
    integral_bs_err(i,1) = ((LW_bs(i,2)/LW_bs(i,1))^2 + (h_bs(i,2)/h_bs(i,1))^2)^0.5 * integral_bs(i,1);
    integral_sqerror(i,1) = (integral_bs_err(i)/integral_bs(i))^2;
end

populations_frac_bs = integral_bs./sum(integral_bs);

for i=1
    populations_frac_bs_err(i,1) = populations_frac_bs(i)*((integral_bs_err(i)/integral_bs(i))^2 + (sqrt(sum(integral_sqerror))/sum(integral_bs))^2)^0.5;
end




CS_fit = [CS(:,1) CS_bs(:,2)]
LW_fit = [LW(:,1) LW_bs(:,2)]
integral_fit = [integral integral_bs_err]
populations_fit = [populations_frac populations_frac_bs_err]


results_fit_total = [CS_fit LW_fit integral_fit./1000 populations_fit]







































%% 4-state fitting

data=data_extract;
clear xData
clear yData
    
for i=1:length(data_extract(:,1)),
    xData(i,1) = data(i,1)+62.2; %shift to set centre between 2 peaks as 0 Hz
 
end

yData = data(:,2)./max(data(:,2));

syms x y
x1 = -4;
x2 = +4;


% Set up fittype and options.
ft = fittype( 'h1/(1+((CS1-x)/LW1)^2)+ h2/(1+((CS2-x)/LW2)^2) + h3/(1+((CS3-x)/LW3)^2) + h4/(1+((CS4-x)/LW4)^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% % 
opts.Lower = [-0.5 -0.5 0 0  0 0 0 0   0 0 0 0]; 
opts.StartPoint = [-0.4 -0.4 0.4 0.4    0.1 0.05 0.1 0.05   1 1 1 1]; %usual

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
fitvalues = coeffvalues(fitresult);         %results of fit

% Plot fit with data.
figure
subplot(2,1,1)
hold on
h = plot(xData-62.2, yData.*max(data(:,2)),'color',[0.1 0.1 0.1].*5)
f1=fitvalues(9).*max(data(:,2))/(1+((fitvalues(1)-62.2-x)/fitvalues(5))^2);
f2=fitvalues(10).*max(data(:,2))/(1+((fitvalues(2)-62.2-x)/fitvalues(6))^2);
f3=fitvalues(11).*max(data(:,2))/(1+((fitvalues(3)-62.2-x)/fitvalues(7))^2);
f4=fitvalues(12).*max(data(:,2))/(1+((fitvalues(4)-62.2-x)/fitvalues(8))^2);
ftotal = f1+f2+f3+f4;
fplot(f1,'b','linewidth',0.5)
fplot(f2,'b','linewidth',0.5)
fplot(f3,'b','linewidth',0.5)
fplot(f4,'b','linewidth',0.5)
fplot(ftotal,'r','linewidth',1)
set(gca,'xdir','reverse');
xlim([-64 -60]) %usual

% Label axes
xlabel x
ylabel y
grid on


%residuals plot
for i=1:length(yData),
    fitted(i,1) = fitvalues(9)/(1+((fitvalues(1)-xData(i))/fitvalues(5))^2)+fitvalues(10)/(1+((fitvalues(2)-xData(i))/fitvalues(6))^2)+fitvalues(11)/(1+((fitvalues(3)-xData(i))/fitvalues(7))^2)+fitvalues(12)/(1+((fitvalues(4)-xData(i))/fitvalues(8))^2); 
    residual(i,1) = yData(i,1)-fitted(i,1);
end
subplot(2,1,2)
hold on
h = plot(xData-62.2, residual.*max(data(:,2)),'color',[0.1 0.1 0.1].*5);
title('residuals')
set(gca,'xdir','reverse');
xlim([-64 -60])
% Label axes
xlabel x
ylabel y
grid on


clear fiterror

confidence = confint(fitresult);
[m n] = size(confidence);
for i=1:n,x
    fiterror(1,i) = abs(confidence(2,i)-confidence(1,i))/2;       %error of fit values
end

fit_data = fitvalues';
fit_data_err =  fiterror';

CS = [fit_data(1:4)-62.2 fit_data_err(1:4)]
LW = abs([fit_data(5:8) fit_data_err(5:8)].*470.611*2)
h = [fit_data(9:12) fit_data_err(9:12)].*max(data(:,2));

for i=1:4,
    integral(i,1) = LW(i)*h(i);
end

populations_frac = integral./sum(integral);

results = [CS(:,1) LW(:,1) populations_frac]

%%
%calculate initial parameters for bootstrapping error analysis
LB = 10;
amp_initial = h(:,1)./max(data(:,2));
freq_initial = (CS(:,1)+62.2).*(470.611*2*pi)./1000;
R2_initial = (LW(:,1)-LB).*pi;
p_initial = [amp_initial freq_initial R2_initial] %multiply R2 by 1000 again
%%
%error analysis by bootstrapping

xData_bs = xData.*470.611*2*pi; %convert ppm to s-1
yData_bs = yData; 
[pfit, pfitErr] = fit_mc_4states(xData_bs, yData_bs);

%% plot fit
fit_bs = reshape(pfit, [3 4])';
fit_bs_err = reshape(pfitErr, [3 4])';
%convert to units
h_bs = [fit_bs(:,1).*max(data(:,2)) fit_bs_err(:,1).*max(data(:,2))]
CS_bs = [fit_bs(:,2)*1000./(470.611*2*pi)-62.2 fit_bs_err(:,2)*1000./(470.611*2*pi)]
LW_bs = [fit_bs(:,3)./pi+40 fit_bs_err(:,3)./pi];
for i=1:4,
    integral_bs(i,1) = LW_bs(i,1)*h(i,1);
    integral_bs_err(i,1) = ((LW_bs(i,2)/LW_bs(i,1))^2 + (h_bs(i,2)/h_bs(i,1))^2)^0.5 * integral_bs(i,1);
    integral_sqerror(i,1) = (integral_bs_err(i)/integral_bs(i))^2;
end
populations_frac_bs = integral_bs./sum(integral_bs);
for i=1:4
    populations_frac_bs_err(i,1) = populations_frac_bs(i)*((integral_bs_err(i)/integral_bs(i))^2 + (sqrt(sum(integral_sqerror))/sum(integral_bs))^2)^0.5;
end