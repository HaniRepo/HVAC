clear all;
close all;
clc;

%% Load Telemetry Data

% Identification Data
data4id.Comp_freq     = xlsread('Data\Comp_Freq.csv','B9220:B217381'); 
data4id.Tw_wex_out    = xlsread('Data\PT.csv','D9220:D217381');  % PT5 
data4id.Tw_wex_in     = xlsread('Data\PT.csv','E9220:E217381');  % PT6
data4id.Tw_wex_delta  = xlsread('Data\PT.csv','H9220:H217381');  % PT5 - PT6
data4id.Air_Temp      = xlsread('Data\MT.csv','G9220:G217381');  % MT6   
data4id.TR_wex_out    = xlsread('Data\MT.csv','D9220:D217381');  % MT4

% Validation Data
data4vl.Comp_freq     = xlsread('Data\Comp_Freq.csv','B95000:B117381');
data4vl.Tw_wex_out    = xlsread('Data\PT.csv','D95000:D117381');  
data4vl.Tw_wex_in     = xlsread('Data\PT.csv','E95000:E117381');
data4vl.Tw_wex_delta  = xlsread('Data\PT.csv','H95000:H117381');
data4vl.TR_wex_out    = xlsread('Data\MT.csv','D95000:D117381');  
data4vl.Air_Temp      = xlsread('Data\MT.csv','G95000:G117381');

%% System Identification Parameters
sampleTime = 0.05; % Sampling time
filter_wex_test = 0; % Set to 1 for filtered validation data
opt_method = 1; % Optimization method selection
nx = 4;  % Number of state variables

%% Run System Identification
output_sys_m = sysid(sampleTime, filter_wex_test, opt_method, nx, data4id, data4vl);

%% System Identification Function
type function sysid(sampleTime, filter_wex_test, opt_method, nx, data4id, data4vl)

function innova_sys_m = sysid(sampleTime, filter_wex_test, opt_method, nx, data4id, data4vl)

% Extract data from structures
Comp_freq    = data4id.Comp_freq;
Tw_wex_out   = data4id.Tw_wex_out;    
Tw_wex_in    = data4id.Tw_wex_in;     
Tw_wex_delta = data4id.Tw_wex_delta; 
Air_Temp     = data4id.Air_Temp;     
TR_wex_out   = data4id.TR_wex_out;  

Comp_freq2   = data4vl.Comp_freq;     
Tw_wex_out2  = data4vl.Tw_wex_out;    
Tw_wex_in2   = data4vl.Tw_wex_in;  
Tw_wex_delta2= data4vl.Tw_wex_delta;  
TR_wex_out2  = data4vl.TR_wex_out;    
Air_Temp2    = data4vl.Air_Temp;        

%% Data Filtering
for i = 1:length(Comp_freq)
   if Comp_freq(i) == 0
       Tw_wex_delta(i) = 0;
       Tw_wex_out(i) = Tw_wex_in(i);
   end
end

if filter_wex_test == 1
    for i = 1:length(Comp_freq2)
        if Comp_freq2(i) == 0
            Tw_wex_delta2(i) = 0;
            Tw_wex_out2(i) = Tw_wex_in2(i);
        end
    end
end

%% System Identification
hvac_a = iddata(Tw_wex_out, [Comp_freq , Tw_wex_in, Air_Temp], sampleTime);
hvac_a.OutputName = {'Tw_wex_out'};
hvac_a.InputName  = {'Comp Freq' , 'Tw_wex_in', 'Air temp'};

mpf = fact(hvac_a, sampleTime, nx, opt_method);

%% Plot Results
figure;
plot(Comp_freq, Tw_wex_delta, '.' );
xlabel('Comp Freq'); ylabel('Water Temp WEx Delta');
grid on;

%% System Results
Final_pred_err  = mpf.Report.Fit.FPE; 
Mean_sq_err     = mpf.Report.Fit.MSE;
FitPercent      = mpf.Report.Fit.FitPercent;

%% Controllability & Observability Checks
if rank(ctrb(mpf.A,mpf.B)) == nx
    disp('The obtained system is Controllable');
end

if rank(obsv(mpf.A,mpf.C)) == nx
    disp('The obtained system is Observable');    
end  

%% Display Results
disp('Final Prediction Error:'); disp(Final_pred_err);
disp('Mean Squared Error:'); disp(Mean_sq_err);
disp('Fit Percentage:'); disp(FitPercent);

%% Output Function Results
innova_sys_m.a = Tw_wex_delta;
innova_sys_m.b = Tw_wex_out2;
innova_sys_m.Final_pred_err = Final_pred_err;
innova_sys_m.Mean_sq_err = Mean_sq_err;
innova_sys_m.FitPercent = FitPercent;

end

%% Optimization Function
function mpf = fact(hvac_a, sampleTime, nx, opt_method)
switch opt_method
    case 1
        opt = ssestOptions('InitializeMethod', 'n4sid', 'Focus', 'prediction', 'SearchMethod', 'auto');
    case 2
        opt = ssestOptions('InitializeMethod', 'n4sid', 'Focus', 'prediction', 'SearchMethod', 'gn');
    case 3
        opt = ssestOptions('InitializeMethod', 'n4sid', 'Focus', 'prediction', 'SearchMethod', 'lm');
    case 4
        opt = ssestOptions('InitializeMethod', 'n4sid', 'Focus', 'prediction', 'SearchMethod', 'grad');
end

mpf = ssest(hvac_a, nx, 'Ts', sampleTime, opt);
end
