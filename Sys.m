clear all;
close all;
clc;

%%  Telemetry data


% Comp_freq     = xlsread('Comp_Freq.csv','B44500:B46500') ;
% Tw_wex_out    = xlsread('Data\PT.csv','D344500:D46500') ;  %PT5
% Tw_wex_in     = xlsread('PT.csv','E44500:E46500') ;  %PT6
% Tw_wex_delta  = xlsread('PT.csv','H44500:H46500') ;  % PT5 - PT6



% External (ambient) air temperature MT6 G

%Data for Identification%
data4id.Comp_freq     = xlsread('Data\Comp_Freq.csv','B9220:B217381') ; 
data4id.Tw_wex_out    = xlsread('Data\PT.csv','D9220:D217381') ;  %PT5 
data4id.Tw_wex_in     = xlsread('Data\PT.csv','E9220:E217381') ;  %PT6
data4id.Tw_wex_delta  = xlsread('Data\PT.csv','H9220:H217381') ;  % PT5 - PT6
data4id.Air_Temp      = xlsread('Data\MT.csv','G9220:G217381') ;  % MT6   
data4id.TR_wex_out    = xlsread('Data\MT.csv','D9220:D217381') ;  %MT4


% Comp_freq2     = xlsread('Comp_Freq.csv','B46511:B47912');
% Tw_wex_out2    = xlsread('PT.csv','D46511:D47912');  %PT5
% Tw_wex_in2     = xlsread('PT.csv','E46511:E47912');  %PT6
% Tw_wex_delta2  = xlsread('PT.csv','H46511:H47912');  % PT5 - PT6

%Data for Validation%
data4vl.Comp_freq2     = xlsread('Data\Comp_Freq.csv','B95000:B117381');
data4vl.Tw_wex_out2    = xlsread('Data\PT.csv','D95000:D117381');   %PT5
data4vl.Tw_wex_in2     = xlsread('Data\PT.csv','E95000:E117381');   %PT6
data4vl.Tw_wex_delta2  = xlsread('Data\PT.csv','H95000:H117381');   % PT5 - PT6
data4vl.TR_wex_out2    = xlsread('Data\MT.csv','D95000:D117381');   %MT4
data4vl.Air_Temp2      = xlsread('Data\MT.csv','G95000:G117381') ;  % MT6   


%% Params
sampleTime = 0.05; % Ts
filter_wex_test = 0 ; % Set the new filter for validation part to 1 if you want to have filtered data

opt_method = 1; % See below
%change it to <2> for Gauss-Newton least squares search
%change it to <3> for Levenberg-Marquardt least squares search 
%change it to <4> for Steepest descent least squares search
%change it to <1> for the default process if you see errors

nx=4;  % number of the state variables


output_sys_m = sysid(sampleTime, filter_wex_test, opt_method, nx, data4id, data4vl);




%% Functions

function innova_sys_m = sysid(sampleTime, filter_wex_test, opt_method, nx, data4id, data4vl)
%%
Comp_freq    = data4id.Comp_freq;
Tw_wex_out   = data4id.Tw_wex_out;    
Tw_wex_in    = data4id.Tw_wex_in;     
Tw_wex_delta = data4id.Tw_wex_delta; 
Air_Temp     = data4id.Air_Temp;     
TR_wex_out   = data4id.TR_wex_out;  

Comp_freq2   = data4vl.Comp_freq2;     
Tw_wex_out2  = data4vl.Tw_wex_out2;    
Tw_wex_in2   = data4vl.Tw_wex_in2;  
Tw_wex_delta2= data4vl.Tw_wex_delta2;  
TR_wex_out2  = data4vl.TR_wex_out2;    
Air_Temp2    = data4vl.Air_Temp2;        

%% Internal Fnc Params


diftemp_init = 0; % simulink system intitial point




%% Data Filtering

for i=1:length(Comp_freq)
   if Comp_freq(i)==0
       Tw_wex_delta(i) = 0;
       Tw_wex_out(i) = Tw_wex_in(i);
   end
end


%new filter to be tested...
if filter_wex_test == 1
      for i=1:length(Comp_freq2)
          if Comp_freq2(i)==0
              Tw_wex_delta2(i) = 0;
              Tw_wex_out2(i) = Tw_wex_in2(i);
          end
      end
end




% 
% for i=1:length(Tw_wex_delta)
%     if Tw_wex_delta(i)<0
%        Tw_wex_delta(i) = 0;
%        Tw_wex_out(i) = Tw_wex_in(i);
%     end
% end



%%   SYSId


hvac_a = iddata(Tw_wex_out, [Comp_freq , Tw_wex_in, Air_Temp], sampleTime);
hvac_a.OutputName = {'Tw_wex_out'};
hvac_a.InputName  = {'Comp Freq' , 'Tw_wex_in', 'Air temp'};



%hvac_21 = iddata([MT4,MT5], Comp, 0.05);
%hvac_21.OutputName = {'temp_MT4';'temp_MT5'};
%hvac_21.InputName  = {'comp'};
%figure()
%plot(hvac_21(:,1,1));
%figure()
%plot(hvac_21(:,2,1));
% mpf = ssest(hvac_a(1:60));
% hvac_23 = iddata([MT4,MT5], [Comp, MT6, PT6], 0.05);
% hvac_23.OutputName = {'temp_MT4';'temp_MT5'};
% hvac_23.InputName  = {'comp', 'temp_MT6','temp_PT6'};
% np=5; % Estimate a fifth-order state-space model
% mp5 = ssest(hvac_23(1:60),np);
% np=2; % Estimate a second-order state-space model
% mp2 = ssest(hvac_23(1:60),np);


mpf = fact(hvac_a, sampleTime, nx, opt_method);



%% Plot

figure(1)
plot(Comp_freq, Tw_wex_delta, '.' );
xlabel('Comp Freq'); ylabel('Water Temp WEx Delta')
grid on
%hold on


%% Primary compare





% b1= round(length(Comp_freq)/3);
% b2= round(length(Comp_freq)/2);
%    
% 
% Comp_freq_smpl = Comp_freq(b1:b2);
% Tw_wex_in_smpl = Tw_wex_in(b1:b2);
% Air_Temp_smpl  = Air_Temp(b1:b2);
% dim=size(Comp_freq_smpl); dim1=dim(1); 
% 
% timevec = sampleTime *(0:dim1-1);
% exmpl.time = timevec';
% 
% 
% exmpl.signals.values =   [Comp_freq_smpl , Tw_wex_in_smpl, Air_Temp_smpl]; %[Comp_freq , Tw_wex_in]
% exmpl.signals.dimensions =   3;
% 
% stop_time = length(Comp_freq_smpl) * sampleTime;
% 
% 
% 
% %assignin('base', 'sim_time', stop_time);
% %assignin('base', 'sim_data', exmpl);
% 
% 
% 
% 
% 
% open('modelsim.slx');
% simout = sim('modelsim.slx', 'ReturnWorkspaceOutputs','on');
% 
% out = simout ;
% 
% 
% %%%%%
% 
% 
% figure(3)
% plot(out.diftempxmpl, 'r.' , 'LineWidth',2,  'MarkerSize',10 );
% hold on
% 
% 
% plot(Tw_wex_out(b1:b2), 'k.' , 'LineWidth',2,  'MarkerSize',10 );
% xlabel('t (min)'); ylabel('Water Temp WEx outlet')
% legend('Identified','data')





%%  The Compare


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a1=1; a2= length(Comp_freq2)-1;


Comp_freq_smpl2 = Comp_freq2(a1:a2);
Tw_wex_in_smpl2 = Tw_wex_in2(a1:a2);
Air_Temp_smpl2   = Air_Temp2(a1:a2);

dim=size(Comp_freq_smpl2); dim1=dim(1); 

timevec = sampleTime*(0:dim1-1);
exmpl.time = timevec';

exmpl.signals.values =   [Comp_freq_smpl2 , Tw_wex_in_smpl2, Air_Temp_smpl2];
exmpl.signals.dimensions =   3;


% diftemp_init =  Tw_wex_out2(a1);

stop_time= length(Comp_freq_smpl2) * sampleTime;


assignin('base', 'stop_time', stop_time)
assignin('base', 'exmpl', exmpl);
assignin('base', 'mpf', mpf);
assignin('base', 'sampleTime', sampleTime);


open('modelsim.slx');
simout = sim('modelsim.slx', 'ReturnWorkspaceOutputs','on');

out = simout;

% figure(6)
% plot(out.freqxmpl, out.diftempxmpl, '.' );

%close all %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7)

set(gcf,'color','w');

plot(out.diftempxmpl, 'r.', 'LineWidth',2,  'MarkerSize',10 );
hold on

plot(Tw_wex_out2(a1:a2), 'k.' , 'LineWidth',2,  'MarkerSize',10  );
xlabel('t (min)'); ylabel('Water Temp WEx outlet')
legend('Identified','data')
axis([700  1900 -5 70])
grid on



%%  Sys results

Fianl_pred_err      = mpf.Report.Fit.FPE; 
Mean_sq_err         = mpf.Report.Fit.MSE;
FitPercent          = mpf.Report.Fit.FitPercent;
mpf_NoiseVariance   = mpf.NoiseVariance;
%mpf.Report.InitialState = 'zero';
mpf_InputDelay      = mpf.InputDelay;
mpf_OutputDelay     = mpf.OutputDelay;
mpf.TimeUnit        = 'min';
eig(mpf.A);


%% Controllability check

if rank(ctrb(mpf.A,mpf.B))== nx
    disp('The obtained system is Controllable')
    
end


%% Observability check

if rank(obsv(mpf.A,mpf.C))== nx
    disp('The obtained system is Observable')    
end  


%%   Command window display
disp ('Final Prediction Error is');
disp(Fianl_pred_err);
disp('');
disp ('Mean squared error is');
disp(Mean_sq_err);
disp('');
disp ('Curve Fitting Percentage is');
disp(FitPercent);
disp('');


%% output of function


innova_sys_m.a = out.diftempxmpl; 
innova_sys_m.b = Tw_wex_out2(a1:a2);
innova_sys_m.Fianl_pred_err = Fianl_pred_err;
innova_sys_m.Mean_sq_err    = Mean_sq_err;
innova_sys_m.FitPercent     = FitPercent;   



end




function mpf = fact( hvac_a, sampleTime, nx, opt_method)
    
switch opt_method
    case 1
opt = ssestOptions('InitializeMethod' , 'n4sid', 'InitialState', 'backcast', ...
      'N4Weight', 'auto', 'Focus', 'prediction', 'EnforceStability', false, ...
      'EstimateCovariance', true, 'display', 'off', 'InputOffset', [],    ...
      'OutputOffset', [], 'OutputWeight', [], 'SearchMethod', 'auto');
    case 2 %  Gauss-Newton least squares search
opt = ssestOptions('InitializeMethod' , 'n4sid', 'InitialState', 'backcast', ...
      'N4Weight', 'auto', 'Focus', 'prediction', 'EnforceStability', false, ...
      'EstimateCovariance', true, 'InputOffset', [],   'OutputOffset',  ...
      [], 'OutputWeight', [], 'SearchMethod', 'gn');        
    case 3  % Levenberg-Marquardt least squares search
opt = ssestOptions('InitializeMethod' , 'n4sid', 'InitialState', 'backcast', ...
      'N4Weight', 'auto', 'Focus', 'prediction', 'EnforceStability', false, ...
      'EstimateCovariance', true, 'InputOffset', [], 'display', 'off',  ...
      'OutputOffset', [], 'OutputWeight', [], 'SearchMethod', 'lm');        
    case 4 % Steepest descent least squares search
opt = ssestOptions('InitialState', 'backcast', 'Focus', 'prediction', ...
      'InitializeMethod' , 'n4sid', 'N4Weight', 'auto', 'EnforceStability', false, ...
      'EstimateCovariance', true, 'display', 'off', 'InputOffset', [],    ...
      'OutputOffset', [], 'OutputWeight', [], 'SearchMethod', 'grad');        
end

  
%'backcast' 'zero'

InM = 'auto' ;          % InM = 'n4sid';
N4W = 'auto' ;          % N4W = 'SSARX';
SearchMethod = 'auto';  
% SearchMethod = 'gn';    Gauss-Newton least squares
% SearchMethod = 'lm';    Levenberg-Marquardt least squares
% SearchMethod = 'grad';  Steepest descent least squares 

mpf_option = ssestOptions('InitializeMethod' , InM , 'InitialState', 'backcast', ...
      'N4Weight', N4W, 'Focus', 'prediction', 'EnforceStability', false, ...
      'EstimateCovariance', true, 'display', 'off', 'InputOffset', [],    ...
      'OutputOffset', [], 'OutputWeight', [], 'SearchMethod', SearchMethod);
mpf_option.Regularization.Lambda;
  
option = ssestOptions('Focus', 'prediction', 'InitializeMethod' , InM , ...
         'InitialState', 'backcast', 'EstimateCovariance', true);
option.Regularization.Lambda;

% opt.Regularization.Lambda = .1;

sys_option = ssestOptions();
sys_option.Regularization.Lambda                         = 0;     % 10 the bias versus variance tradeoff
sys_option.Regularization.R                              = 1;     %  Weighting matrix
%sys_option.Regularization.Nominal                       = 0;     %   
sys_option.SearchOptions.Tolerance                       = 0.01;  % Minimum percentage difference ... 
%between the current value of the loss function and its expected improvement after the next iteration
sys_option.SearchOptions.MaxIterations                   = 20;    % Maximum number of 
%iterations during loss-function minimization
sys_option.SearchOptions.Advanced.GnPinvConstant         = 10000; % Jacobian matrix singular value threshold
sys_option.SearchOptions.Advanced.InitialGnaTolerance    = 0.0001;% Initial value of gamma 
sys_option.SearchOptions.Advanced.LMStartValue           = 0.001; %	Starting value of search-direction 
%length d in the Levenberg-Marquardt method
sys_option.SearchOptions.Advanced.LMStep                 = 2;     % Size of the Levenberg-Marquardt step
sys_option.SearchOptions.Advanced.MaxBisections          = 25;    % Max number of bisections used for line search
sys_option.SearchOptions.Advanced.MaxFunctionEvaluations = Inf;   %Maximum number of calls to the model file
sys_option.SearchOptions.Advanced.MinParameterChange     = 0;     % Smallest parameter update allowed per iteration
sys_option.SearchOptions.Advanced.RelativeImprovement    = 0;     % Relative improvement threshold
sys_option.SearchOptions.Advanced.StepReduction          = 2;     % Step reduction factor
mpf = ssest(hvac_a, nx ,'Ts', sampleTime, opt);


end

