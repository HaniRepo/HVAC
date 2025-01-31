function innova_sys_m = sysid(sampleTime, filter_wex_test, opt_method, nx, data4id, data4vl,m,n)
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


figure(10 * m + n)
%figure(7)
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


innova_sys_m.a =out.diftempxmpl;    % Comes from simulink
innova_sys_m.b = Tw_wex_out2(a1:a2); % comes from Data
innova_sys_m.Fianl_pred_err = Fianl_pred_err;
innova_sys_m.Mean_sq_err    = Mean_sq_err;
innova_sys_m.FitPercent     = FitPercent;   


end
