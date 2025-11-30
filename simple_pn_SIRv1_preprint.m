 % This is the main PN SIRS Model code for

% This code will not work if you do not 
% first run code for simple_pn_SIRv1_pdf.m and COMMON_PRE.m 

%% Descrtiption and Directions
% First, choose section to run by changing Run_Section on line 54 to desired number.
% Then go to section and change variables listed in description to desired
% level. With Run_Section in find, hitting command+G on mac will skip between each section.
% 1 = run basic SIRS PN model with 1 set of parameters and ode set to same parameters
    % beta, gamma, delta, tau (called timidivi, time steps per unit time),
    % Initial Conditions, a total time for simulations are what are params
    % to be set
    % go to line 53 to change parameters and look over section

% 2 = run SIRS PN model over range of parameters and compares RRMSE
    % beta, gamma can be any size grid between [0,1] (using bglength and spacer,
    % delta needs to be grid of size 5, tau (called timidivi, time steps per unit time),
    % Initial Conditions, a total time for simulations are params to be set
    % Go to line __ to change parameters and look over section

% 3 = finds RRMSE for various spread of parameter values
    % Set tau_spread, where tau is spread of parameters between [0,1]
    % different from how tau is defined in the paper
    % Go to line 797 to change parameters and look over section

% 4 = comparing different types of rounding methods for single param value
    % all params to be set are same as section 1
    % Go to line 2093 to change parameters and look over section

% 5 = run SIRS PN model over range of parameters and Rounding methods and compares RRMSE
    % all params set same as Run_Section 2
    % Go to line __ to change parameters and look over section

% 6 = Load parameters from supercomputer runs to compare mean RRMSE of
    % different time steps
    % Go to line 4569 to change parameters and look over section

% 7 = Computation time plot
    % Go to line 4627 to change parameters and look over section

% 8 = Comparing outside PN runs

% 9 =

% 10= parellel runs manual recombination (used for popscalar runes
%   Go to line 5400


clear all; clc; 
global global_info m0Susceptible m0Infected m0Recovered %beta
tic %for timing program run time

Run_Section=1;

%% Run_Section 1, run 1 SIRS PN model and compare to ODE
if Run_Section==1
clear all; clc; 
global global_info m0Susceptible m0Infected m0Recovered %beta
tic

% TL=readtable('exampleStochastic_Branching_1direct_0.1_0.1_0.001.csv');
% Spike_S=table2array(TL(1:100,6));
% Spike_I=table2array(TL(1:100,2));
% Spike_R=table2array(TL(1:100,4));

% popscaler_ar=1:1:10;
% timidivi_ar=[1, 10:10:80];
%test
popscaler_ar=1:1:2;
timidivi_ar=[1, 2];

% %1
% popscaler_ar=1:1:5;
% timidivi_ar=[1, 10:10:80];
% %2
% popscaler_ar=6;
% timidivi_ar=[1, 10:10:80];
% %3
% popscaler_ar=7;
% timidivi_ar=[1, 10:10:80];
% %4
% popscaler_ar=8;
% timidivi_ar=[1, 10:10:80];
% %5
% popscaler_ar=9;
% timidivi_ar=[1, 10:10:80];
% %6
% popscaler_ar=10;
% timidivi_ar=[1, 10:10:80];

err_S=zeros(length(popscaler_ar), length(timidivi_ar));
err_I=zeros(length(popscaler_ar), length(timidivi_ar));
err_R=zeros(length(popscaler_ar), length(timidivi_ar));
for popscalerfor=1:1:length(popscaler_ar)
    for timidivifor=1:1:length(timidivi_ar)
%popscaler=1;
%timidivi=1;                    %tau in paper
popscaler=popscaler_ar(popscalerfor);
timidivi=timidivi_ar(timidivifor);
delta=0.001; beta=0.1; gamma=0.1;       %parameter values 
%delta=0.1; beta=0.2; gamma=0.3;


try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
% in base settings
MAX_ITERATIONS = 100*timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

% Initial Population
S_0=1000*popscaler;
I_0=10*popscaler;
R_0=10*popscaler;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        global_info.tInf_pInf = zero;

        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
        % pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-ceil(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        % tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-ceil(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
   
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        % global_info.pSus_tInf = ceil(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        % global_info.tInf_pInf = ceil(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);

        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
        % pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-ceil(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        % tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-ceil(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
   
    else 
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        % global_info.pSus_tInf = ceil(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        % global_info.tInf_pInf = ceil(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);

        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        % pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-ceil(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        % tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-ceil(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    % global_info.pInf_tRec = ceil(gamma_T * (m0Infected) + pInf_tRec_resEr);
    % global_info.tRec_pRec = ceil(gamma_T * (m0Infected) + tRec_pRec_resEr);
    % global_info.pRec_tSus = ceil(delta_T * (m0Recovered) + pRec_tSus_resEr);
    % global_info.tSus_pSus = ceil(delta_T * (m0Recovered) + tSus_pSus_resEr);
    
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    % pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-ceil(gamma_T * (m0Infected) + pInf_tRec_resEr);
    % tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-ceil(gamma_T * (m0Infected) + tRec_pRec_resEr);
    % pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-ceil(delta_T * (m0Recovered) + pRec_tSus_resEr);
    % tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-ceil(delta_T * (m0Recovered) + tSus_pSus_resEr);


    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', int2str(delta), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

MAX_ITERATIONS = 100*timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

%Initial Population
S_0=1000*popscaler;
I_0=10*popscaler;
R_0=10*popscaler;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop  
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        %beta2=(1/(m0Infected+??????));
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end



y0=[S_0/popscaler;I_0/popscaler;R_0/popscaler];           %ICs for ODE
delta_ODE=delta; 
beta_ODE=beta; 
gamma_ODE=gamma;
params=[delta_ODE;beta_ODE;gamma_ODE];
[t,y]=ode15s(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);    %ODE model run
%[t,y]=ode89(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);

% RRMSE Calculation
 err_S(popscalerfor,timidivifor) = (rmse(MSusceptible(1:timidivi:MAX_ITERATIONS)/popscaler , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_I(popscalerfor,timidivifor) = (rmse(MInfected(1:timidivi:MAX_ITERATIONS)/popscaler , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 err_R(popscalerfor,timidivifor) = (rmse(MRecovered(1:timidivi:MAX_ITERATIONS)/popscaler , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;
    end
end

save('RRMSE_heat_pop_tau_1', 'err_S', 'err_I', 'err_R');
%% --- VISUALIZATION SECTION ---

% Define plotting parameters
% Setup Axis Limits (Use your hardcoded values or calculate dynamic max)
casxismin = 0;
% Calculate max error across all matrices to set a common scale (or use fixed 44)
% global_max_err = max([max(err_S(:)), max(err_I(:)), max(err_R(:))]);
% casxismax = ceil(global_max_err); 
 casxismax = 44; % Uncomment to use your fixed limit

% Define Custom Colormap
Cus_map=[(ones(2,3).*[0 0.5 0]);(ones(8,3).*[0 1 0]);(ones(40,3).*[1 1 0]);(ones(50,3).*[1 0.75 0]);...
    (ones(100,3).*[1 0.4 0]);(ones(240,3).*[1 0 0])];

% Create Figure
fig = figure('Name', 'RRMSE Analysis', 'Color', 'w', 'Position', [100, 100, 1200, 400]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'RRMSE Comparison for SIR Populations', 'FontSize', 16, 'FontWeight', 'bold');

% Create Grid for plotting
% X axis = Time Division, Y axis = Pop Scaler
[X, Y] = meshgrid(timidivi_ar, popscaler_ar);

% --- Plot 1: Susceptible Error ---
nexttile
s1 = surf(X, Y, err_S);
title('Susceptible Error (err\_S)')
ylabel('Pop Scaler'); 
xlabel('Time Steps'); 
zlabel('RRMSE (%)');
s1.FaceColor = 'interp';
s1.EdgeColor = 'interp'; % or 'none' for smoother look
view(2); % Top-down heatmap view
colormap(Cus_map);
clim([casxismin, casxismax]);
grid off; box on;

% --- Plot 2: Infected Error ---
nexttile
s2 = surf(X, Y, err_I);
title('Infected Error (err\_I)')
ylabel('Pop Scaler'); 
xlabel('Time Steps'); 
zlabel('RRMSE (%)');
s2.FaceColor = 'interp';
s2.EdgeColor = 'interp';
view(2);
colormap(Cus_map);
clim([casxismin, casxismax]);
grid off; box on;

% --- Plot 3: Recovered Error ---
nexttile
s3 = surf(X, Y, err_R);
title('Recovered Error (err\_R)')
ylabel('Pop Scaler'); 
xlabel('Time Steps'); 
zlabel('RRMSE (%)');
s3.FaceColor = 'interp';
s3.EdgeColor = 'interp';
view(2);
colormap(Cus_map);
clim([casxismin, casxismax]);
grid off; box on;

% Add shared colorbar
h = colorbar;
h.Layout.Tile = 'east';
h.Label.String = "RRMSE (%)";
h.Label.Rotation = 270;
h.Label.VerticalAlignment = "bottom";
% % Plot Results
% figure()
% plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible/popscaler , 'b--o','LineWidth',8);
% hold on
% plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MInfected/popscaler ,'c--o','LineWidth',8);
% plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MRecovered/popscaler , 'r--o','LineWidth',8);
% 
% plot(t,y(:,1), 'y','LineWidth',5)
% plot(t,y(:,2),'m','LineWidth',5)
% plot(t,y(:,3), 'k','LineWidth',5);
% plot(t,Spike_S, 'r','LineWidth',2)
% plot(t,Spike_I,'g','LineWidth',2)
% plot(t,Spike_R, 'b','LineWidth',2)
% legend('Susceptible_{PN}','Infected_{PN}','Recovered_{PN}',...
%     'Susceptible_{ODE}','Infected_{ODE}','Recovered_{ODE}', ...
%     'Susceptible_{Sk}','Infected_{Sk}','Recovered_{Sk}');
% title("SIR Petri Net vs ODE", 'FontSize', 24)
% subtitle("beta="+beta+', gamma='+gamma+' delta='+delta+", Susc. MSE="+err_S...
%    +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
% ylabel('Tokens (population)'); xlabel('Time'); 
% hold off
% toc

%end of section 1



















%% Run Section 2, runs PN and ODE model over range of params
elseif Run_Section==2

timidivi_col_S=[];
timidivi_col_I=[];
timidivi_col_R=[];

timidivi=1;                        %tau in paper

bglength=10;                        %size of grid for beta and gamma

spacer=linspace(0,1,bglength);      %grid of params for beta and gamma

MSE_s=[zeros(bglength,bglength,5)]; %storage array for place s (susceptible) RRMSE
MSE_i=[zeros(bglength,bglength,5)]; %storage array for place s (susceptible) RRMSE
MSE_r=[zeros(bglength,bglength,5)]; %storage array for place s (susceptible) RRMSE

 deltactr=1;
  for delta=[0,logspace(-3,0,4)]    %grid of params for delta (must be of size 5 for figs)
     betactr=1;
     for beta=spacer   
         gammactr=1;
         for gamma=spacer

try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;   %scale params based on tau

MAX_ITERATIONS = 100*timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
   
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
   
    %Find arc weights for each iteration using standard+residual rounding
    %method
    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

       
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

MAX_ITERATIONS = 100*timidivi;

%Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=0.1','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected \delta=0.1','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=0.1','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=1','fontsize',12)
s.FaceColor='interp';
end
end

%ODE function run
y0=[S_0;I_0;R_0];
delta_ODE=delta; 
beta_ODE=beta; 
gamma_ODE=gamma;
params=[delta_ODE;beta_ODE;gamma_ODE];
[t,y]=ode15s(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);


% RRMSE Calculation
 err_S = (rmse(MSusceptible(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_I = (rmse(MInfected(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 err_R = (rmse(MRecovered(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;
 % % 
 % % timidivi_col_S=[timidivi_col_S err_S];
 % % timidivi_col_I=[timidivi_col_I err_I];
 % % timidivi_col_R=[timidivi_col_R err_R];
% %end

disp(['still going: ', num2str(deltactr)]);

        % MSE_s(gammactr,betactr)=err_S;
        % MSE_i(gammactr,betactr)=err_I;
        % MSE_iBeta(gammactr,betactr)=beta;
        % MSE_iGamma(gammactr,betactr)=gamma;
        % MSE_r(gammactr,betactr)=err_R;

        %RRMSE storage
        MSE_s(gammactr,betactr,deltactr)=err_S;
        MSE_i(gammactr,betactr,deltactr)=err_I;
        MSE_r(gammactr,betactr,deltactr)=err_R;

           gammactr=gammactr+1;
         end
        betactr=betactr+1;
     end
     deltactr=deltactr+1;
  end

%Plot the results of grid param value RRMSE vs ODE
% casxismin=min([min(MSE_s(:)),min(MSE_i(:)),min(MSE_r(:))]);
% trial=1;
% if trial==1
% casxismax=max([max(MSE_s(:)),max(MSE_i(:)),max(MSE_r(:))]);
% trial=trial+1;
% else
% end

casxismin=0;
%Max RRMSE was found across all runs and then set
casxismax=44; 

%custom color map
Cus_map=[(ones(2,3).*[0 0.5 0]);(ones(8,3).*[0 1 0]);(ones(40,3).*[1 1 0]);(ones(50,3).*[1 0.75 0]);...
    (ones(100,3).*[1 0.4 0]);(ones(240,3).*[1 0 0])];


figure()
%title("RRMSE for parameter Ranges with 1 pn step per 1 ode unit time", 'FontSize', 14)
tiledlayout(5,3);
nexttile
s=surf(spacer,spacer,MSE_s(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.001','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.001','fontsize',12')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.001','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,3));
ylabel('\gamma (rate of recovery per unit time)','fontsize',24); 
%xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.01','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.01','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.01','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,5));
%ylabel('gamma'); zlabel('RRMSE');
xlabel('\beta (rate of infection per unit of time)','FontSize',24); 
title('Infected \delta=1','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=1','fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');
%h=axes('Position', [0.85 0.15 0.05 0.7]);% axes('Position', [0.85 0.15 0.05 0.7]),   %'Location', 'westoutside',
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;

MSE_s
MSE_i
MSE_r
%save all variables
save(strcat(num2str(timidivi),'paramsb.mat'))

%end of section 2



















%% Run Section 3
%not currently used in paper so code not updated
elseif Run_Section==3
timidivi_col_S=[];
timidivi_col_I=[];
timidivi_col_R=[];
%timidivi_timestep=1:6:79;
%for timidivi=1:6:79

%timidivi=1;
%for timidivi=[2]
timidivi=20;
% The following section is used for finding the RRMSE for various
%parameter values in comparison to the respective ODE, this is limited two
%comparison of 3 parameters in this configuration, for visual display.
  bglength=1;
% tau_S=[];
% tau_I=[];
% tau_R=[];
% tau_spread=linspace(0,1,10);

% for tau=linspace(0,1,10)
%spacer=[0 logspace(-5,0,bglength-1)];
%spacer=[0 logspace(-5,-1,bglength) linspace(0.2,0.8,bglength) (ones(1,bglength)-logspace(-1,-5,bglength)) 1];
% spacer=linspace(0,tau,bglength);
spacer=linspace(0,1,bglength);
% gridsize=size(spacer,2);
% 
% MSE_s=[zeros(gridsize,gridsize)];
% MSE_i=[zeros(gridsize,gridsize)];
% MSE_r=[zeros(gridsize,gridsize)];
MSE_s=[zeros(bglength,bglength,5)];
MSE_i=[zeros(bglength,bglength,5)];
MSE_r=[zeros(bglength,bglength,5)];

% MSE_S_c=zeros(1,5);
% MSE_S_r=zeros(1,5);
% MSE_S_v=zeros(1,5);
% MSE_S_C=zeros(bglength*bglength,5);
% %0.1
%deltar=[0,logspace(-3,0,4)];
 % % % % % deltactr=1;
 % % % % %  for delta=[0,logspace(-3,0,4)]
 % % % % %     betactr=1;
 % % % % %     %for beta=linspace(beta_ODE*0.9,beta_ODE*1.1,bglength) %beta_ODE = .0008;     gamma_ODE = 0.08;
 % % % % %     for beta=spacer   
 % % % % %         %beta=spacer(betaindex);
 % % % % %         gammactr=1;
 % % % % %         %for gamma=linspace(gamma_ODE*0.9,gamma_ODE*1.1,bglength)
 % % % % %         for gamma=spacer
             %gamma=spacer(gammaindex);
delta=1; beta=1; gamma=1;
% delta1=1/timidivi; beta1=1/timidivi; gamma1=1/timidivi;
for Rounder_test=1:1:4
    if Rounder_test==1
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;
%beta = 0.1;     gamma = 0.2;    delta = 0.1;
%beta = 0.2;     gamma = 0.5;    delta = 0.4;
%beta = 3.1000e-04;     gamma =  0.0809;    
%delta = 0.00;


%Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
    %  Arc weights
    %This below is for splitting time
    % if Number_of_iterations==1
    %     beta=beta1;
    %     delta=delta1;
    %     gamma=gamma1;
    % else
    %     % beta=1;
    %     % delta=1;
    %     % gamma=1;
    % end
    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        %beta2=(1/(m0Infected+??????));
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    %{ 
%One way to handle rounding error
    if m0Infected<=(1/beta)

    if round(beta * m0Susceptible * m0Infected)==0 && (beta * m0Susceptible * m0Infected)~=0 && MinfiringS==0 
        timestepstofireS=round((1/(beta * m0Susceptible * m0Infected)));
        tstf_counterS=tstf_counterS+1;
        MinfiringS=1;
        global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
        (beta * m0Susceptible * m0Infected)
    elseif MinfiringS==0
        global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
        (beta * m0Susceptible * m0Infected)
    elseif MinfiringS==1 && tstf_counterS<timestepstofireS
        tstf_counterS=tstf_counterS+1;
        global_info.pSus_tInf = 0;
        global_info.tInf_pInf = 0;
        (beta * m0Susceptible * m0Infected)
    elseif MinfiringS==1 && tstf_counterS==timestepstofireS
        global_info.pSus_tInf = 1;
        global_info.tInf_pInf = 1;
        (beta * m0Susceptible * m0Infected)
        MinfiringS=0;
        tstf_counterS=0;
    else
        error("logic loop broken")
    end
 
    
    else

        beta2=1/m0Infected;
    if round(beta2 * m0Susceptible * m0Infected)==0 && (beta2 * m0Susceptible * m0Infected)~=0 && MinfiringS==0 
        timestepstofireS=round(1/(beta2 * m0Susceptible * m0Infected));
        tstf_counterS=tstf_counterS+1;
        MinfiringS=1;
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected);
    elseif MinfiringS==0
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected);
    elseif MinfiringS==1 && tstf_counterS<timestepstofireS
        tstf_counterS=tstf_counterS+1;
        global_info.pSus_tInf = 0;
        global_info.tInf_pInf = 0;
    elseif MinfiringS==1 && tstf_counterS==timestepstofireS
        global_info.pSus_tInf = 1;
        global_info.tInf_pInf = 1;
        MinfiringS=0;
        tstf_counterS=0;
    else
        error("logic loop broken")
    end

    end
   

    % global_info.pInf_tInf = ceil(beta * m0Susceptible * m0Infected);
    if round(gamma * m0Infected)==0 && (gamma * m0Infected)~=0 && MinfiringI==0 
        timestepstofireI=round(1/(gamma * m0Infected));
        tstf_counterI=tstf_counterI+1;
        MinfiringI=1;
        global_info.pInf_tRec = round(gamma * m0Infected);
        global_info.tRec_pRec = round(gamma * m0Infected);
    elseif MinfiringI==0
        global_info.pInf_tRec = round(gamma * m0Infected);
        global_info.tRec_pRec = round(gamma * m0Infected);
    elseif MinfiringI==1 && tstf_counterI<timestepstofireI
        tstf_counterI=tstf_counterI+1;
        global_info.pInf_tRec = 0;
        global_info.tRec_pRec = 0;
    elseif MinfiringI==1 && tstf_counterI==timestepstofireI
        global_info.pInf_tRec =1;
        global_info.tRec_pRec = 1;
        MinfiringI=0;
        tstf_counterI=0;
    else
        error("logic loop broken")
    end


    if round(delta * m0Recovered)==0 && (delta * m0Recovered)~=0 && MinfiringR==0 
        timestepstofireR=round(1/(delta * m0Recovered));
        tstf_counterR=tstf_counterR+1;
        MinfiringR=1;
        global_info.pRec_tSus = round(delta * m0Recovered);
        global_info.tSus_pSus = round(delta * m0Recovered);
    elseif MinfiringR==0
        global_info.pRec_tSus = round(delta * m0Recovered);
        global_info.tSus_pSus = round(delta * m0Recovered);
    elseif MinfiringR==1 && tstf_counterR<timestepstofireR
        tstf_counterR=tstf_counterR+1;
        global_info.pRec_tSus = 0;
        global_info.tSus_pSus = 0;
    elseif MinfiringR==1 && tstf_counterI==timestepstofireI
        global_info.pRec_tSus = 1;
        global_info.tSus_pSus = 1;
        MinfiringR=0;
        tstf_counterR=0;
    else
        error("logic loop broken")
    end
    %}


    % if floor(beta * m0Susceptible * m0Infected)>=0
    % global_info.pSus_tInf = floor(beta * m0Susceptible * m0Infected);
    % global_info.tInf_pInf = floor(beta * m0Susceptible * m0Infected);
    % else
    % global_info.pSus_tInf = 0;
    % global_info.tInf_pInf = 0;
    % end
    % if floor(gamma * m0Infected)>=0
    % global_info.pInf_tRec = floor(gamma * m0Infected);
    % global_info.tRec_pRec = floor(gamma * m0Infected);
    % else
    % global_info.pInf_tRec = 0;
    % global_info.tRec_pRec = 0;
    % end
    % if floor(delta * m0Recovered)>=0
    % global_info.pRec_tSus = floor(delta * m0Recovered);
    % global_info.tSus_pSus = floor(delta * m0Recovered);
    % else
    % global_info.pRec_tSus = 0;
    % global_info.tSus_pSus = 0;    
    % end
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        %beta2=(1/(m0Infected+??????));
        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end


MS_ER=MSusceptible;
MI_ER=MInfected;
MR_ER=MRecovered;


    elseif Rounder_test==2 %Floor
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;
%beta = 0.1;     gamma = 0.2;    delta = 0.1;
%beta = 0.2;     gamma = 0.5;    delta = 0.4;
%beta = 3.1000e-04;     gamma =  0.0809;    
%delta = 0.00;


% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 


    if floor(beta * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma * m0Infected);
    global_info.tRec_pRec = floor(gamma * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta * m0Recovered);
    global_info.tSus_pSus = floor(delta * m0Recovered);
    else
    global_info.pRec_tSus = 0;
    global_info.tSus_pSus = 0;    
    end
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch  % Floor Catch

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

%Initial Population
S_0=1000;
I_0=10;
R_0=10;
m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if floor(beta * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma * m0Infected);
    global_info.tRec_pRec = floor(gamma * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta * m0Recovered);
    global_info.tSus_pSus = floor(delta * m0Recovered);
    else
    global_info.pRec_tSus = 0;
    global_info.tSus_pSus = 0;    
    end

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'r--o','LineWidth',2);
MS_Floor=MSusceptible;
MI_Floor=MInfected;
MR_Floor=MRecovered;
    


%%CEILING
elseif Rounder_test==3 %Ceil
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;

% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
   


    global_info.pSus_tInf = ceil(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = ceil(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
   
    global_info.pRec_tSus = ceil(delta * m0Recovered);
    global_info.tSus_pSus = ceil(delta * m0Recovered);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop 
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

   global_info.pSus_tInf = ceil(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = ceil(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
   
    global_info.pRec_tSus = ceil(delta * m0Recovered);
    global_info.tSus_pSus = ceil(delta * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
MS_Ceil=MSusceptible;
MI_Ceil=MInfected;
MR_Ceil=MRecovered;




%%standar round
elseif Rounder_test==4 %Round standard
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
%delta=1; beta=1; gamma=1;

% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    delta = 0.001;

%% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered;

    global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = round(gamma * m0Infected);
    global_info.tRec_pRec = ceil(gamma * m0Infected);
   
   
    global_info.pRec_tSus = round(delta * m0Recovered);
    global_info.tSus_pSus = round(delta * m0Recovered);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

% Initial Population
S_0=1000;
I_0=10;
R_0=10;

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

   global_info.pSus_tInf = round(beta * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta * m0Susceptible * m0Infected);
   
   
    global_info.pInf_tRec = round(gamma * m0Infected);
    global_info.tRec_pRec = round(gamma * m0Infected);
   
   
    global_info.pRec_tSus = round(delta * m0Recovered);
    global_info.tSus_pSus = round(delta * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
MS_SR=MSusceptible;
MI_SR=MInfected;
MR_SR=MRecovered;
    end
end

y0=[S_0;I_0;R_0];
delta_ODE=delta; 
beta_ODE=beta; 
gamma_ODE=gamma;
params=[delta_ODE;beta_ODE;gamma_ODE];
[t,y]=ode15s(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);


figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'b--o','LineWidth',2);
hold on
plot(t,y(:,1), 'y','LineWidth',2)
hold off
% RRMSE
 err_S = (rmse(MSusceptible(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_I = (rmse(MInfected(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 err_R = (rmse(MRecovered(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;
 % % 
 % % timidivi_col_S=[timidivi_col_S err_S];
 % % timidivi_col_I=[timidivi_col_I err_I];
 % % timidivi_col_R=[timidivi_col_R err_R];
% %end


 % Plot times steps vs rrmse %%%%%%%%%%%%%%%%%%
% figure()
% plot(timidivi_timestep, timidivi_col_S, 'b--o','LineWidth',2);
% hold on
% plot(timidivi_timestep, timidivi_col_I,'c--o','LineWidth',2);
% plot(timidivi_timestep, timidivi_col_R, 'r--o','LineWidth',2);
% legend('Susceptible_{RRMSE}','Infected_{RRMS}','Recovered_{RRMSE}');
% title("SIR Petri Net vs ODE RRMSE based on times steps", 'FontSize', 24)
% subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi)
% xlabel('PN Time steps (per 1 ODE time interval)'); ylabel('RRMSE');
% hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_Ceil, 'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_ER, 'b','LineWidth',4);
plot(t,y(:,1), 'r','LineWidth',2); 
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Susceptible", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_Ceil,'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_ER, 'b','LineWidth',4);
plot(t,y(:,2), 'r','LineWidth',2);
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Infected", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_Ceil,'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_ER, 'b','LineWidth',4);
plot(t,y(:,3), 'r','LineWidth',2);  
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Recovered", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off
% % disp(['still going: ', num2str(deltactr)]);

% Plot the results of single param value run vs ODE %%%%%%%%%%%%%%%%%%
 figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'b--o','LineWidth',2);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MInfected,'c--o','LineWidth',2);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MRecovered, 'r--o','LineWidth',2);

plot(t,y(:,1), 'y','LineWidth',2)
plot(t,y(:,2),'m','LineWidth',2)
plot(t,y(:,3), 'k','LineWidth',2);
legend('Susceptible_{PN}','Infected_{PN}','Recovered_{PN}',...
    'Susceptible_{ODE}','Infected_{ODE}','Recovered_{ODE}');
title("SIR Petri Net vs ODE", 'FontSize', 24)
subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi+", Susc. MSE="+err_S...
   +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
hold off
% 
% figure()
% plot([1,1+1/timidivi,2:1:100], MSusceptible, 'b--o','LineWidth',2);
% hold on
% plot([1,1+1/timidivi,2:1:100], MInfected,'c--o','LineWidth',2);
% plot([1,1+1/timidivi,2:1:100], MRecovered, 'r--o','LineWidth',2);
% 
% plot(t,y(:,1), 'y','LineWidth',2)
% plot(t,y(:,2),'m','LineWidth',2)
% plot(t,y(:,3), 'k','LineWidth',2);
% legend('Susceptible_{PN}','Infected_{PN}','Recovered_{PN}',...
%     'Susceptible_{ODE}','Infected_{ODE}','Recovered_{ODE}');
% title("SIR Petri Net vs ODE", 'FontSize', 24)
% subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi+", Susc. MSE="+err_S...
%    +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
% hold off


        % MSE_s(gammactr,betactr)=err_S;
        % MSE_i(gammactr,betactr)=err_I;
        % MSE_iBeta(gammactr,betactr)=beta;
        % MSE_iGamma(gammactr,betactr)=gamma;
        % MSE_r(gammactr,betactr)=err_R;

        MSE_s(gammactr,betactr,deltactr)=err_S;
        MSE_i(gammactr,betactr,deltactr)=err_I;
        MSE_r(gammactr,betactr,deltactr)=err_R;
        % 
% MSE_S_c(deltactr)=zeros(1,5);
% MSE_S_r(deltactr)=zeros(1,5);
% MSE_S_v(deltactr)=zeros(1,5);
% MSE_S_C(dr)=[MSE_S_C(dr) err_S];
% % 
  % % % % %          gammactr=gammactr+1;
  % % % % %        end
  % % % % %       betactr=betactr+1;
  % % % % %    end
  % % % % %    deltactr=deltactr+1;
  % % % % % end
  
  
%   tau_S=[tau_S (sum(MSE_s,"all")/bglength^3)];
%   tau_I=[tau_I (sum(MSE_i,"all")/bglength^3)];
%   tau_R=[tau_R (sum(MSE_r,"all")/bglength^3)];
% end
% figure()
% plot(tau_spread, tau_S, 'b--o','LineWidth',2);
% hold on
% plot(tau_spread, tau_I, 'c--o','LineWidth',2);
% plot(tau_spread, tau_R, 'r--o','LineWidth',2);
% legend('Susceptible Mean RRMSE','Infected Mean RRMSE','Recovered Mean RRMSE');
% title("Mean Error for different Parameter ranges (0,\tau)", 'FontSize', 24)
% xlabel("\tau (parameters \gamma,\beta,\delta range [0,\tau])");ylabel('RRMSE Mean')
% subtitle("beta="+beta+', gamma='+gamma+' delta='+delta+", Susc. MSE="+err_S...
%     +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
%  hold off


% minValue = min(MSE_i(:));
% % Find all (row, column) pairs where M = the min value.
% [rows, columns] = find(MSE_i == minValue)
% % Print them out:
% for k = 1 : length(rows)
% 	fprintf('M equals the min value of %f at row = %d, column = %d.\n', ...
% 		minValue, rows(k), columns(k));
% end
% MSE_iBeta(rows(k),columns(k))
%  MSE_iGamma(rows(k),columns(k))

% Plot the results of grid param value RRMSE vs ODE
% casxismin=min([min(MSE_s(:)),min(MSE_i(:)),min(MSE_r(:))]);
% trial=1;
% if trial==1
% casxismax=max([max(MSE_s(:)),max(MSE_i(:)),max(MSE_r(:))]);
% trial=trial+1;
% else
% end

casxismin=0;
casxismax=44;
Cus_map=[(ones(2,3).*[0 0.5 0]);(ones(8,3).*[0 1 0]);(ones(40,3).*[1 1 0]);(ones(50,3).*[1 0.75 0]);...
    (ones(100,3).*[1 0.4 0]);(ones(240,3).*[1 0 0])];


figure()
%title("RRMSE for parameter Ranges with 1 pn step per 1 ode unit time", 'FontSize', 14)
tiledlayout(5,3);
nexttile
s=surf(spacer,spacer,MSE_s(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible Delta=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected Delta=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered Delta=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible Delta=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected Delta=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered Delta=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible Delta=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected Delta=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered Delta=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible Delta=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected Delta=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered Delta=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible Delta=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected Delta=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered Delta=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');
%h=axes('Position', [0.85 0.15 0.05 0.7]);% axes('Position', [0.85 0.15 0.05 0.7]),   %'Location', 'westoutside',
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;

MSE_s
MSE_i
MSE_r

%end

%end of Run Section 3














%% Run Section 4, comparing different types of rounding methods for single param value
elseif Run_Section==4

timidivi=20;                %tau in paper

  bglength=1;

spacer=linspace(0,1,bglength);

MSE_s=[zeros(bglength,bglength,5)];
MSE_i=[zeros(bglength,bglength,5)];
MSE_r=[zeros(bglength,bglength,5)];


delta=0.001; beta=0.1; gamma=0.9;   %parameter values

% Initial Population
S_0=1000;
I_0=10;
R_0=10;

MAX_ITERATIONS = 100*timidivi;  % Max iterations for Petri Net and ODE

for Rounder_test=1:1:4      %loop for each rounding methoded tested
    
    if Rounder_test==1      % Rounding method standard+residual
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop 
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
    
    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    
    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch % Rounding method standard+residual

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    
    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end

%store values of standard+residual run
MS_ER=MSusceptible;
MI_ER=MInfected;
MR_ER=MRecovered;




    elseif Rounder_test==2 % Rounding method Floor
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if floor(beta * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta_T * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma_T * m0Infected);
    global_info.tRec_pRec = floor(gamma_T * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta_T * m0Recovered);
    global_info.tSus_pSus = floor(delta_T * m0Recovered);
    else
    global_info.pRec_tSus_T = 0;
    global_info.tSus_pSus_T = 0;    
    end  

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch  % Rounding method Floor Catch

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
   
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if floor(beta * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta_T * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma_T * m0Infected);
    global_info.tRec_pRec = floor(gamma_T * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta_T * m0Recovered);
    global_info.tSus_pSus = floor(delta_T * m0Recovered);
    else
    global_info.pRec_tSus = 0;
    global_info.tSus_pSus = 0;    
    end

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
%plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'r--o','LineWidth',2);
%store values of floor function
MS_Floor=MSusceptible;
MI_Floor=MInfected;
MR_Floor=MRecovered;
    


%%CEILING
elseif Rounder_test==3 % Rounding method Ceiling

    try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    global_info.pSus_tInf = ceil(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta_T * m0Susceptible * m0Infected);
   
    global_info.pInf_tRec = ceil(gamma_T * m0Infected);
    global_info.tRec_pRec = ceil(gamma_T * m0Infected);
   
    global_info.pRec_tSus = ceil(delta_T * m0Recovered);
    global_info.tSus_pSus = ceil(delta_T * m0Recovered);
    
    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
    catch % Rounding method Ceiling

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop 
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

   global_info.pSus_tInf = ceil(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta_T * m0Susceptible * m0Infected);
  
    global_info.pInf_tRec = ceil(gamma_T * m0Infected);
    global_info.tRec_pRec = ceil(gamma_T * m0Infected);
   
    global_info.pRec_tSus = ceil(delta_T * m0Recovered);
    global_info.tSus_pSus = ceil(delta_T * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
    end
    %store value from ceiling function
MS_Ceil=MSusceptible;
MI_Ceil=MInfected;
MR_Ceil=MRecovered;






%standard round
elseif Rounder_test==4 % Rounding method standard

try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop  
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered;

    global_info.pSus_tInf = round(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta_T * m0Susceptible * m0Infected);
   
    global_info.pInf_tRec = round(gamma_T * m0Infected);
    global_info.tRec_pRec = ceil(gamma_T * m0Infected);
   
    global_info.pRec_tSus = round(delta_T * m0Recovered);
    global_info.tSus_pSus = round(delta_T * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch % Rounding method standard catch

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop 
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    global_info.pSus_tInf = round(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta_T * m0Susceptible * m0Infected);

    global_info.pInf_tRec = round(gamma_T * m0Infected);
    global_info.tRec_pRec = round(gamma_T * m0Infected);

    global_info.pRec_tSus = round(delta_T * m0Recovered);
    global_info.tSus_pSus = round(delta_T * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
%store standard rounding values
MS_SR=MSusceptible;
MI_SR=MInfected;
MR_SR=MRecovered;
    end
end

%Run ode model

y0=[S_0;I_0;R_0];
delta_ODE=delta; 
beta_ODE=beta; 
gamma_ODE=gamma;
params=[delta_ODE;beta_ODE;gamma_ODE];
[t,y]=ode15s(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);

% %% RRMSE
%  err_S = (rmse(MSusceptible(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
%  err_I = (rmse(MInfected(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
%  err_R = (rmse(MRecovered(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_Ceil, 'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MS_ER, 'b','LineWidth',4);
plot(t,y(:,1), 'r','LineWidth',2); 
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Susceptible", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_Ceil,'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MI_ER, 'b','LineWidth',4);
plot(t,y(:,2), 'r','LineWidth',2);
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Infected", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

figure()
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_Ceil,'color', [.15 .15 .15],'LineWidth',12);
hold on
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_Floor,'color', [.65 .65 .65],'LineWidth',8);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_SR, 'color', [.90 .90 .90],'LineWidth',4);
plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MR_ER, 'b','LineWidth',4);
plot(t,y(:,3), 'r','LineWidth',2);  
lgd=legend('Ceiling','Floor','Standard','Standard+Residuals','ODE');legend('Location','best');
lgd.Title.String = 'Arc Weight Rounding Method';
title("SIR PN Rounding Methods vs ODE Recovered", 'FontSize', 24)
subtitle("gamma (\gamma )="+gamma+", beta (\beta)="+beta+", delta (\delta)="+delta)
xlabel('Time')
ylabel('Population') 
xlim([0,105])
hold off

%end of Run Section 4



















%% Run Section 5, run SIRS PN model over range of parameters and Rounding methods and compares RRMSE
% not currently used in paper so code not updated
elseif Run_Section==5
timidivi_col_S=[];
timidivi_col_I=[];
timidivi_col_R=[];

timidivi=2;

bglength=10;

spacer=linspace(0,1,bglength);

MSE_s_ER=[zeros(bglength,bglength,5)];
MSE_i_ER=[zeros(bglength,bglength,5)];
MSE_r_ER=[zeros(bglength,bglength,5)];

MSE_s_Floor=[zeros(bglength,bglength,5)];
MSE_i_Floor=[zeros(bglength,bglength,5)];
MSE_r_Floor=[zeros(bglength,bglength,5)];

MSE_s_Ceil=[zeros(bglength,bglength,5)];
MSE_i_Ceil=[zeros(bglength,bglength,5)];
MSE_r_Ceil=[zeros(bglength,bglength,5)];

MSE_s_SR=[zeros(bglength,bglength,5)];
MSE_i_SR=[zeros(bglength,bglength,5)];
MSE_r_SR=[zeros(bglength,bglength,5)];

MSE_s_Ceil_R=[zeros(bglength,bglength,5)];
MSE_i_Ceil_R=[zeros(bglength,bglength,5)];
MSE_r_Ceil_R=[zeros(bglength,bglength,5)];

MAX_ITERATIONS = 100*timidivi;
S_0=1000;                       %initial token value of place S (initial population of susceptible)
I_0=10;                         %initial token value of place I (initial population of infected)
R_0=10;                         %initial token value of place R (initial population of recovered)

 deltactr=1;
  for delta=[0,logspace(-3,0,4)]
     betactr=1;
     %for beta=linspace(beta_ODE*0.9,beta_ODE*1.1,bglength) %beta_ODE = .0008;     gamma_ODE = 0.08;
     for beta=spacer   
         %beta=spacer(betaindex);
         gammactr=1;
         %for gamma=linspace(gamma_ODE*0.9,gamma_ODE*1.1,bglength)
         for gamma=spacer
             %gamma=spacer(gammaindex);

% delta1=1/timidivi; beta1=1/timidivi; gamma1=1/timidivi;
for Rounder_test=1:1:4      %loop for each rounding methoded tested
    
    if Rounder_test==1      % Rounding method standard+residual
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop 
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
    
    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    
    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch % Rounding method standard+residual

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if beta_T==0 || m0Infected==0
        global_info.pSus_tInf = zero;
        pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = zero;
        tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    elseif m0Infected<=(1/(beta_T))
        global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
        pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
        tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    else 

        beta2=1/(m0Infected+1);
        global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
        pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
        global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
        tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    
    end

    global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    global_info.pRec_tSus = round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * (m0Recovered) + pRec_tSus_resEr)-round(delta_T * (m0Recovered) + pRec_tSus_resEr);
    global_info.tSus_pSus = round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * (m0Recovered)+ tSus_pSus_resEr)-round(delta_T * (m0Recovered) + tSus_pSus_resEr);
    

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    
    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end

%store values of standard+residual run
MS_ER=MSusceptible;
MI_ER=MInfected;
MR_ER=MRecovered;


  
    
    elseif Rounder_test==2 % Rounding method Floor
try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if floor(beta_T * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta_T * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma_T * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma_T * m0Infected);
    global_info.tRec_pRec = floor(gamma_T * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta_T * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta_T * m0Recovered);
    global_info.tSus_pSus = floor(delta_T * m0Recovered);
    else
    global_info.pRec_tSus = 0;
    global_info.tSus_pSus = 0;    
    end  

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch  % Rounding method Floor Catch

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop   
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
   
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    if floor(beta_T * m0Susceptible * m0Infected)>=0
    global_info.pSus_tInf = floor(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = floor(beta_T * m0Susceptible * m0Infected);
    else
    global_info.pSus_tInf = 0;
    global_info.tInf_pInf = 0;
    end
    if floor(gamma_T * m0Infected)>=0
    global_info.pInf_tRec = floor(gamma_T * m0Infected);
    global_info.tRec_pRec = floor(gamma_T * m0Infected);
    else
    global_info.pInf_tRec = 0;
    global_info.tRec_pRec = 0;
    end
    if floor(delta_T * m0Recovered)>=0
    global_info.pRec_tSus = floor(delta_T * m0Recovered);
    global_info.tSus_pSus = floor(delta_T * m0Recovered);
    else
    global_info.pRec_tSus = 0;
    global_info.tSus_pSus = 0;    
    end

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
%plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'r--o','LineWidth',2);
%store values of floor function
MS_Floor=MSusceptible;
MI_Floor=MInfected;
MR_Floor=MRecovered;
    




%%CEILING
elseif Rounder_test==3 % Rounding method Ceiling

    try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    global_info.pSus_tInf = ceil(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta_T * m0Susceptible * m0Infected);
   
    global_info.pInf_tRec = ceil(gamma_T * m0Infected);
    global_info.tRec_pRec = ceil(gamma_T * m0Infected);
   
    global_info.pRec_tSus = ceil(delta_T * m0Recovered);
    global_info.tSus_pSus = ceil(delta_T * m0Recovered);
    
    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
    catch % Rounding method Ceiling

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop 
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

   global_info.pSus_tInf = ceil(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = ceil(beta_T * m0Susceptible * m0Infected);
  
    global_info.pInf_tRec = ceil(gamma_T * m0Infected);
    global_info.tRec_pRec = ceil(gamma_T * m0Infected);
   
    global_info.pRec_tSus = ceil(delta_T * m0Recovered);
    global_info.tSus_pSus = ceil(delta_T * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
    end
    %store value from ceiling function
MS_Ceil=MSusceptible;
MI_Ceil=MInfected;
MR_Ceil=MRecovered;







%standard round
elseif Rounder_test==4 % Rounding method standard

try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop  
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered;

    global_info.pSus_tInf = round(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta_T * m0Susceptible * m0Infected);
   
    global_info.pInf_tRec = round(gamma_T * m0Infected);
    global_info.tRec_pRec = ceil(gamma_T * m0Infected);
   
    global_info.pRec_tSus = round(delta_T * m0Recovered);
    global_info.tSus_pSus = round(delta_T * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
catch % Rounding method standard catch

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

% the loop 
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    global_info.pSus_tInf = round(beta_T * m0Susceptible * m0Infected);
    global_info.tInf_pInf = round(beta_T * m0Susceptible * m0Infected);

    global_info.pInf_tRec = round(gamma_T * m0Infected);
    global_info.tRec_pRec = round(gamma_T * m0Infected);

    global_info.pRec_tSus = round(delta_T * m0Recovered);
    global_info.tSus_pSus = round(delta_T * m0Recovered);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end
%store standard rounding values
MS_SR=MSusceptible;
MI_SR=MInfected;
MR_SR=MRecovered;




%%CEILING + Residual
elseif Rounder_test==5 % Rounding method Ceiling

    try
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    global_info.pSus_tInf = ceil(beta_T * m0Susceptible * m0Infected+pSus_tInf_resEr);
    pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-ceil(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
    global_info.tInf_pInf = ceil(beta_T * m0Susceptible * m0Infected+tInf_pInf_resEr);
    tInf_pInf_resEr=(beta_T * m0Susceptible * m0Infected+tInf_pInf_resEr)-ceil(beta_T * m0Susceptible * m0Infected+tInf_pInf_resEr);
   
    global_info.pInf_tRec = ceil(gamma_T * m0Infected + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * m0Infected+pInf_tRec_resEr)-ceil(gamma_T * m0Infected+pInf_tRec_resEr);
    global_info.tRec_pRec = ceil(gamma_T * m0Infected + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * m0Infected+ tRec_pRec_resEr)-ceil(gamma_T * m0Infected+ tRec_pRec_resEr);
   
    global_info.pRec_tSus = ceil(delta_T * m0Recovered+pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * m0Recovered+pRec_tSus_resEr)-ceil(delta_T * m0Recovered+pRec_tSus_resEr);
    global_info.tSus_pSus = ceil(delta_T * m0Recovered+tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * m0Recovered+tSus_pSus_resEr)-ceil(delta_T * m0Recovered+tSus_pSus_resEr);
    
    pns = pnstruct('simple_pn_SIRv1_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
    catch % Rounding method Ceiling

    disp(['Error at delta=', num2str(delta), ' , gamma=',num2str(gamma), ' , beta=',num2str(beta), ' , timesteps=',num2str(timidivi)]);
    
    if delta==0
    delta=delta; 
    else
    delta=delta-10^-7; 
    end
   
    if beta==0
    beta=beta; 
    else
    beta=beta-10^-7;
    end
    
    if gamma==0
    gamma=gamma;
    else
    gamma=gamma-10^-7;
    end
delta_T=delta/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 

m0Susceptible = S_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%the loop 
MinfiringS=0; tstf_counterS=0; MinfiringI=0; tstf_counterI=0; MinfiringR=0;tstf_counterR=0;

pSus_tInf_resEr=0; tInf_pInf_resEr=0;
pInf_tRec_resEr=0; tRec_pRec_resEr=0;
pRec_tSus_resEr=0; tSus_pSus_resEr=0;
zero=0;

Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', num2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    global_info.pSus_tInf = ceil(beta_T * m0Susceptible * m0Infected+pSus_tInf_resEr);
    pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-ceil(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
    global_info.tInf_pInf = ceil(beta_T * m0Susceptible * m0Infected+tInf_pInf_resEr);
    tInf_pInf_resEr=(beta_T * m0Susceptible * m0Infected+tInf_pInf_resEr)-ceil(beta_T * m0Susceptible * m0Infected+tInf_pInf_resEr);
   
    global_info.pInf_tRec = ceil(gamma_T * m0Infected + pInf_tRec_resEr);
    pInf_tRec_resEr=(gamma_T * m0Infected+pInf_tRec_resEr)-ceil(gamma_T * m0Infected+pInf_tRec_resEr);
    global_info.tRec_pRec = ceil(gamma_T * m0Infected + tRec_pRec_resEr);
    tRec_pRec_resEr=(gamma_T * m0Infected+ tRec_pRec_resEr)-ceil(gamma_T * m0Infected+ tRec_pRec_resEr);
   
    global_info.pRec_tSus = ceil(delta_T * m0Recovered+pRec_tSus_resEr);
    pRec_tSus_resEr=(delta_T * m0Recovered+pRec_tSus_resEr)-ceil(delta_T * m0Recovered+pRec_tSus_resEr);
    global_info.tSus_pSus = ceil(delta_T * m0Recovered+tSus_pSus_resEr);
    tSus_pSus_resEr=(delta_T * m0Recovered+tSus_pSus_resEr)-ceil(delta_T * m0Recovered+tSus_pSus_resEr);

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs
    %figure()
    %plotp(Sim_Results, {'pSusceptible','pInfected','pRecovered'});

    m0Susceptible = ntokens('pSusceptible');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
    end
    %store value from ceiling function
MS_Ceil_R=MSusceptible;
MI_Ceil_R=MInfected;
MR_Ceil_R=MRecovered;


    end
end

%ODE Model run
y0=[S_0;I_0;R_0];
delta_ODE=delta; 
beta_ODE=beta; 
gamma_ODE=gamma;
params=[delta_ODE;beta_ODE;gamma_ODE];
[t,y]=ode15s(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);

% RRMSE
%RRMSE for Standard+Residual found 
 err_S_ER = (rmse(MS_ER(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_I_ER = (rmse(MI_ER(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 err_R_ER = (rmse(MR_ER(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;

 %RRMSE for Floor function rounding found 
 err_S_Floor = (rmse(MS_Floor(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_I_Floor = (rmse(MI_Floor(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 err_R_Floor = (rmse(MR_Floor(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;

 %RRMSE for Ceiling function rounding found 
 err_S_Ceil = (rmse(MS_Ceil(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_I_Ceil = (rmse(MI_Ceil(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 err_R_Ceil = (rmse(MR_Ceil(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;

 %RRMSE for Standard rounding found 
 err_S_SR = (rmse(MS_SR(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_I_SR = (rmse(MI_SR(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 err_R_SR = (rmse(MR_SR(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;

 %RRMSE for Ceilin+Residual rounding found
 err_S_Ceil_R = (rmse(MS_Ceil_R(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_I_Ceil_R = (rmse(MI_Ceil_R(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 err_R_Ceil_R = (rmse(MR_Ceil_R(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;

disp(['still going: ', num2str(deltactr)]);

%RRMSE Stored
        MSE_s_ER(gammactr,betactr,deltactr)=err_S_ER;
        MSE_i_ER(gammactr,betactr,deltactr)=err_I_ER;
        MSE_r_ER(gammactr,betactr,deltactr)=err_R_ER;

        MSE_s_Floor(gammactr,betactr,deltactr)=err_S_Floor;
        MSE_i_Floor(gammactr,betactr,deltactr)=err_I_Floor;
        MSE_r_Floor(gammactr,betactr,deltactr)=err_R_Floor;

        MSE_s_Ceil(gammactr,betactr,deltactr)=err_S_Ceil;
        MSE_i_Ceil(gammactr,betactr,deltactr)=err_I_Ceil;
        MSE_r_Ceil(gammactr,betactr,deltactr)=err_R_Ceil;

        MSE_s_SR(gammactr,betactr,deltactr)=err_S_SR;
        MSE_i_SR(gammactr,betactr,deltactr)=err_I_SR;
        MSE_r_SR(gammactr,betactr,deltactr)=err_R_SR;
        
        MSE_s_Ceil_R(gammactr,betactr,deltactr)=err_S_Ceil_R;
        MSE_i_Ceil_R(gammactr,betactr,deltactr)=err_I_Ceil_R;
        MSE_r_Ceil_R(gammactr,betactr,deltactr)=err_R_Ceil_R;

 
           gammactr=gammactr+1;
         end
        betactr=betactr+1;
     end
     deltactr=deltactr+1;
  end
  




% Plot the results of grid param value RRMSE vs ODE
% casxismin=min([min(MSE_s(:)),min(MSE_i(:)),min(MSE_r(:))]);

 %casxismax=max([max(MSE_s_ER(:)),max(MSE_i_ER(:)),max(MSE_r_ER(:)),...
 %    max(MSE_s_Floor(:)),max(MSE_i_Floor(:)),max(MSE_r_Floor(:)),...
 %    max(MSE_s_Ceil(:)),max(MSE_i_Ceil(:)),max(MSE_r_Ceil(:)),...
 %    max(MSE_s_SR(:)),max(MSE_i_SR(:)),max(MSE_r_SR(:))]);

casxismin=0;
casxismax=44;          %use if max is known to better customize color map
Cus_map=[(ones(2,3).*[0 0.5 0]);(ones(8,3).*[0 1 0]);(ones(40,3).*[1 1 0]);(ones(50,3).*[1 0.75 0]);...
    (ones(100,3).*[1 0.4 0]);(ones(240,3).*[1 0 0])];

%Standard+Residuals Figure
figure()
tiledlayout(5,3);
nexttile
s=surf(spacer,spacer,MSE_s_ER(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('S+R Susceptible \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_ER(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_ER(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_ER(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_ER(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_ER(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_ER(:,:,3));
ylabel('\gamma (rate of recovery per unit of time)','Fontsize',24); 
%xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_ER(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_ER(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_ER(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_ER(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_ER(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_ER(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_ER(:,:,5));
%ylabel('gamma'); 
xlabel('\beta (rate of infection per unit of time)','Fontsize',24); %zlabel('RRMSE');
title('Infected \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_ER(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');
%h=axes('Position', [0.85 0.15 0.05 0.7]);% axes('Position', [0.85 0.15 0.05 0.7]),   %'Location', 'westoutside',
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;



%Floor figure
figure()
tiledlayout(5,3);
nexttile
s=surf(spacer,spacer,MSE_s_Floor(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Floor Susceptible \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Floor(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Floor(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Floor(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Floor(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Floor(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Floor(:,:,3));
ylabel('\gamma (rate of recovery per unit of time)','Fontsize',24); 
%xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Floor(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Floor(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Floor(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Floor(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Floor(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Floor(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Floor(:,:,5));
%ylabel('gamma'); 
xlabel('\beta (rate of infection per unit of time)','Fontsize',24); %zlabel('RRMSE');
title('Infected \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Floor(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');


%Ceiling Figure
figure()
tiledlayout(5,3);
nexttile
s=surf(spacer,spacer,MSE_s_Ceil(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Ceiling Susceptible \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Ceil(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Ceil(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Ceil(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Ceil(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Ceil(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Ceil(:,:,3));
ylabel('\gamma (rate of recovery per unit of time)','Fontsize',24); 
%xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Ceil(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Ceil(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Ceil(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Ceil(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Ceil(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Ceil(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Ceil(:,:,5));
%ylabel('gamma'); 
xlabel('\beta (rate of infection per unit of time)','Fontsize',24); %zlabel('RRMSE');
title('Infected \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Ceil(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');



%Standard Figure
figure()
tiledlayout(5,3);
nexttile
s=surf(spacer,spacer,MSE_s_SR(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Stand Susceptible \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_SR(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_SR(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_SR(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_SR(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_SR(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_SR(:,:,3));
ylabel('\gamma (rate of recovery per unit of time)','Fontsize',24); 
%xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_SR(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_SR(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_SR(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_SR(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_SR(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_SR(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_SR(:,:,5));
%ylabel('gamma'); 
xlabel('\beta (rate of infection per unit of time)','Fontsize',24); %zlabel('RRMSE');
title('Infected \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_SR(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');


%Ceiling + ResidualFigure
figure()
tiledlayout(5,3);
nexttile
s=surf(spacer,spacer,MSE_s_Ceil_R(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Ceiling+R Susceptible \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Ceil_R(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Ceil_R(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Ceil_R(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Ceil_R(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Ceil_R(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Ceil_R(:,:,3));
ylabel('\gamma (rate of recovery per unit of time)','Fontsize',24); 
%xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Ceil_R(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Ceil_R(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Ceil_R(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Ceil_R(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Ceil_R(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s_Ceil_R(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i_Ceil_R(:,:,5));
%ylabel('gamma'); 
xlabel('\beta (rate of infection per unit of time)','Fontsize',24); %zlabel('RRMSE');
title('Infected \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r_Ceil_R(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');

save("roundingerrorparams_5var_20tau.mat")

%end

%end of Run Section 5












%% Run_Section 6, loads params from different time steps, showing RRMSE with dif. time steps

elseif Run_Section==6

load("1stepvariables.mat")
MSE_S_1_M=mean(MSE_s,"all");
MSE_I_1_M=mean(MSE_i,"all");
MSE_R_1_M=mean(MSE_r,"all");

load("20params.mat")
MSE_S_20_M=mean(MSE_s,"all");
MSE_I_20_M=mean(MSE_i,"all");
MSE_R_20_M=mean(MSE_r,"all");

load("40params.mat")
MSE_S_40_M=mean(MSE_s,"all");
MSE_I_40_M=mean(MSE_i,"all");
MSE_R_40_M=mean(MSE_r,"all");

load("60params.mat")
MSE_S_60_M=mean(MSE_s,"all");
MSE_I_60_M=mean(MSE_i,"all");
MSE_R_60_M=mean(MSE_r,"all");

load("80values.mat")
MSE_S_80_M=mean(MSE_s,"all");
MSE_I_80_M=mean(MSE_i,"all");
MSE_R_80_M=mean(MSE_r,"all");

timidivi_col_S=[MSE_S_1_M,MSE_S_20_M,MSE_S_40_M,MSE_S_60_M,MSE_S_80_M];
timidivi_col_I=[MSE_I_1_M,MSE_I_20_M,MSE_I_40_M,MSE_I_60_M,MSE_I_80_M];
timidivi_col_R=[MSE_R_1_M,MSE_R_20_M,MSE_R_40_M,MSE_R_60_M,MSE_R_80_M];
figure()
plot([1,20,40,60,80], timidivi_col_S, 'b--o','LineWidth',2);
hold on
plot([1,20,40,60,80], timidivi_col_I,'c--o','LineWidth',2);
plot([1,20,40,60,80], timidivi_col_R, 'r--o','LineWidth',2);
legend('Susceptible_{RRMSE}','Infected_{RRMS}','Recovered_{RRMSE}');
title("SIR Petri Net vs ODE RRMSE based on times steps", 'FontSize', 24)
%subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' delta='+delta*timidivi)
xlabel('PN Time steps (per 1 ODE time interval)'); ylabel('RRMSE');
hold off

%End of Run Section 6















%% Run_Section 7, plots comutation time of single run of dif. time steps
elseif Run_Section==7

%values obtained from finding mean run time of runs both on laptop and SOL
%supercomputer

% figure()
% plot([1,20,40,60,80], [10.64,160.51,284.77,605.17,994.09], 'r-o','LineWidth',2);
% hold on
% legend('Computation Time');
% title("Computation Rime for PN Depending on Times Steps (\tau)", 'FontSize', 24)
% subtitle("This is mean run time for 1 PN model")
% xlabel('Ratio of Petri Net Time Steps to ODE Unit Time (\tau)'); ylabel('Mean Computation Time For 1 PN Simulation (seconds)');
% hold off

figure()
plot([1,2,3,4,5], [605.17,2317.87,5299.18,9101.68,13261.52], 'r-o','LineWidth',2);
hold on
legend('Computation Time');
title("Effect of Population Scalar on Computation Time", 'FontSize', 24)
xlabel('Population Scalar'); ylabel('Mean Computation Time For 1 PN Simulation (seconds)');
hold off

%end of Run section 7

















%% Run_Section 8, takes in CSV files from Spike run and compares to ode 
elseif Run_Section==8


    timidivi=1;       %tau in paper
MAX_ITERATIONS = 100*timidivi;
    timescalar=1;
    beta_array=[0, 0.1111, 0.2222, 0.3333, 0.4444, 0.5556, 0.6667, 0.7778, 0.8889, 1.0000];
    beta_array_sc=beta_array/timescalar;
    gamma_array=[0, 0.1111, 0.2222, 0.3333, 0.4444, 0.5556, 0.6667, 0.7778, 0.8889, 1.0000];
    gamma_array_sc=gamma_array/timescalar;
    delta_array=[0, 0.0010, 0.0100, 0.1000, 1.0000];
     delta_array_sc= delta_array/timescalar;
    S_0=1000;
    I_0=10;
    R_0=10;
    
    tau=100; %if tere is a relevant tau in spike code
    deltactr=1;
  for dai=1:1:5    %grid of params for delta (must be of size 5 for figs)
     betactr=1;
     for bai=1:1:10   
         gammactr=1;
        for gai=1:1:10
                %load Spike runs
                
                
                % if dai==2
                % loadeddata=readmatrix("/Users/treckell/Downloads/LY1_LY2_new_11_27_24/LY2_new_11_27_24/tauLeaping/Stochastic_SIRS_LY2tauLeaping_"...
                %     +num2str(beta_array_sc(bai))+"_"+num2str(gamma_array_sc(gai))+"_"+num2str(delta_array_sc(dai),'%.6f')...
                %     +".csv");
                % 
                % else
                %     loadeddata=readmatrix("/Users/treckell/Downloads/LY1_LY2_new_11_27_24/LY2_new_11_27_24/tauLeaping/Stochastic_SIRS_LY2tauLeaping_"...
                %     +num2str(beta_array_sc(bai))+"_"+num2str(gamma_array_sc(gai))+"_"+num2str(delta_array_sc(dai))...
                %     +".csv");
                % end
                
                % loadeddata=readmatrix("/Users/treckell/Downloads/New_Stoch_meths_LAYOUT_1_1000trials_/Direct_stoch_meth_results_1000trials_layout_1/Stochastic_SIRS_1000_trials_direct_"...
                %     +num2str(beta_array(bai),'%.4f')+"_"+num2str(gamma_array(gai),'%.4f')+"_"+num2str(delta_array(dai),'%.4f')...
                %     +".csv");
                
                % loadeddata=readmatrix("/Users/treckell/Downloads/Run tau 1/SIRS-SPN_SOStochastic_Branching_stochastic_"...
                %     +num2str(beta_array(bai))+"_"+num2str(gamma_array(gai))+"_"+num2str(delta_array(dai))...
                %     +"_"+num2str(10000)+".csv");
               
                % MS_loaded_pre=loadeddata(:,4);
                % MS_loaded= MS_loaded_pre(1:tau:end);
                % 
                % MI_loaded_pre=loadeddata(:,2);
                % MI_loaded= MI_loaded_pre(1:tau:end);
                % 
                % MR_loaded_pre=loadeddata(:,3);
                % MR_loaded=MR_loaded_pre(1:tau:end);

% loadeddata=readmatrix("/Users/treckell/Downloads/Layout_2_new_runs_direct/Layout_2_new_runs_direct/Stochastic_SIRS_LY2direct_"...
%                     +num2str(beta_array_sc(bai))+"_"+num2str(gamma_array_sc(gai))+"_"+num2str(delta_array_sc(dai))...
%                     +".csv");
% 
%                 % MS_loaded_pre=loadeddata(:,6)/100;
%                 % MS_loaded= MS_loaded_pre(1:100:10000);
%                 MS_loaded=loadeddata(1:100:10000,6)/100;
%                 % MI_loaded_pre=loadeddata(:,2)/100;
%                 % MI_loaded= MI_loaded_pre(1:100:10000);
%                 MI_loaded=loadeddata(1:100:10000,2)/100;
%                 % MR_loaded_pre=loadeddata(:,4)/100;
%                 % MR_loaded=MR_loaded_pre(1:100:10000);
%                 MR_loaded=loadeddata(1:100:10000,4)/100;

               loadeddata=readtable("/Users/treckell/Downloads/Layout_2_all_methods_pop_scalar10__12_05_24/Layout_2_all_methods_pop_scalar10__12_05_24_/deltaLeaping/Stochastic_SIRS_LY2deltaLeaping_"...
                    +num2str(beta_array_sc(bai))+"_"+num2str(gamma_array_sc(gai))+"_"+num2str(delta_array_sc(dai))...
                    +".csv");

                MS_loaded=table2array(loadeddata(1:100:10000,6))/10;
                MI_loaded=table2array(loadeddata(1:100:10000,2))/10;
                MR_loaded=table2array(loadeddata(1:100:10000,4))/10;

                y0=[S_0;I_0;R_0];
                delta_ODE=delta_array(dai); 
                beta_ODE=beta_array(bai); 
                gamma_ODE=gamma_array(gai);
                params=[delta_ODE;beta_ODE;gamma_ODE];
                [t,y]=ode89(@SIR,[1:1:100],y0,[],params);


                % RRMSE
                %RRMSE for Standard+Residual found 
                err_S_loaded = (rmse(MS_loaded(1:1:MAX_ITERATIONS) , y(:,1))/sqrt(sumsqr(y(:,1))))*100;
                 % err_SL2_pop100 = (rmse(Spike_SL2_pop100(1:1:MAX_ITERATIONS) , y(:,1))/sqrt(sumsqr(y(:,1))))*100;
                err_I_loaded = (rmse(MI_loaded(1:1:MAX_ITERATIONS) , y(:,2))/sqrt(sumsqr( y(:,2))))*100;
                err_R_loaded = (rmse(MR_loaded(1:1:MAX_ITERATIONS) , y(:,3))/sqrt(sumsqr(y(:,3))))*100;

if bai==7 && gai==7 && dai==4
    disp([num2str(err_S_loaded)]);
    disp([num2str(err_I_loaded)]);
    disp([num2str(err_R_loaded)]);
    MS_loaded_7_7_4=MS_loaded;
    ytes=y;
end
                disp(['still going: ', num2str(deltactr)]);

                % err_SL2_pop100 = (rmse(Spike_SL2_pop100(1:1:MAX_ITERATIONS) , y(:,1))/sqrt(sumsqr(y(:,1))))*100;
                % err_IL2_pop100 = (rmse(Spike_IL2_pop100(1:1:MAX_ITERATIONS) , y(:,2))/sqrt(sumsqr( y(:,2))))*100;
                % err_RL2_pop100 = (rmse(Spike_RL2_pop100(1:1:MAX_ITERATIONS) , y(:,3))/sqrt(sumsqr(y(:,3))))*100;




                %RRMSE Stored
                MSE_s_loaded(gammactr,betactr,deltactr)=err_S_loaded;
                MSE_i_loaded(gammactr,betactr,deltactr)=err_I_loaded;
                MSE_r_loaded(gammactr,betactr,deltactr)=err_R_loaded;
            % S=y(:,1);
            % I=y(:,2);
            % R=y(:,3);
            % SIRS=[S,I,R];
            %     writematrix(SIRS,"SIRS_ODE_"+num2str(beta_array(bai))+"_"+num2str(gamma_array(gai))+"_"+num2str(delta_array(dai))+".csv")
                % save("SIRS_ODE_"+num2str(beta_array(bai))+"_"+num2str(gamma_array(gai))+"_"+num2str(delta_array(dai))+".txt","S","I","R");
            gammactr=gammactr+1;
         end
        betactr=betactr+1;
     end
     deltactr=deltactr+1;
  end

casxismin=0;
casxismax=44;          %use if max is known to better customize color map
Cus_map=[(ones(2,3).*[0 0.5 0]);(ones(8,3).*[0 1 0]);(ones(40,3).*[1 1 0]);(ones(50,3).*[1 0.75 0]);...
    (ones(100,3).*[1 0.4 0]);(ones(240,3).*[1 0 0])];

%Loaded Figure
figure()
tiledlayout(5,3);
nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,3));
ylabel('\gamma (rate of recovery per unit of time)','Fontsize',24); 
%xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,5));
%ylabel('gamma'); 
xlabel('\beta (rate of infection per unit of time)','Fontsize',24); %zlabel('RRMSE');
title('Infected \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');


 %End of Section 8














%% Run_Section 9, takes in single delta run from Gpensim runs to combine into 1 figure 

elseif Run_Section==9
    beta_array=[0, 0.1111, 0.2222, 0.3333, 0.4444, 0.5556, 0.6667, 0.7778, 0.8889, 1.0000];
    gamma_array=[0, 0.1111, 0.2222, 0.3333, 0.4444, 0.5556, 0.6667, 0.7778, 0.8889, 1.0000];
    delta_array=[0, 0.0010, 0.0100, 0.1000, 1.0000];
    MSE_s_loaded=zeros(10,10,5);
    MSE_i_loaded=zeros(10,10,5);
    MSE_r_loaded=zeros(10,10,5);

    load('0_Popscalar_param.mat')
    MSE_s_loaded(:,:,1)=MSE_s(:,:,1);
    MSE_i_loaded(:,:,1)=MSE_i(:,:,1);
    MSE_r_loaded(:,:,1)=MSE_r(:,:,1);
    load('0.001_Popscalar_param.mat')
    MSE_s_loaded(:,:,2)=MSE_s(:,:,1);
    MSE_i_loaded(:,:,2)=MSE_i(:,:,1);
    MSE_r_loaded(:,:,2)=MSE_r(:,:,1);
    load('0.01_Popscalar_param.mat')
    MSE_s_loaded(:,:,3)=MSE_s(:,:,1);
    MSE_i_loaded(:,:,3)=MSE_i(:,:,1);
    MSE_r_loaded(:,:,3)=MSE_r(:,:,1);
    load('0.1_Popscalar_param.mat')
    MSE_s_loaded(:,:,4)=MSE_s(:,:,1);
    MSE_i_loaded(:,:,4)=MSE_i(:,:,1);
    MSE_r_loaded(:,:,4)=MSE_r(:,:,1);
    load('1_Popscalar_param.mat')
    MSE_s_loaded(:,:,5)=MSE_s(:,:,1);
    MSE_i_loaded(:,:,5)=MSE_i(:,:,1);
    MSE_r_loaded(:,:,5)=MSE_r(:,:,1);


casxismin=0;
casxismax=44;          %use if max is known to better customize color map
Cus_map=[(ones(2,3).*[0 0.5 0]);(ones(8,3).*[0 1 0]);(ones(40,3).*[1 1 0]);(ones(50,3).*[1 0.75 0]);...
    (ones(100,3).*[1 0.4 0]);(ones(240,3).*[1 0 0])];



%Loaded Figure
figure()
tiledlayout(5,3);
nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,3));
ylabel('\gamma (rate of recovery per unit of time)','Fontsize',24); 
%xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,5));
%ylabel('gamma'); 
xlabel('\beta (rate of infection per unit of time)','Fontsize',24); %zlabel('RRMSE');
title('Infected \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');

 %end of Section 9







elseif Run_Section==10
    beta_array=[0, 0.1111, 0.2222, 0.3333, 0.4444, 0.5556, 0.6667, 0.7778, 0.8889, 1.0000];
    gamma_array=[0, 0.1111, 0.2222, 0.3333, 0.4444, 0.5556, 0.6667, 0.7778, 0.8889, 1.0000];
    delta_array=[0, 0.0010, 0.0100, 0.1000, 1.0000];
    MSE_s_loaded=zeros(10,10,5);
    MSE_i_loaded=zeros(10,10,5);
    MSE_r_loaded=zeros(10,10,5);

    % load('0_Popscalar_param_0_to_0_2222.mat');
    % MSE_s_loaded(:,1:3,1)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,1:3,1)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,1:3,1)=MSE_r(:,1:3,1);
    % load('0_Popscalar_param_0_33_to_0_55.mat');
    % MSE_s_loaded(:,4:6,1)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,4:6,1)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,4:6,1)=MSE_r(:,1:3,1);
    % load('0_Popscalar_param_0_66_to_0_77.mat');
    % MSE_s_loaded(:,7:8,1)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,7:8,1)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,7:8,1)=MSE_r(:,1:2,1);
    % load('0_Popscalar_param_0_88_to_1.mat');
    % MSE_s_loaded(:,9:10,1)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,9:10,1)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,9:10,1)=MSE_r(:,1:2,1);
    % 
    % load('0.001_Popscalar_param_0_to_0_2222.mat');
    % MSE_s_loaded(:,1:3,2)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,1:3,2)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,1:3,2)=MSE_r(:,1:3,1);
    % load('0.001_Popscalar_param_0_33_to_0_55.mat');
    % MSE_s_loaded(:,4:6,2)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,4:6,2)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,4:6,2)=MSE_r(:,1:3,1);
    % load('0.001_Popscalar_param_0_66_to_0_77.mat');
    % MSE_s_loaded(:,7:8,2)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,7:8,2)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,7:8,2)=MSE_r(:,1:2,1);
    % load('0.001_Popscalar_param_0_88_to_1.mat');
    % MSE_s_loaded(:,9:10,2)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,9:10,2)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,9:10,2)=MSE_r(:,1:2,1);
    % 
    % 
    % load('0.01_Popscalar_param_0_to_0_2222.mat');
    % MSE_s_loaded(:,1:3,3)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,1:3,3)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,1:3,3)=MSE_r(:,1:3,1);
    % load('0.01_Popscalar_param_0_33_to_0_55.mat');
    % MSE_s_loaded(:,4:6,3)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,4:6,3)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,4:6,3)=MSE_r(:,1:3,1);
    % load('0.01_Popscalar_param_0_66_to_0_77.mat');
    % MSE_s_loaded(:,7:8,3)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,7:8,3)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,7:8,3)=MSE_r(:,1:2,1);
    % load('0.01_Popscalar_param_0_88_to_1.mat');
    % MSE_s_loaded(:,9:10,3)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,9:10,3)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,9:10,3)=MSE_r(:,1:2,1);
    % 
    % load('0.1_Popscalar_param_0_to_0_2222.mat');
    % MSE_s_loaded(:,1:3,4)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,1:3,4)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,1:3,4)=MSE_r(:,1:3,1);
    % load('0.1_Popscalar_param_0_33_to_0_55.mat');
    % MSE_s_loaded(:,4:6,4)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,4:6,4)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,4:6,4)=MSE_r(:,1:3,1);
    % load('0.1_Popscalar_param_0_66_to_0_77.mat');
    % MSE_s_loaded(:,7:8,4)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,7:8,4)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,7:8,4)=MSE_r(:,1:2,1);
    % load('0.1_Popscalar_param_0_88_to_1.mat');
    % MSE_s_loaded(:,9:10,4)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,9:10,4)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,9:10,4)=MSE_r(:,1:2,1);
    % 
    % load('1_Popscalar_param_0_to_0_2222.mat');
    % MSE_s_loaded(:,1:3,5)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,1:3,5)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,1:3,5)=MSE_r(:,1:3,1);
    % load('1_Popscalar_param_0_33_to_0_55.mat');
    % MSE_s_loaded(:,4:6,5)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,4:6,5)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,4:6,5)=MSE_r(:,1:3,1);
    % load('1_Popscalar_param_0_66_to_0_77.mat');
    % MSE_s_loaded(:,7:8,5)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,7:8,5)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,7:8,5)=MSE_r(:,1:2,1);
    % load('1_Popscalar_param_0_88_to_1.mat');
    % MSE_s_loaded(:,9:10,5)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,9:10,5)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,9:10,5)=MSE_r(:,1:2,1);
   
    % %Load for pop scalar 3
    % load('0_Popscalar_3_param_0_to_0_22.mat');
    % MSE_s_loaded(:,1:3,1)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,1:3,1)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,1:3,1)=MSE_r(:,1:3,1);
    % load('0_Popscalar_3_param_0_33_to_0_55.mat');
    % MSE_s_loaded(:,4:6,1)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,4:6,1)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,4:6,1)=MSE_r(:,1:3,1);
    % load('0_Popscalar_3_param_0_66_to_0_77.mat');
    % MSE_s_loaded(:,7:8,1)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,7:8,1)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,7:8,1)=MSE_r(:,1:2,1);
    % load('0_Popscalar_3_param_0_88_to_1.mat');
    % MSE_s_loaded(:,9:10,1)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,9:10,1)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,9:10,1)=MSE_r(:,1:2,1);
    % 
    % %load('0.001_Popscalar_3_param_0_to_0_22.mat');
    % load('0.0099999_Popscalar_3_param_0_to_0_22.mat')
    % MSE_s_loaded(:,1:3,2)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,1:3,2)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,1:3,2)=MSE_r(:,1:3,1);
    % load('0.001_Popscalar_3_param_0_33_to_0_55.mat');
    % MSE_s_loaded(:,4:6,2)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,4:6,2)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,4:6,2)=MSE_r(:,1:3,1);
    % load('0.001_Popscalar_3_param_0_66_to_0_77.mat');
    % MSE_s_loaded(:,7:8,2)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,7:8,2)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,7:8,2)=MSE_r(:,1:2,1);
    % load('0.001_Popscalar_3_param_0_88_to_1.mat');
    % MSE_s_loaded(:,9:10,2)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,9:10,2)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,9:10,2)=MSE_r(:,1:2,1);
    % 
    % 
    % %load('0.01_Popscalar_3_param_0_to_0_22.mat');
    % load('0.0009999_Popscalar_3_param_0_to_0_22.mat')
    % MSE_s_loaded(:,1:3,3)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,1:3,3)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,1:3,3)=MSE_r(:,1:3,1);
    % load('0.01_Popscalar_3_param_0_33_to_0_55.mat');
    % MSE_s_loaded(:,4:6,3)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,4:6,3)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,4:6,3)=MSE_r(:,1:3,1);
    % load('0.01_Popscalar_3_param_0_66_to_0_77.mat');
    % MSE_s_loaded(:,7:8,3)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,7:8,3)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,7:8,3)=MSE_r(:,1:2,1);
    % load('0.01_Popscalar_3_param_0_88_to_1.mat');
    % MSE_s_loaded(:,9:10,3)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,9:10,3)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,9:10,3)=MSE_r(:,1:2,1);
    % 
    % load('0.1_Popscalar_3_param_0_to_0_22.mat'); %
    % MSE_s_loaded(:,1:3,4)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,1:3,4)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,1:3,4)=MSE_r(:,1:3,1);
    % load('0.1_Popscalar_3_param_0_33_to_0_55.mat');%
    % MSE_s_loaded(:,4:6,4)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,4:6,4)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,4:6,4)=MSE_r(:,1:3,1);
    % load('0.1_Popscalar_3_param_0_66_to_0_77.mat');
    % MSE_s_loaded(:,7:8,4)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,7:8,4)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,7:8,4)=MSE_r(:,1:2,1);
    % load('0.1_Popscalar_3_param_0_88_to_1.mat');
    % MSE_s_loaded(:,9:10,4)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,9:10,4)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,9:10,4)=MSE_r(:,1:2,1);
    % 
    % load('1_Popscalar_3_param_0_to_0_22.mat');
    % MSE_s_loaded(:,1:3,5)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,1:3,5)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,1:3,5)=MSE_r(:,1:3,1);
    % load('1_Popscalar_3_param_0_33_to_0_55.mat');
    % MSE_s_loaded(:,4:6,5)=MSE_s(:,1:3,1);
    % MSE_i_loaded(:,4:6,5)=MSE_i(:,1:3,1);
    % MSE_r_loaded(:,4:6,5)=MSE_r(:,1:3,1);
    % load('1_Popscalar_3_param_0_66_to_0_77.mat');
    % MSE_s_loaded(:,7:8,5)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,7:8,5)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,7:8,5)=MSE_r(:,1:2,1);
    % load('1_Popscalar_3_param_0_88_to_1.mat');
    % MSE_s_loaded(:,9:10,5)=MSE_s(:,1:2,1);
    % MSE_i_loaded(:,9:10,5)=MSE_i(:,1:2,1);
    % MSE_r_loaded(:,9:10,5)=MSE_r(:,1:2,1);


% % Load for pop scalar 2
%     load('0_Popscalar_param.mat');
%     MSE_s_loaded(:,:,1)=MSE_s(:,:,1);
%     MSE_i_loaded(:,:,1)=MSE_i(:,:,1);
%     MSE_r_loaded(:,:,1)=MSE_r(:,:,1);
% 
%     load('0.001_Popscalar_param.mat');
%     MSE_s_loaded(:,:,2)=MSE_s(:,:,1);
%     MSE_i_loaded(:,:,2)=MSE_i(:,:,1);
%     MSE_r_loaded(:,:,2)=MSE_r(:,:,1);
% 
%     load('0.01_Popscalar_param.mat');
%     MSE_s_loaded(:,:,3)=MSE_s(:,:,1);
%     MSE_i_loaded(:,:,3)=MSE_i(:,:,1);
%     MSE_r_loaded(:,:,3)=MSE_r(:,:,1);
% 
%     load('0.1_Popscalar_param.mat');
%     MSE_s_loaded(:,:,4)=MSE_s(:,:,1);
%     MSE_i_loaded(:,:,4)=MSE_i(:,:,1);
%     MSE_r_loaded(:,:,4)=MSE_r(:,:,1);
% 
%     load('1_Popscalar_param.mat');
%     MSE_s_loaded(:,:,5)=MSE_s(:,:,1);
%     MSE_i_loaded(:,:,5)=MSE_i(:,:,1);
%     MSE_r_loaded(:,:,5)=MSE_r(:,:,1);

%Load for popscalar 5
    load('0_0_Popscalar_5_param.mat');
    MSE_s_loaded(:,1,1)=MSE_s(:,1,1);
    MSE_i_loaded(:,1,1)=MSE_i(:,1,1);
    MSE_r_loaded(:,1,1)=MSE_r(:,1,1);
    load('0_0.1111_Popscalar_5_param.mat');
    MSE_s_loaded(:,2,1)=MSE_s(:,1,1);
    MSE_i_loaded(:,2,1)=MSE_i(:,1,1);
    MSE_r_loaded(:,2,1)=MSE_r(:,1,1);
    load('0_0.2222_Popscalar_5_param.mat');
    MSE_s_loaded(:,3,1)=MSE_s(:,1,1);
    MSE_i_loaded(:,3,1)=MSE_i(:,1,1);
    MSE_r_loaded(:,3,1)=MSE_r(:,1,1);
    load('0_0.3333_Popscalar_5_param.mat');
    MSE_s_loaded(:,4,1)=MSE_s(:,1,1);
    MSE_i_loaded(:,4,1)=MSE_i(:,1,1);
    MSE_r_loaded(:,4,1)=MSE_r(:,1,1);
    load('0_0.4444_Popscalar_5_param.mat');
    MSE_s_loaded(:,5,1)=MSE_s(:,1,1);
    MSE_i_loaded(:,5,1)=MSE_i(:,1,1);
    MSE_r_loaded(:,5,1)=MSE_r(:,1,1);
    load('0_0.5556_Popscalar_5_param.mat');
    MSE_s_loaded(:,6,1)=MSE_s(:,1,1);
    MSE_i_loaded(:,6,1)=MSE_i(:,1,1);
    MSE_r_loaded(:,6,1)=MSE_r(:,1,1);
    load('0_0.6667_Popscalar_5_param.mat');
    MSE_s_loaded(:,7,1)=MSE_s(:,1,1);
    MSE_i_loaded(:,7,1)=MSE_i(:,1,1);
    MSE_r_loaded(:,7,1)=MSE_r(:,1,1);
    load('0_0.6667_Popscalar_5_param.mat');
    % load('0_0.7778_Popscalar_5_param.mat');
    MSE_s_loaded(:,8,1)=MSE_s(:,1,1);
    MSE_i_loaded(:,8,1)=MSE_i(:,1,1);
    MSE_r_loaded(:,8,1)=MSE_r(:,1,1);
    load('0_0.8889_Popscalar_5_param.mat');
    MSE_s_loaded(:,9,1)=MSE_s(:,1,1);
    MSE_i_loaded(:,9,1)=MSE_i(:,1,1);
    MSE_r_loaded(:,9,1)=MSE_r(:,1,1);
    load('0_1_Popscalar_5_param.mat');
    MSE_s_loaded(:,10,1)=MSE_s(:,1,1);
    MSE_i_loaded(:,10,1)=MSE_i(:,1,1);
    MSE_r_loaded(:,10,1)=MSE_r(:,1,1);

    load('0.001_0_Popscalar_5_param.mat');
    MSE_s_loaded(:,1,2)=MSE_s(:,1,1);
    MSE_i_loaded(:,1,2)=MSE_i(:,1,1);
    MSE_r_loaded(:,1,2)=MSE_r(:,1,1);
    load('0.001_0.1111_Popscalar_5_param.mat');
    MSE_s_loaded(:,2,2)=MSE_s(:,1,1);
    MSE_i_loaded(:,2,2)=MSE_i(:,1,1);
    MSE_r_loaded(:,2,2)=MSE_r(:,1,1);
    load('0.001_0.2222_Popscalar_5_param.mat');
    MSE_s_loaded(:,3,2)=MSE_s(:,1,1);
    MSE_i_loaded(:,3,2)=MSE_i(:,1,1);
    MSE_r_loaded(:,3,2)=MSE_r(:,1,1);
    load('0.001_0.3333_Popscalar_5_param.mat');
    MSE_s_loaded(:,4,2)=MSE_s(:,1,1);
    MSE_i_loaded(:,4,2)=MSE_i(:,1,1);
    MSE_r_loaded(:,4,2)=MSE_r(:,1,1);
    load('0.001_0.4444_Popscalar_5_param.mat');
    MSE_s_loaded(:,5,2)=MSE_s(:,1,1);
    MSE_i_loaded(:,5,2)=MSE_i(:,1,1);
    MSE_r_loaded(:,5,2)=MSE_r(:,1,1);
    load('0.001_0.5556_Popscalar_5_param.mat');
    MSE_s_loaded(:,6,2)=MSE_s(:,1,1);
    MSE_i_loaded(:,6,2)=MSE_i(:,1,1);
    MSE_r_loaded(:,6,2)=MSE_r(:,1,1);
    load('0.001_0.6667_Popscalar_5_param.mat');
    MSE_s_loaded(:,7,2)=MSE_s(:,1,1);
    MSE_i_loaded(:,7,2)=MSE_i(:,1,1);
    MSE_r_loaded(:,7,2)=MSE_r(:,1,1);
    load('0.001_0.7778_Popscalar_5_param.mat');
    MSE_s_loaded(:,8,2)=MSE_s(:,1,1);
    MSE_i_loaded(:,8,2)=MSE_i(:,1,1);
    MSE_r_loaded(:,8,2)=MSE_r(:,1,1);
    load('0.001_0.8889_Popscalar_5_param.mat');
    MSE_s_loaded(:,9,2)=MSE_s(:,1,1);
    MSE_i_loaded(:,9,2)=MSE_i(:,1,1);
    MSE_r_loaded(:,9,2)=MSE_r(:,1,1);
    load('0.001_1_Popscalar_5_param.mat');
    MSE_s_loaded(:,10,2)=MSE_s(:,1,1);
    MSE_i_loaded(:,10,2)=MSE_i(:,1,1);
    MSE_r_loaded(:,10,2)=MSE_r(:,1,1);

     load('0.01_0_Popscalar_5_param.mat');
    MSE_s_loaded(:,1,3)=MSE_s(:,1,1);
    MSE_i_loaded(:,1,3)=MSE_i(:,1,1);
    MSE_r_loaded(:,1,3)=MSE_r(:,1,1);
    load('0.01_0.1111_Popscalar_5_param.mat');
    MSE_s_loaded(:,2,3)=MSE_s(:,1,1);
    MSE_i_loaded(:,2,3)=MSE_i(:,1,1);
    MSE_r_loaded(:,2,3)=MSE_r(:,1,1);
    load('0.01_0.2222_Popscalar_5_param.mat');
    MSE_s_loaded(:,3,3)=MSE_s(:,1,1);
    MSE_i_loaded(:,3,3)=MSE_i(:,1,1);
    MSE_r_loaded(:,3,3)=MSE_r(:,1,1);
    load('0.01_0.3333_Popscalar_5_param.mat');
    MSE_s_loaded(:,4,3)=MSE_s(:,1,1);
    MSE_i_loaded(:,4,3)=MSE_i(:,1,1);
    MSE_r_loaded(:,4,3)=MSE_r(:,1,1);
    load('0.01_0.4444_Popscalar_5_param.mat');
    MSE_s_loaded(:,5,3)=MSE_s(:,1,1);
    MSE_i_loaded(:,5,3)=MSE_i(:,1,1);
    MSE_r_loaded(:,5,3)=MSE_r(:,1,1);
    load('0.01_0.5556_Popscalar_5_param.mat');
    MSE_s_loaded(:,6,3)=MSE_s(:,1,1);
    MSE_i_loaded(:,6,3)=MSE_i(:,1,1);
    MSE_r_loaded(:,6,3)=MSE_r(:,1,1);
    load('0.01_0.6667_Popscalar_5_param.mat');
    MSE_s_loaded(:,7,3)=MSE_s(:,1,1);
    MSE_i_loaded(:,7,3)=MSE_i(:,1,1);
    MSE_r_loaded(:,7,3)=MSE_r(:,1,1);
    load('0.01_0.7778_Popscalar_5_param.mat');
    MSE_s_loaded(:,8,3)=MSE_s(:,1,1);
    MSE_i_loaded(:,8,3)=MSE_i(:,1,1);
    MSE_r_loaded(:,8,3)=MSE_r(:,1,1);
    load('0.01_0.8889_Popscalar_5_param.mat');
    MSE_s_loaded(:,9,3)=MSE_s(:,1,1);
    MSE_i_loaded(:,9,3)=MSE_i(:,1,1);
    MSE_r_loaded(:,9,3)=MSE_r(:,1,1);
    load('0.01_1_Popscalar_5_param.mat');
    MSE_s_loaded(:,10,3)=MSE_s(:,1,1);
    MSE_i_loaded(:,10,3)=MSE_i(:,1,1);
    MSE_r_loaded(:,10,3)=MSE_r(:,1,1);

     load('0.1_0_Popscalar_5_param.mat');
    MSE_s_loaded(:,1,4)=MSE_s(:,1,1);
    MSE_i_loaded(:,1,4)=MSE_i(:,1,1);
    MSE_r_loaded(:,1,4)=MSE_r(:,1,1);
    load('0.1_0.1111_Popscalar_5_param.mat');
    MSE_s_loaded(:,2,4)=MSE_s(:,1,1);
    MSE_i_loaded(:,2,4)=MSE_i(:,1,1);
    MSE_r_loaded(:,2,4)=MSE_r(:,1,1);
    load('0.1_0.2222_Popscalar_5_param.mat');
    MSE_s_loaded(:,3,4)=MSE_s(:,1,1);
    MSE_i_loaded(:,3,4)=MSE_i(:,1,1);
    MSE_r_loaded(:,3,4)=MSE_r(:,1,1);
    load('0.1_0.3333_Popscalar_5_param.mat');
    MSE_s_loaded(:,4,4)=MSE_s(:,1,1);
    MSE_i_loaded(:,4,4)=MSE_i(:,1,1);
    MSE_r_loaded(:,4,4)=MSE_r(:,1,1);
    load('0.1_0.4444_Popscalar_5_param.mat');
    MSE_s_loaded(:,5,4)=MSE_s(:,1,1);
    MSE_i_loaded(:,5,4)=MSE_i(:,1,1);
    MSE_r_loaded(:,5,4)=MSE_r(:,1,1);
    load('0.1_0.5556_Popscalar_5_param.mat');
    MSE_s_loaded(:,6,4)=MSE_s(:,1,1);
    MSE_i_loaded(:,6,4)=MSE_i(:,1,1);
    MSE_r_loaded(:,6,4)=MSE_r(:,1,1);
    load('0.1_0.6667_Popscalar_5_param.mat');
    MSE_s_loaded(:,7,4)=MSE_s(:,1,1);
    MSE_i_loaded(:,7,4)=MSE_i(:,1,1);
    MSE_r_loaded(:,7,4)=MSE_r(:,1,1);
    load('0.1_0.7778_Popscalar_5_param.mat');
    MSE_s_loaded(:,8,4)=MSE_s(:,1,1);
    MSE_i_loaded(:,8,4)=MSE_i(:,1,1);
    MSE_r_loaded(:,8,4)=MSE_r(:,1,1);
    load('0.1_0.8889_Popscalar_5_param.mat');
    MSE_s_loaded(:,9,4)=MSE_s(:,1,1);
    MSE_i_loaded(:,9,4)=MSE_i(:,1,1);
    MSE_r_loaded(:,9,4)=MSE_r(:,1,1);
    load('0.1_1_Popscalar_5_param.mat');
    MSE_s_loaded(:,10,4)=MSE_s(:,1,1);
    MSE_i_loaded(:,10,4)=MSE_i(:,1,1);
    MSE_r_loaded(:,10,4)=MSE_r(:,1,1);

     load('1_0_Popscalar_5_param.mat');
    MSE_s_loaded(:,1,5)=MSE_s(:,1,1);
    MSE_i_loaded(:,1,5)=MSE_i(:,1,1);
    MSE_r_loaded(:,1,5)=MSE_r(:,1,1);
    load('1_0.1111_Popscalar_5_param.mat');
    MSE_s_loaded(:,2,5)=MSE_s(:,1,1);
    MSE_i_loaded(:,2,5)=MSE_i(:,1,1);
    MSE_r_loaded(:,2,5)=MSE_r(:,1,1);
    load('1_0.2222_Popscalar_5_param.mat');
    MSE_s_loaded(:,3,5)=MSE_s(:,1,1);
    MSE_i_loaded(:,3,5)=MSE_i(:,1,1);
    MSE_r_loaded(:,3,5)=MSE_r(:,1,1);
    load('1_0.3333_Popscalar_5_param.mat');
    MSE_s_loaded(:,4,5)=MSE_s(:,1,1);
    MSE_i_loaded(:,4,5)=MSE_i(:,1,1);
    MSE_r_loaded(:,4,5)=MSE_r(:,1,1);
    load('1_0.4444_Popscalar_5_param.mat');
    MSE_s_loaded(:,5,5)=MSE_s(:,1,1);
    MSE_i_loaded(:,5,5)=MSE_i(:,1,1);
    MSE_r_loaded(:,5,5)=MSE_r(:,1,1);
    load('1_0.5556_Popscalar_5_param.mat');
    MSE_s_loaded(:,6,5)=MSE_s(:,1,1);
    MSE_i_loaded(:,6,5)=MSE_i(:,1,1);
    MSE_r_loaded(:,6,5)=MSE_r(:,1,1);
    load('1_0.6667_Popscalar_5_param.mat');
    MSE_s_loaded(:,7,5)=MSE_s(:,1,1);
    MSE_i_loaded(:,7,5)=MSE_i(:,1,1);
    MSE_r_loaded(:,7,5)=MSE_r(:,1,1);
    load('1_0.7778_Popscalar_5_param.mat');
    MSE_s_loaded(:,8,5)=MSE_s(:,1,1);
    MSE_i_loaded(:,8,5)=MSE_i(:,1,1);
    MSE_r_loaded(:,8,5)=MSE_r(:,1,1);
    load('1_0.8889_Popscalar_5_param.mat');
    MSE_s_loaded(:,9,5)=MSE_s(:,1,1);
    MSE_i_loaded(:,9,5)=MSE_i(:,1,1);
    MSE_r_loaded(:,9,5)=MSE_r(:,1,1);
    load('1_1_Popscalar_5_param.mat');
    MSE_s_loaded(:,10,5)=MSE_s(:,1,1);
    MSE_i_loaded(:,10,5)=MSE_i(:,1,1);
    MSE_r_loaded(:,10,5)=MSE_r(:,1,1);



casxismin=0;
casxismax=10;          %use if max is known to better customize color map
Cus_map=[(ones(1,3).*[0 0.5 0]);(ones(9,3).*[0 1 0]);(ones(10,3).*[1 1 0]);(ones(10,3).*[1 0.75 0]);...
    (ones(10,3).*[1 0.4 0]);(ones(60,3).*[1 0 0])];
%Loaded Figure
sclr_i=2;
sclr_r=100;
figure()
tiledlayout(5,1);
nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,1)/sclr_i+(rand(size(MSE_i_loaded(:,:,1)))/(sclr_r)));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('\sigma=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,2)/sclr_i+(rand(size(MSE_i_loaded(:,:,2)))/(sclr_r-16)));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('\sigma=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,3)/sclr_i+(rand(size(MSE_i_loaded(:,:,3)))/(sclr_r-24)));
ylabel('\gamma (rate of recovery per unit of time)','Fontsize',24);
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('\sigma=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,4)/sclr_i+(rand(size(MSE_i_loaded(:,:,4)))/(sclr_r-32)));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('\sigma=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,5)/sclr_i+(rand(size(MSE_i_loaded(:,:,5)))/(sclr_r-40)));
%ylabel('gamma'); 
xlabel('\beta (rate of infection per unit of time)','Fontsize',24); %zlabel('RRMSE');
title('\sigma=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');



casxismin=0;
casxismax=44;          %use if max is known to better customize color map
Cus_map=[(ones(2,3).*[0 0.5 0]);(ones(8,3).*[0 1 0]);(ones(40,3).*[1 1 0]);(ones(50,3).*[1 0.75 0]);...
    (ones(100,3).*[1 0.4 0]);(ones(240,3).*[1 0 0])];

%Loaded Figure
figure()
tiledlayout(5,3);
nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,1));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,2));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.001','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,3));
ylabel('\gamma (rate of recovery per unit of time)','Fontsize',24); 
%xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,3));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered \delta=0.01','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,4));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=0.1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_s_loaded(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_i_loaded(:,:,5));
%ylabel('gamma'); 
xlabel('\beta (rate of infection per unit of time)','Fontsize',24); %zlabel('RRMSE');
title('Infected \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(beta_array,gamma_array,MSE_r_loaded(:,:,5));
%ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered \delta=1','Fontsize',12)
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
view(2)
h=colorbar('Position',[0.931120331950207,0.11,0.014914910340541,0.816052332195677]);
title(h, 'RRMSE (%)');

 %End of section 10











%% Run_Section 11, run 1 SIRS PN model and compare to ODE
elseif Run_Section==11

tic

bai=7;
gai=7;
dai=4;
timescalar=1;
    beta_array=[0, 0.1111, 0.2222, 0.3333, 0.4444, 0.5556, 0.6667, 0.7778, 0.8889, 1.0000];
    beta_array_sc=beta_array/timescalar;
    gamma_array=[0, 0.1111, 0.2222, 0.3333, 0.4444, 0.5556, 0.6667, 0.7778, 0.8889, 1.0000];
    gamma_array_sc=gamma_array/timescalar;
    delta_array=[0, 0.0010, 0.0100, 0.1000, 1.0000];
     delta_array_sc= delta_array/timescalar;
% TL=readtable("/Users/treckell/Downloads/New_Stoch_meths_LAYOUT_1_1000trials_/Direct_stoch_meth_results_1000trials_layout_1/Stochastic_SIRS_1000_trials_direct_"...
%                     +num2str(beta_array(bai))+"_"+num2str(gamma_array(gai))+"_"+num2str(delta_array(dai))...
%                     +".csv");
TLL1_100=readtable("/Users/treckell/Downloads/Stochastic_SIRS_direct_with_timescaling/Stochastic_SIRS_int_0_100_4000_LY1direct_0.00250000000_0.00250000000_0.00002500000"...
                    +".csv");
TLL1_4000=readtable("/Users/treckell/Downloads/Stochastic_SIRS_direct_with_timescaling/Stochastic_SIRS_int_0_4000_4000_LY1direct_0.00250000000_0.00250000000_0.00002500000"...
                    +".csv");
TLL2_100=readtable("/Users/treckell/Downloads/Stochastic_SIRS_direct_with_timescaling/Stochastic_SIRS_int_0_100_4000_LY2direct_0.00250000000_0.00250000000_0.00002500000"...
                    +".csv");
TLL2_4000=readtable("/Users/treckell/Downloads/Stochastic_SIRS_direct_with_timescaling/Stochastic_SIRS_int_0_4000_4000_LY2direct_0.00250000000_0.00250000000_0.00002500000"...
                    +".csv");
TLL2_pop100=readtable("/Users/treckell/Downloads/Layout_2_new_runs_direct/Layout_2_new_runs_direct/Stochastic_SIRS_LY2direct_"...
                    +num2str(beta_array_sc(bai))+"_"+num2str(gamma_array_sc(gai))+"_"+num2str(delta_array_sc(dai))...
                    +".csv");

Spike_SL2_pop100=table2array(TLL2_pop100(1:100:10000,6))/100;
Spike_IL2_pop100=table2array(TLL2_pop100(1:100:10000,2))/100;
Spike_RL2_pop100=table2array(TLL2_pop100(1:100:10000,4))/100;


Spike_SL1_100=table2array(TLL1_100(1:100,6))/4;
Spike_IL1_100=table2array(TLL1_100(1:100,2))/4;
Spike_RL1_100=table2array(TLL1_100(1:100,4))/4;


Spike_SL1_4000=table2array(TLL1_4000(1:40:4000,6))/4;
Spike_IL1_4000=table2array(TLL1_4000(1:40:4000,2))/4;
Spike_RL1_4000=table2array(TLL1_4000(1:40:4000,4))/4;

Spike_SL2_100=table2array(TLL2_100(1:100,6))/4;
Spike_IL2_100=table2array(TLL2_100(1:100,2))/4;
Spike_RL2_100=table2array(TLL2_100(1:100,4))/4;

Spike_SL2_4000=table2array(TLL2_4000(1:40:4000,6))/4;
Spike_IL2_4000=table2array(TLL2_4000(1:40:4000,2))/4;
Spike_RL2_4000=table2array(TLL2_4000(1:40:4000,4))/4;

popscaler=1;
timidivi=1;       %tau in paper
MAX_ITERATIONS = 100*timidivi;
S_0=1000*popscaler;
I_0=10*popscaler;
R_0=10*popscaler;
%delta=delta_array(dai); beta=beta_array(bai); gamma=gamma_array(gai);       %parameter values 
%delta=0.001; beta=0.1; gamma=0.1;

tau=1; %if tere is a relevant tau in spike code

% loadeddata=readmatrix("/Users/treckell/Downloads/New_Stoch_meths_LAYOUT_1_1000trials_/tauLeaping_stoch_meth_results_1000trials_layout_1/Stochastic_SIRS_1000_trials_tauLeaping_"...
%                     +num2str(beta_array(bai))+"_"+num2str(gamma_array(gai))+"_"+num2str(delta_array(dai))...
%                     +".csv");
% 
%                 % loadeddata=readmatrix("/Users/treckell/Downloads/New_Stoch_meths_LAYOUT_1_1000trials_/Direct_stoch_meth_results_1000trials_layout_1/Stochastic_SIRS_1000_trials_direct_"...
%                 %     +num2str(beta_array(bai),'%.4f')+"_"+num2str(gamma_array(gai),'%.4f')+"_"+num2str(delta_array(dai),'%.4f')...
%                 %     +".csv");
% 
%                 % loadeddata=readmatrix("/Users/treckell/Downloads/Run tau 1/SIRS-SPN_SOStochastic_Branching_stochastic_"...
%                 %     +num2str(beta_array(bai))+"_"+num2str(gamma_array(gai))+"_"+num2str(delta_array(dai))...
%                 %     +"_"+num2str(10000)+".csv");
% 
%                 % MS_loaded_pre=loadeddata(:,4);
%                 % MS_loaded= MS_loaded_pre(1:tau:end);
%                 % 
%                 % MI_loaded_pre=loadeddata(:,2);
%                 % MI_loaded= MI_loaded_pre(1:tau:end);
%                 % 
%                 % MR_loaded_pre=loadeddata(:,3);
%                 % MR_loaded=MR_loaded_pre(1:tau:end);
% 
%                 MS_loaded_pre=loadeddata(:,6);
%                 MS_loaded= MS_loaded_pre(1:tau:end);
%                 Spike_S=table2array(TL(1:100,6));
% Spike_I=table2array(TL(1:100,2));
% Spike_R=table2array(TL(1:100,4));
% 
%                 MI_loaded_pre=loadeddata(:,2);
%                 MI_loaded= MI_loaded_pre(1:tau:end);
%                 Spike_I=table2array(TL(1:100,2));
% 
%                 MR_loaded_pre=loadeddata(:,4);
%                 MR_loaded=MR_loaded_pre(1:tau:end);

                y0=[S_0;I_0;R_0];
                delta_ODE=delta_array(dai); 
                beta_ODE=beta_array(bai); 
                gamma_ODE=gamma_array(gai);
                params=[delta_ODE;beta_ODE;gamma_ODE];
                [t,y]=ode89(@SIR,[1:1:100],y0,[],params);



% 
% y0=[S_0/popscaler;I_0/popscaler;R_0/popscaler];           %ICs for ODE
% delta_ODE=delta; 
% beta_ODE=beta; 
% gamma_ODE=gamma;
% params=[delta_ODE;beta_ODE;gamma_ODE];
% [t,y]=ode15s(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);    %ODE model run
% %[t,y]=ode89(@SIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);

% RRMSE Calculation
 % err_S = (rmse(MSusceptible(1:timidivi:MAX_ITERATIONS)/popscaler , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 % err_I = (rmse(MInfected(1:timidivi:MAX_ITERATIONS)/popscaler , y(:,2)')/sqrt(sumsqr( y(:,2)')))*100;
 % err_R = (rmse(MRecovered(1:timidivi:MAX_ITERATIONS)/popscaler , y(:,3)')/sqrt(sumsqr(y(:,3)')))*100;

 err_SL1_100 = (rmse(Spike_SL1_100(1:1:MAX_ITERATIONS) , y(:,1))/sqrt(sumsqr(y(:,1))))*100;
 err_IL1_100 = (rmse(Spike_IL1_100(1:1:MAX_ITERATIONS) , y(:,2))/sqrt(sumsqr( y(:,2))))*100;
 err_RL1_100 = (rmse(Spike_RL1_100(1:1:MAX_ITERATIONS) , y(:,3))/sqrt(sumsqr(y(:,3))))*100;

 err_SL1_4000 = (rmse(Spike_SL1_4000(1:1:MAX_ITERATIONS), y(:,1))/sqrt(sumsqr(y(:,1))))*100;
 err_IL1_4000 = (rmse(Spike_IL1_4000(1:1:MAX_ITERATIONS), y(:,2))/sqrt(sumsqr( y(:,2))))*100;
 err_RL1_4000 = (rmse(Spike_RL1_4000(1:1:MAX_ITERATIONS), y(:,3))/sqrt(sumsqr(y(:,3))))*100;

 err_SL2_100 = (rmse(Spike_SL2_100(1:1:MAX_ITERATIONS), y(:,1))/sqrt(sumsqr(y(:,1))))*100;
 err_IL2_100 = (rmse(Spike_IL2_100(1:1:MAX_ITERATIONS), y(:,2))/sqrt(sumsqr( y(:,2))))*100;
 err_RL2_100 = (rmse(Spike_RL2_100(1:1:MAX_ITERATIONS), y(:,3))/sqrt(sumsqr(y(:,3))))*100;

 err_SL2_4000 = (rmse(Spike_SL2_4000(1:1:MAX_ITERATIONS), y(:,1))/sqrt(sumsqr(y(:,1))))*100;
 err_IL2_4000 = (rmse(Spike_IL2_4000(1:1:MAX_ITERATIONS), y(:,2))/sqrt(sumsqr( y(:,2))))*100;
 err_RL2_4000 = (rmse(Spike_RL2_4000(1:1:MAX_ITERATIONS), y(:,3))/sqrt(sumsqr(y(:,3))))*100;
 
 err_SL2_pop100 = (rmse(Spike_SL2_pop100(1:1:MAX_ITERATIONS) , y(:,1))/sqrt(sumsqr(y(:,1))))*100;
 err_IL2_pop100 = (rmse(Spike_IL2_pop100(1:1:MAX_ITERATIONS) , y(:,2))/sqrt(sumsqr( y(:,2))))*100;
 err_RL2_pop100 = (rmse(Spike_RL2_pop100(1:1:MAX_ITERATIONS) , y(:,3))/sqrt(sumsqr(y(:,3))))*100;

% Plot Results
% figure()
% plot(t,y(:,1), 'k','LineWidth',10) 
% hold on
% plot(t,Spike_SL1_100, 'r','LineWidth',6)
% plot(t,Spike_SL1_4000,'g','LineWidth',2)
% plot(t,Spike_SL2_100, 'b','LineWidth',6)
% plot(t,Spike_SL2_4000,'m','LineWidth',2)
% legend('ODE','L1 100','L1 4000', ...
%     'L2 1000','L2 4000');
% title("Susceptible Direct SIR Petri Net vs ODE", 'FontSize', 24)
% subtitle("beta="+beta+', gamma='+gamma+', delta='+delta+', L1 100 RRMSE='+err_SL1_100+...
%     ', L1 4000 RRMSE='+err_SL1_4000+', L2 100 RRMSE='+err_SL2_100+', L2 4000 RRMSE='+err_SL2_4000,'FontSize', 6)
% ylabel('Tokens (population)'); xlabel('Time'); 
% hold off
% 
% 
% figure()
% plot(t,y(:,2), 'k','LineWidth',10) 
% hold on
% plot(t,Spike_IL1_100, 'r','LineWidth',6)
% plot(t,Spike_IL1_4000,'g','LineWidth',2)
% plot(t,Spike_IL2_100, 'b','LineWidth',6)
% plot(t,Spike_IL2_4000,'m','LineWidth',2)
% legend('ODE','L1 100','L1 4000', ...
%     'L2 1000','L2 4000');
% title("Infected Direct SIR Petri Net vs ODE", 'FontSize', 24)
% subtitle("beta="+beta+', gamma='+gamma+', delta='+delta+', L1 100 RRMSE='+err_IL1_100+...
%     ', L1 4000 RRMSE='+err_IL1_4000+', L2 100 RRMSE='+err_IL2_100+', L2 4000 RRMSE='+err_IL2_4000,'FontSize', 6)
% ylabel('Tokens (population)'); xlabel('Time'); 
% hold off
% 
% figure()
% plot(t,y(:,3), 'k','LineWidth',10) 
% hold on
% plot(t,Spike_RL1_100, 'r','LineWidth',6)
% plot(t,Spike_RL1_4000,'g','LineWidth',2)
% plot(t,Spike_RL2_100, 'b','LineWidth',6)
% plot(t,Spike_RL2_4000,'m','LineWidth',2)
% legend('ODE','L1 100','L1 4000', ...
%     'L2 1000','L2 4000');
% title("Recovered Direct SIR Petri Net vs ODE", 'FontSize', 24)
% subtitle("beta="+beta+', gamma='+gamma+', delta='+delta+', L1 100 RRMSE='+err_RL1_100+...
%     ', L1 4000 RRMSE='+err_RL1_4000+', L2 100 RRMSE='+err_RL2_100+', L2 4000 RRMSE='+err_RL2_4000,'FontSize', 6)
% ylabel('Tokens (population)'); xlabel('Time'); 
% hold off


figure()
plot(t,y(:,1), 'k','LineWidth',5) 
hold on
plot(t,Spike_SL2_pop100, 'r','LineWidth',3)
legend('ODE','L2 pop 100');
title("Susceptible Direct SIR PN pop 100 vs ODE", 'FontSize', 24)
subtitle("beta="+beta_ODE+', gamma='+gamma_ODE+', delta='+delta_ODE+', L2 pop 100 RRMSE='+err_SL2_pop100,'FontSize', 6)
ylabel('Tokens (population)'); xlabel('Time'); 
hold off


figure()
plot(t,y(:,2), 'k','LineWidth',5) 
hold on
plot(t,Spike_IL2_pop100, 'r','LineWidth',3)
legend('ODE','L2 pop 100');
title("Infected Direct SIR PN pop 100 vs ODE", 'FontSize', 24)
subtitle("beta="+beta_ODE+', gamma='+gamma_ODE+', delta='+delta_ODE+', L2 pop 100 RRMSE='+err_IL2_pop100,'FontSize', 6)
ylabel('Tokens (population)'); xlabel('Time'); 
hold off


figure()
plot(t,y(:,3), 'k','LineWidth',5) 
hold on
plot(t,Spike_RL2_pop100, 'r','LineWidth',3)
legend('ODE','L2 pop 100');
title("Recovered Direct SIR PN pop 100 vs ODE", 'FontSize', 24)
subtitle("beta="+beta_ODE+', gamma='+gamma_ODE+', delta='+delta_ODE+', L2 pop 100 RRMSE='+err_RL2_pop100,'FontSize', 6)
ylabel('Tokens (population)'); xlabel('Time'); 
hold off

end

%end of section 11


toc



%% ODE Function
function dSIR = SIR(t,a,params)
S = a(1); I = a(2); R = a(3);
delta=params(1); beta=params(2); gamma=params(3);
dS=delta*R-beta*S*I;
dI=beta*S*I-gamma*I;
dR=gamma*I-delta*R;
dSIR=[dS;dI;dR];

end


%% CALCR0 Calculate basic reproduction number (Eq 29)
function R0 = calcR0(w,r)
R0 = r/((exp(r) - 1)*sum(w.*exp(-r*(1:length(w)))));
end


%% CALCR Calculate instantaneous reproduction number (Eq 34)
function R = calcR(w,I)
K = length(w);
N = length(I);
R = NaN(N,1);
for n = 1:N
    s = 0;
    for k = 1:min(n,K)
        s = s + w(k)*I(n+1-k);
    end
    
    if s > 0
        R(n) = I(n)/s;
    end
end
end


%% CALCRC Calculate case reproduction number (Eq 37)
function [Rc,R] = calcRc(w,I)
R = calcR(w,I);
K = length(w);
N = length(I);
Rc = NaN(N,1);
for n = 1:N
    s = 0;
    for k = 1:min(N-n+1,K)
        s = s + w(k)*R(n+k-1);
    end
    Rc(n) = s;
end
end

% function y = optimization_func(x, params)
% 
%     y = params(1) * (params(2) * x);
% end
