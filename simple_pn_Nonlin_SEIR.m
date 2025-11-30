% Nonlin SEIRv1: A Simple Example for variable arc weights
% the main file to run simulation 
% run code first simple_pn_Nonlin_SEIR_pdf.m and COMMON_PRE.m 
 

clear all; clc; 
global global_info m0Susceptible m0Exposed m0Infected m0Recovered %beta
tic

timidivi_col_S=[];
timidivi_col_E=[];
timidivi_col_I=[];
timidivi_col_R=[];
%timidivi_timestep=1:6:79;
%for timidivi=1:6:79
ode_run_counter=1;
%timidivi=1;
%for timidivi=[2]
timidivi=80;
%% The following section is used for finding the RRMSE for various
%parameter values in comparison to the respective ODE, this is limited two
%comparison of 3 parameters in this configuration, for visual display.
  bglength=10;
  popscalar=1;
  %bglength=1;
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

MSE_s=[zeros(bglength,bglength,5)];
MSE_e=[zeros(bglength,bglength,5)];
MSE_i=[zeros(bglength,bglength,5)];
MSE_r=[zeros(bglength,bglength,5)];

mu=0.01;
alpha=0.01;
sigma=0.05;
gamma=0.05;
beta=0.05;
% MSE_S_c=zeros(1,5);
% MSE_S_r=zeros(1,5);
% MSE_S_v=zeros(1,5);
% MSE_S_C=zeros(bglength*bglength,5);
% %0.1
%sigmar=[0,logspace(-3,0,4)];
 sigmactr=1;
  for sigma=[0,logspace(-3,0,4)]
     betactr=1;
     %for beta=linspace(beta_ODE*0.9,beta_ODE*1.1,bglength) %beta_ODE = .0008;     gamma_ODE = 0.08;
     for beta=spacer   
         %beta=spacer(betaindex);
         gammactr=1;
         %for gamma=linspace(gamma_ODE*0.9,gamma_ODE*1.1,bglength)
         for gamma=spacer
             %gamma=spacer(gammaindex);
%sigma=0; beta=0; gamma=0.5;
% sigma1=1/timidivi; beta1=1/timidivi; gamma1=1/timidivi;
try
%alpha_T=alpha/timidivi; 
mu_T=mu/timidivi; sigma_T=sigma/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;
alpha_T=alpha;
%sigma=1; beta=1; gamma=1;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MExposed = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS); 


%% Parameter values  %%%%%%%%%%%%%%%%%%%%%
%Used for setting specific parameter values (vs loop above for grid of
%param values)

%beta = 0.001;     gamma = 0.05;    sigma = 0.001;
%beta = 0.1;     gamma = 0.2;    sigma = 0.1;
%beta = 0.2;     gamma = 0.5;    sigma = 0.4;
%beta = 3.1000e-04;     gamma =  0.0809;    
%sigma = 0.00;


%% Initial Population
S_0=1000*popscalar;
E_0=20*popscalar;
I_0=10*popscalar;
R_0=10*popscalar;

m0Susceptible = S_0; 
m0Exposed = E_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
t_1_pSus_resEr=0; pSus_t_2_resEr=0; 
pExp_t_4_resEr=0; pInf_t_6_resEr=0;
pRec_t_8_resEr=0; pSus_t_3_resEr=0;
t_3_pExp_resEr=0; pExp_t_5_resEr=0;
t_5_pInf_resEr=0; pInf_t_7_resEr=0;
t_7_pRec_resEr=0; 
zero=0;

Number_of_iterations = 1;

while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MExposed(Number_of_iterations) = m0Exposed; 
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 
    %  Arc weights
    %This below is for splitting time
    % if Number_of_iterations==1
    %     beta=beta1;
    %     sigma=sigma1;
    %     gamma=gamma1;
    % else
    %     % beta=1;
    %     % sigma=1;
    %     % gamma=1;
    % end
    % % if beta_T==0 || m0Infected==0
    % %     global_info.pSus_tInf = zero;
    % %     pSus_tInf_resEr=(beta_T * m0Susceptible *  m0Infected + pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
    % %     global_info.tInf_pInf = zero;
    % %     tInf_pInf_resEr=(beta_T * m0Susceptible *  m0Infected+tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    % % elseif m0Infected<=(1/(beta_T))
    % %     global_info.pSus_tInf = round(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr);
    % %     pSus_tInf_resEr=(beta_T * m0Susceptible *   m0Infected +pSus_tInf_resEr)-round(beta_T * m0Susceptible * m0Infected + pSus_tInf_resEr);
    % %     global_info.tInf_pInf = round(beta_T * m0Susceptible *   m0Infected + tInf_pInf_resEr);
    % %     tInf_pInf_resEr=(beta_T * m0Susceptible *   m0Infected +tInf_pInf_resEr)-round(beta_T * m0Susceptible * m0Infected + tInf_pInf_resEr);
    % % else 
    % % 
    % %     %beta2=(1/(m0Infected+??????));
    % %     beta2=1/(m0Infected+1);
    % %     global_info.pSus_tInf = round(beta2 * m0Susceptible * m0Infected + pSus_tInf_resEr);
    % %     pSus_tInf_resEr=(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr)-round(beta2 * m0Susceptible * m0Infected +pSus_tInf_resEr);
    % %     global_info.tInf_pInf = round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    % %     tInf_pInf_resEr=(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr)-round(beta2 * m0Susceptible * m0Infected +tInf_pInf_resEr);
    % % 
    % % end
    % % 
    % % global_info.pInf_tRec = round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    % % pInf_tRec_resEr=(gamma_T * (m0Infected) + pInf_tRec_resEr)-round(gamma_T * (m0Infected) + pInf_tRec_resEr);
    % % global_info.tRec_pRec = round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    % % tRec_pRec_resEr=(gamma_T * (m0Infected) + tRec_pRec_resEr)-round(gamma_T * (m0Infected) + tRec_pRec_resEr);
    % % global_info.pRec_tSus = round(sigma_T * (m0Recovered) + pRec_tSus_resEr);
    % % pRec_tSus_resEr=(sigma_T * (m0Recovered) + pRec_tSus_resEr)-round(sigma_T * (m0Recovered) + pRec_tSus_resEr);
    % % global_info.tSus_pSus = round(sigma_T * (m0Recovered) + tSus_pSus_resEr);
    % % tSus_pSus_resEr=(sigma_T * (m0Recovered)+ tSus_pSus_resEr)-round(sigma_T * (m0Recovered) + tSus_pSus_resEr);
    I_n = double(m0Infected) / popscalar; 
    S_n = double(m0Susceptible) / popscalar;
    nm = beta * S_n * I_n;
    dm = 1 + (alpha * I_n * I_n);
    rate_density = nm / dm;
    expected_transitions = rate_density* popscalar/timidivi;
    calc_value = expected_transitions + pSus_t_3_resEr;
    integer_fire = round(calc_value);

    if beta_T==0 || m0Infected==0
        global_info.pSus_t_3=zero;
        pSus_t_3_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr);
        global_info.t_3_pExp=zero;
        t_3_pExp_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr); 
    elseif (beta_T*m0Infected)<=(1+alpha_T*m0Infected*m0Infected)
        global_info.pSus_t_3 = integer_fire;
        pSus_t_3_resEr = calc_value - integer_fire;
        global_info.t_3_pExp = integer_fire;
        t_3_pExp_resEr = calc_value - integer_fire;
        
        % global_info.pSus_t_3=round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr);
        % pSus_t_3_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr);
        % global_info.t_3_pExp=round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr);
        % t_3_pExp_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr);
    else 
        beta2=1/(m0Infected+1);
        nm = beta2 * S_n * I_n;
        dm = 1 + (alpha * I_n * I_n);
        rate_density = nm / dm;
        expected_transitions = rate_density* popscalar/timidivi;
        calc_value = expected_transitions + pSus_t_3_resEr;
        integer_fire = round(calc_value);
         global_info.pSus_t_3 = integer_fire;
        pSus_t_3_resEr = calc_value - integer_fire;
        global_info.t_3_pExp = integer_fire;
        t_3_pExp_resEr = calc_value - integer_fire;
    end
    global_info.t_1_pSus=round(mu_T + t_1_pSus_resEr);
    t_1_pSus_resEr=(mu_T+t_1_pSus_resEr)-round(mu_T+t_1_pSus_resEr);
    global_info.pSus_t_2=round(mu_T* m0Susceptible + pSus_t_2_resEr);
    pSus_t_2_resEr=(mu_T* m0Susceptible + pSus_t_2_resEr)-round(mu_T* m0Susceptible + pSus_t_2_resEr); 
    global_info.pExp_t_4=round(mu_T * m0Exposed + pExp_t_4_resEr);
    pExp_t_4_resEr=(mu_T * m0Exposed + pExp_t_4_resEr)-round(mu_T * m0Exposed + pExp_t_4_resEr); 
    global_info.pInf_t_6=round(mu_T * m0Infected + pInf_t_6_resEr);
    pInf_t_6_resEr=(mu_T * m0Infected + pInf_t_6_resEr)-round(mu_T * m0Infected + pInf_t_6_resEr);
    global_info.pRec_t_8=round(mu_T * m0Recovered+pRec_t_8_resEr);
    pRec_t_8_resEr=(mu_T * m0Recovered+pRec_t_8_resEr)-round(mu_T * m0Recovered+pRec_t_8_resEr); 
    global_info.pSus_t_3=round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr);
    pSus_t_3_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr);
    global_info.t_3_pExp=round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr);
    t_3_pExp_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr); 
    global_info.pExp_t_5=round(sigma_T * m0Exposed + pExp_t_5_resEr);
    pExp_t_5_resEr=(sigma_T * m0Exposed + pExp_t_5_resEr)-round(sigma_T * m0Exposed + pExp_t_5_resEr);
    global_info.t_5_pInf=round(sigma_T * m0Exposed + t_5_pInf_resEr);
    t_5_pInf_resEr=(sigma_T * m0Exposed + t_5_pInf_resEr)-round(sigma_T * m0Exposed + t_5_pInf_resEr); 
    global_info.pInf_t_7=round(gamma_T * m0Infected + pInf_t_7_resEr);
    pInf_t_7_resEr=(gamma_T * m0Infected + pInf_t_7_resEr)-round(gamma_T * m0Infected + pInf_t_7_resEr);
    global_info.t_7_pRec=round(gamma_T * m0Infected + t_7_pRec_resEr);
    t_7_pRec_resEr=(gamma_T * m0Infected + t_7_pRec_resEr)-round(gamma_T * m0Infected + t_7_pRec_resEr); 

 %    global_info.t_1_pSus,... 
 % 'pSusceptible','t_2', global_info.pSus_t_2,... 
 % 'pExposed','t_4', global_info.pExp_t_4,... 
 % 'pInfected','t_6', global_info.pInf_t_6,... 
 % 'pRecovered','t_8',global_info.pRec_t_8,... 
 % 'pSusceptible','t_3',global_info.pSus_t_3,... 
 % 't_3','pExposed',global_info.t_3_pExp,... 
 % 'pExposed','t_5',global_info.pExp_t_5,...
 % 't_5','pInfected',global_info.t_5_pInf,... 
 % 'pInfected','t_7',global_info.pInf_t_7,... 
 % 't_7','pRecovered',global_info.t_7_pRec,... 

    pns = pnstruct('simple_pn_Nonlin_SEIR_pdf');
    dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
          'pExposed',m0Exposed,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs    

    m0Susceptible = ntokens('pSusceptible');
    m0Exposed = ntokens('pExposed');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1
end
catch

    disp(['Error at sigma=', int2str(sigma), ' , gamma=',int2str(gamma), ' , beta=',int2str(beta), ' , timesteps=',int2str(timidivi)]);
    
    if sigma==0
    sigma=sigma; 
    else
    sigma=sigma-10^-7; 
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
sigma_T=sigma/timidivi; beta_T=beta/timidivi; gamma_T=gamma/timidivi;

%% Max iterations for Petri Net and ODE, petri net does 1 step per iteration
%in base settings
MAX_ITERATIONS = 100*timidivi;

%%%%% Vectors to store the population after each iteration
MSusceptible = zeros(1, MAX_ITERATIONS); 
MExposed = zeros(1, MAX_ITERATIONS); 
MInfected = zeros(1, MAX_ITERATIONS); 
MRecovered = zeros(1, MAX_ITERATIONS);  

%% Initial Population
S_0=1000*popscalar;
E_0=20*popscalar;
I_0=10*popscalar;
R_0=10*popscalar;

m0Susceptible = S_0; 
m0Exposed = E_0; 
m0Infected = I_0;
m0Recovered = R_0; 

%% the loop   %%%%%%%%%%%%%%%%%%%%%%%%%
t_1_pSus_resEr=0; pSus_t_2_resEr=0; 
pExp_t_4_resEr=0; pInf_t_6_resEr=0;
pRec_t_8_resEr=0; pSus_t_3_resEr=0;
t_3_pExp_resEr=0; pExp_t_5_resEr=0;
t_5_pInf_resEr=0; pInf_t_7_resEr=0;
t_7_pRec_resEr=0; 
zero=0;


Number_of_iterations = 1;
while le(Number_of_iterations, MAX_ITERATIONS)
    %disp(['Iteration: ', int2str(Number_of_iterations)]);
    MSusceptible(Number_of_iterations) = m0Susceptible;
    MExposed(Number_of_iterations) = m0Exposed; 
    MInfected(Number_of_iterations) = m0Infected; 
    MRecovered(Number_of_iterations) = m0Recovered; 

    I_n = double(m0Infected) / popscalar; 
    S_n = double(m0Susceptible) / popscalar;
    nm = beta * S_n * I_n;
    dm = 1 + (alpha * I_n * I_n);
    rate_density = nm / dm;
    expected_transitions = rate_density* popscalar/timidivi;
    calc_value = expected_transitions + pSus_t_3_resEr;
    integer_fire = round(calc_value);

    if beta_T==0 || m0Infected==0
        global_info.pSus_t_3=zero;
        pSus_t_3_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr);
        global_info.t_3_pExp=zero;
        t_3_pExp_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr); 
    elseif (beta_T*m0Infected)<=(1+alpha_T*m0Infected*m0Infected)
        global_info.pSus_t_3 = integer_fire;
        pSus_t_3_resEr = calc_value - integer_fire;
        global_info.t_3_pExp = integer_fire;
        t_3_pExp_resEr = calc_value - integer_fire;
        
        % global_info.pSus_t_3=round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr);
        % pSus_t_3_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr);
        % global_info.t_3_pExp=round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr);
        % t_3_pExp_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr);
    else 
        beta2=1/(m0Infected+1);
        nm = beta2 * S_n * I_n;
        dm = 1 + (alpha * I_n * I_n);
        rate_density = nm / dm;
        expected_transitions = rate_density* popscalar/timidivi;
        calc_value = expected_transitions + pSus_t_3_resEr;
        integer_fire = round(calc_value);
         global_info.pSus_t_3 = integer_fire;
        pSus_t_3_resEr = calc_value - integer_fire;
        global_info.t_3_pExp = integer_fire;
        t_3_pExp_resEr = calc_value - integer_fire;
    end
    
    global_info.t_1_pSus=round(mu_T + t_1_pSus_resEr);
    t_1_pSus_resEr=(mu_T+t_1_pSus_resEr)-round(mu_T+t_1_pSus_resEr);
    global_info.pSus_t_2=round(mu_T* m0Susceptible + pSus_t_2_resEr);
    pSus_t_2_resEr=(mu_T* m0Susceptible + pSus_t_2_resEr)-round(mu_T* m0Susceptible + pSus_t_2_resEr); 
    global_info.pExp_t_4=round(mu_T * m0Exposed + pExp_t_4_resEr);
    pExp_t_4_resEr=(mu_T * m0Exposed + pExp_t_4_resEr)-round(mu_T * m0Exposed + pExp_t_4_resEr); 
    global_info.pInf_t_6=round(mu_T * m0Infected + pInf_t_6_resEr);
    pInf_t_6_resEr=(mu_T * m0Infected + pInf_t_6_resEr)-round(mu_T * m0Infected + pInf_t_6_resEr);
    global_info.pRec_t_8=round(mu_T * m0Recovered+pRec_t_8_resEr);
    pRec_t_8_resEr=(mu_T * m0Recovered+pRec_t_8_resEr)-round(mu_T * m0Recovered+pRec_t_8_resEr); 
    global_info.pSus_t_3=round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr);
    pSus_t_3_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+pSus_t_3_resEr);
    global_info.t_3_pExp=round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr);
    t_3_pExp_resEr=((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr)-round((beta_T* m0Susceptible * m0Infected)/(1+alpha_T* m0Infected* m0Infected)+t_3_pExp_resEr); 
    global_info.pExp_t_5=round(sigma_T * m0Exposed + pExp_t_5_resEr);
    pExp_t_5_resEr=(sigma_T * m0Exposed + pExp_t_5_resEr)-round(sigma_T * m0Exposed + pExp_t_5_resEr);
    global_info.t_5_pInf=round(sigma_T * m0Exposed + t_5_pInf_resEr);
    t_5_pInf_resEr=(sigma_T * m0Exposed + t_5_pInf_resEr)-round(sigma_T * m0Exposed + t_5_pInf_resEr); 
    global_info.pInf_t_7=round(gamma_T * m0Infected + pInf_t_7_resEr);
    pInf_t_7_resEr=(gamma_T * m0Infected + pInf_t_7_resEr)-round(gamma_T * m0Infected + pInf_t_7_resEr);
    global_info.t_7_pRec=round(gamma_T * m0Infected + t_7_pRec_resEr);
    t_7_pRec_resEr=(gamma_T * m0Infected + t_7_pRec_resEr)-round(gamma_T * m0Infected + t_7_pRec_resEr); 

    pns = pnstruct('simple_pn_SIRv1_pdf');
    %dyn.ft=1/timidivi;
    dyn.m0 = {'pSusceptible', m0Susceptible,...
            'pExposed',m0Exposed,...
          'pInfected', m0Infected, ...
          'pRecovered', m0Recovered};
    pni = initialdynamics(pns, dyn); 
    Sim_Results = gpensim(pni); % perform simulation runs

    m0Susceptible = ntokens('pSusceptible');
    m0Exposed= ntokens('pExposed');
    m0Infected = ntokens('pInfected');
    m0Recovered = ntokens('pRecovered');
    Number_of_iterations = Number_of_iterations + 1;
end
end

if ode_run_counter==1
y0=[S_0/popscalar;E_0/popscalar;I_0/popscalar;R_0/popscalar];
sigma_ODE=sigma; 
alpha_ODE=alpha;
beta_ODE=beta; 
gamma_ODE=gamma;
mu_ODE=mu;
params=[sigma_ODE;beta_ODE;gamma_ODE;mu_ODE;alpha_ODE];
[t,y]=ode15s(@NLSEIR,[1:1:(MAX_ITERATIONS/timidivi)],y0,[],params);
else
end

%% RRMSE
 err_S = (rmse(MSusceptible(1:timidivi:MAX_ITERATIONS) , y(:,1)')/sqrt(sumsqr(y(:,1)')))*100;
 err_E = (rmse(MExposed(1:timidivi:MAX_ITERATIONS) , y(:,2)')/sqrt(sumsqr(y(:,2)')))*100;
 err_I = (rmse(MInfected(1:timidivi:MAX_ITERATIONS) , y(:,3)')/sqrt(sumsqr( y(:,3)')))*100;
 err_R = (rmse(MRecovered(1:timidivi:MAX_ITERATIONS) , y(:,4)')/sqrt(sumsqr(y(:,4)')))*100;
 % % 
 % % timidivi_col_S=[timidivi_col_S err_S];
 % % timidivi_col_I=[timidivi_col_I err_I];
 % % timidivi_col_R=[timidivi_col_R err_R];
% %end
 %% Plot times steps vs rrmse %%%%%%%%%%%%%%%%%%
% figure()
% plot(timidivi_timestep, timidivi_col_S, 'b--o','LineWidth',2);
% hold on
% plot(timidivi_timestep, timidivi_col_I,'c--o','LineWidth',2);
% plot(timidivi_timestep, timidivi_col_R, 'r--o','LineWidth',2);
% legend('Susceptible_{RRMSE}','Infected_{RRMS}','Recovered_{RRMSE}');
% title("SIR Petri Net vs ODE RRMSE based on times steps", 'FontSize', 24)
% subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' sigma='+sigma*timidivi)
% xlabel('PN Time steps (per 1 ODE time interval)'); ylabel('RRMSE');
% hold off


%disp(['still going: ', int2str(sigmactr)]);

%% Plot the results of single param value run vs ODE %%%%%%%%%%%%%%%%%%
% %  figure()
% % plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible, 'b--o','LineWidth',2);
% % hold on
% % plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MInfected,'c--o','LineWidth',2);
% % plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MRecovered, 'r--o','LineWidth',2);
% % 
% % plot(t,y(:,1), 'y','LineWidth',2)
% % plot(t,y(:,2),'m','LineWidth',2)
% % plot(t,y(:,3), 'k','LineWidth',2);
% % legend('Susceptible_{PN}','Infected_{PN}','Recovered_{PN}',...
% %     'Susceptible_{ODE}','Infected_{ODE}','Recovered_{ODE}');
% % title("SIR Petri Net vs ODE", 'FontSize', 24)
% % subtitle("beta="+beta*timidivi+', gamma='+gamma*timidivi+' sigma='+sigma*timidivi+", Susc. MSE="+err_S...
% %    +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
% % hold off
% 
% % % % % figure()
% % % % % plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MSusceptible/popscalar, 'b--o','LineWidth',4);
% % % % % hold on
% % % % % plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MExposed/popscalar, 'k--o','LineWidth',4);
% % % % % plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MInfected/popscalar,'c--o','LineWidth',4);
% % % % % plot(((1/timidivi):(1/timidivi):(MAX_ITERATIONS/timidivi))+(1-(1/timidivi)), MRecovered/popscalar, 'r--o','LineWidth',4);
% % % % % % plot([1,1+1/timidivi,2:1:100], MSusceptible, 'b--o','LineWidth',4);
% % % % % % hold on
% % % % % % plot([1,1+1/timidivi,2:1:100], MExposed, 'k--o','LineWidth',4);
% % % % % % plot([1,1+1/timidivi,2:1:100], MInfected,'c--o','LineWidth',4);
% % % % % % plot([1,1+1/timidivi,2:1:100], MRecovered, 'r--o','LineWidth',4);
% % % % % plot(t,y(:,1), 'b','LineWidth',2)
% % % % % plot(t,y(:,2),'k','LineWidth',2)
% % % % % plot(t,y(:,3), 'c','LineWidth',2);
% % % % % plot(t,y(:,4), 'r','LineWidth',2);
% % % % % legend('Susceptible_{PN}','Exposed_{PN}','Infected_{PN}','Recovered_{PN}',...
% % % % %     'Susceptible_{ODE}','Exposed_{ODE}','Infected_{ODE}','Recovered_{ODE}');
% % % % % title("NL SEIR Petri Net vs ODE", 'FontSize', 24)
% % % % % subtitle("S rrmse="+err_S+', E rrmse='+err_E+', I rrmse='+err_I+", R rrmse="+err_R, 'FontSize', 14)
% % % % % %   +", Infe. MSE="+err_I+", Reco. MSE="+err_R, 'FontSize', 14)
% % % % % hold off


        % MSE_s(gammactr,betactr)=err_S;
        % MSE_i(gammactr,betactr)=err_I;
        % MSE_iBeta(gammactr,betactr)=beta;
        % MSE_iGamma(gammactr,betactr)=gamma;
        % MSE_r(gammactr,betactr)=err_R;

        MSE_s(gammactr,betactr,sigmactr)=err_S;
        MSE_e(gammactr,betactr,sigmactr)=err_E;
        MSE_i(gammactr,betactr,sigmactr)=err_I;
        MSE_r(gammactr,betactr,sigmactr)=err_R;
        % 
% MSE_S_c(sigmactr)=zeros(1,5);
% MSE_S_r(sigmactr)=zeros(1,5);
% MSE_S_v(sigmactr)=zeros(1,5);
% MSE_S_C(dr)=[MSE_S_C(dr) err_S];
% % 
           gammactr=gammactr+1;
         end
        betactr=betactr+1;
     end
     sigmactr=sigmactr+1;
  end
  
  
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
% xlabel("\tau (parameters \gamma,\beta,\sigma range [0,\tau])");ylabel('RRMSE Mean')
% subtitle("beta="+beta+', gamma='+gamma+' sigma='+sigma+", Susc. MSE="+err_S...
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

%% Plot the results of grid param value RRMSE vs ODE
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
tiledlayout(5,4);
nexttile
s=surf(spacer,spacer,MSE_s(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible sigma=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
%h.Label.String = "RRMSE (%)";h.Label.Rotation = 0;
view(2)

nexttile
s=surf(spacer,spacer,MSE_e(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Exposed sigma=0')
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
title('Infected sigma=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,1));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered sigma=0')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h, 'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible sigma=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_e(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Exposed sigma=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected sigma=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,2));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered sigma=0.001')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Susceptible sigma=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_e(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Exposed sigma=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Infected sigma=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,3));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE (%)');
title('Recovered sigma=0.01')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible sigma=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_e(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Exposed sigma=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected sigma=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,4));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered sigma=0.1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_s(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Susceptible sigma=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_e(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Exposed sigma=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_i(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Infected sigma=1')
s.FaceColor='interp';
s.EdgeColor='interp';
clim([casxismin,casxismax]);
colormap(Cus_map);
%h=colorbar;title(h,'RRMSE (%)');
view(2)

nexttile
s=surf(spacer,spacer,MSE_r(:,:,5));
ylabel('gamma'); xlabel('beta'); zlabel('RRMSE');
title('Recovered sigma=1')
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
MSE_e
MSE_i
MSE_r
ode_run_counter=ode_run_counter+1;
%end

toc
%% ODE Function
function dNLSEIR = NLSEIR(t,a,params)
S = a(1); E = a(2); I = a(3); R = a(4);

sigma=params(1); beta=params(2); gamma=params(3);
mu=params(4); alpha=params(5);

dS=mu-((beta*S*I)/(1+alpha*I*I))-mu*S;
dE=((beta*S*I)/(1+alpha*I*I))-sigma*E-mu*E;
dI=sigma*E-gamma*I-mu*I;
dR=gamma*I-mu*R;
dNLSEIR=[dS;dE;dI;dR];
end