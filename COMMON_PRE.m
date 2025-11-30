% Common_Pre for SIR_v1
% sets firing conditions
function [fire, transition] = COMMON_PRE(transition)
global global_info m0Susceptible m0Infected m0Recovered beta
%global beta m0Susceptible m0Infected

% in each itertaion, all three transitions are allowed 
%      to fire only once
fire = not(timesfired(transition.name));


ntI = timesfired('tInfect'); 
ntR = timesfired('tRecover'); 
ntS = timesfired('tSuscept');


% if ntI==10 && ntR==10 && ntS==10
%     global_info.STOP_SIMULATION = and(ntI, and(ntR, ntS));
% end

%stop the run after the three transitions have fired once
global_info.STOP_SIMULATION = and(ntI, and(ntR, ntS));

% if beta*m0Susceptible*m0Infected<=m0Susceptible
% % stop the run after the three transitions have fired once
% global_info.STOP_SIMULATION = and(ntI, and(ntR, ntS));
% else
% % stop the run after the three transitions have fired twice
% global_info.STOP_SIMULATION = and((ntI-1), and((ntR-1), (ntS-1)));
% end