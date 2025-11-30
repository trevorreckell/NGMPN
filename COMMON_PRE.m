% NONlinear SEIR common_pre
function [fire, transition] = COMMON_PRE(transition)
global global_info

% % in each itertaion, all transitions are allowed 
% %      to fire only once
% fire = not(timesfired(transition.name));
% 
% nt_1= timesfired('t_1');
% nt_2= timesfired('t_2');
% nt_3= timesfired('t_3');
% nt_4= timesfired('t_4');
% nt_5= timesfired('t_5');
% nt_6= timesfired('t_6');
% nt_7= timesfired('t_7');
% nt_8= timesfired('t_8');
% 
% % stop the run after the transitions have fired once
% global_info.STOP_SIMULATION = and(nt_1, and(nt_2,and(nt_3,and(nt_4,and(nt_5,and(nt_6,and(nt_7,nt_8)))))));
fire = not(timesfired(transition.name));
end
