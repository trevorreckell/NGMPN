% SIR_v1 outlines and SIRS model with the addtional arc
% The _pdf file in general lays out all Places, transition, arc, and arc
% weights, but for variable pdf those need to be specified in main file
% file: ’simple_pn_SIRv1_pdf.m’
% Run this and COMMON_PRE.m before main file ’simple_pn_SIRv1_preprint.m’

function [pns] = simple_pn_SIRv1_pdf()
global global_info

pns.PN_name = 'SIR_v1';

%List of all Places
pns.set_of_Ps = {'pSusceptible','pInfected','pRecovered'};

%list of all Transitions
pns.set_of_Ts = {'tInfect','tRecover','tSuscept'}; 

%list of all Arcs/arc weights
pns.set_of_As = {...
 'pSusceptible','tInfect', global_info.pSus_tInf,... % tInfect
 'tInfect','pInfected', global_info.tInf_pInf,... % tInfect
 'pInfected','tRecover',global_info.pInf_tRec,... % tRecover
 'tRecover','pRecovered',global_info.tRec_pRec,... % tRecover
 'pRecovered','tSuscept',global_info.pRec_tSus,... % tSuscept 
 'tSuscept','pSusceptible',global_info.tSus_pSus ... % tSuscept 
  };    
