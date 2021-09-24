%%
%CICLO BRAYTON
%VARIABLES DEL CICLO BRAYTON
PR=15;
T3_gases=1400+273;
% DATOS
p_amb=101325;
T_amb=288;
p1=p_amb;
T1=T_amb;
p4_gases=101325;
Cp_aire=1006;
Cv_aire=717.462;
gamma_aire=Cp_aire/Cv_aire;
LHV=50047000;
Cp_gases=1100;
gamma_gases=1.4;
Cp_fuel=2253.7;
T_fuel=25+273;
eta_c_iso=0.89;
eta_tg_iso=0.91;
eta_cc=0.99;
% Compresor.
p2=PR*p1;
T2=T1*(1+(1/eta_c_iso)*(PR^((gamma_aire-1)/gamma_aire)-1));
W_c=Cp_aire*(T2-T1);
% Cámara de combustión.
F=(Cp_gases*T3_gases-Cp_aire*T2)/(eta_cc*LHVCp_gases*T3_gases+Cp_fuel*T_fuel);
Q_cc=F*LHV;
% Turbina de gas.
p3_gases=p2;
T4_gases=T3_gases*(1-eta_tg_iso*(1-(p4_gases/p3_gases)^((gamma_gases1)/gamma_gases)));
W_tg=(1+F)*Cp_gases*(T3_gases-T4_gases);
W_net_Brayton=W_tg-W_c;
eta_Brayton=(W_tg-W_c)/Q_cc;