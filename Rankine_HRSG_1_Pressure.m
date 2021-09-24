%CICLO RANKINE Y HRSG DE UN NIVEL DE PRESIÓN
%Variables de entrada del ciclo Rankine y del HRSG
p1_H2O=10000000;
T1_H2O=565+273;
% DATOS
PP_min=5;
T7_gases_min=90+273;
eta_tv_iso=0.9;
eta_b_iso=0.85;


% (1) A la entrada de la turbina de vapor.
h1_H2O=CoolProp.PropsSI ('H', 'P', p1_H2O, 'T', T1_H2O, 'Water');
s1_H2O=CoolProp.PropsSI ('S', 'P', p1_H2O, 'T', T1_H2O, 'Water');
% (2) Entre la turbina de vapor y el condensador.

T2_H2O=T_amb+15;
p2_H2O=CoolProp.PropsSI ('P', 'T', T2_H2O, 'Q', 0, 'Water');
s2i_H2O=s1_H2O;
h2i_H2O=CoolProp.PropsSI ('H', 'P', p2_H2O, 'S', s2i_H2O, 'Water');
h2_H2O=h1_H2O-eta_tv_iso*(h1_H2O-h2i_H2O);
Q2_H2O=CoolProp.PropsSI ('Q', 'P', p2_H2O, 'H', h2_H2O, 'Water');

% (3) Entre el condensador y la bomba de alimentación.
p3_H2O=p2_H2O;
T3_H2O=CoolProp.PropsSI ('T', 'P', p3_H2O, 'Q', 0, 'Water');
h3_H2O=CoolProp.PropsSI ('H', 'P', p3_H2O, 'Q', 0, 'Water');
s3_H2O=CoolProp.PropsSI ('S', 'P', p3_H2O, 'Q', 0, 'Water');
% (4) Entre la bomba y la entrada de la caldera de recuperación de
calor(HRSG).
p4_H2O=p1_H2O;
s4i_H2O=s3_H2O;
h4i_H2O=CoolProp.PropsSI ('H', 'S', s4i_H2O, 'P', p4_H2O, 'Water');
h4_H2O=h3_H2O+(h4i_H2O-h3_H2O)/eta_b_iso;
T4_H2O=CoolProp.PropsSI ('T', 'P', p4_H2O, 'H', h4_H2O, 'Water');

%% HRSG
p4_H2O=p1_H2O;
p5_H2O=p1_H2O;
p6_H2O=p1_H2O;
% (5) Entre el economizador y el calderín.
T5_H2O=CoolProp.PropsSI ('T', 'P', p5_H2O, 'Q', 0, 'Water');
h5_H2O=CoolProp.PropsSI ('H', 'P', p5_H2O, 'Q', 0, 'Water');
% (6) Entre el calderín y el sobrecalentador.
T6_H2O=CoolProp.PropsSI ('T', 'P', p6_H2O, 'Q', 1, 'Water');
h6_H2O=CoolProp.PropsSI ('H', 'P', p6_H2O, 'Q', 1, 'Water');
%%
% Se llevará a cabo un proceso iterativo para el cálculo de la
%relación de vapor de agua y aire (W). Para ello, se inicializará a un
%valor mínimo y evaluaremos para cada iteración si cumple los
%requisitos termodinámicos impuestos:
% -Temperatura mínima de Pinch Point de 5ºC.
% -Temperatura mínima de los gases a la salida de la caldera de
% recuperación de 90ºC
i=1;
N=1000;
W=linspace(0.001,1,1000);
 for i=1:N
 % BALANCES DE ENERGÍA
 % Sobrecalentador
 T5_gases=T4_gases-W(i)/((1+F)*Cp_gases)*(h1_H2O-h6_H2O);
 % Evaporador

 T6_gases=T5_gases-W(i)/((1+F)*Cp_gases)*(h6_H2O-h5_H2O);

 % Economizador

 T7_gases=T6_gases-W(i)/((1+F)*Cp_gases)*(h5_H2O-h4_H2O);
 % Pinch Point
 PP=T6_gases-T5_H2O;
 if (T7_gases>T7_gases_min) && (PP>PP_min)&&
(T4_gases>T5_gases) && (T5_gases>T6_gases) &&
(T6_gases>T7_gases)
 W_sol(i)=W(i);
 T7_gases_sol(i)=T7_gases;
 PP_sol(i)=PP;
 Q_HRSG_sol(i)=F*LHV-(W_tg-W_c);
 Else
 W_sol(i)=NaN;
 T7_gases_sol(i)=NaN;
 PP_sol(i)=NaN;
 Q_HRSG_sol(i)=NaN;
 end
 end
%%
% Calores
 for i=1:N

Q_HRSG_sol(i)=F*LHV-(W_tg-W_c);
 Q_Rankine_sol(i)=(1+F)*Cp_gases*(T4_gases-T7_gases_sol(i));
 Q_perdidas_sol(i)=Q_HRSG_sol(i)-Q_Rankine_sol(i);
 eta_HRSG_sol(i)=Q_Rankine_sol(i)/Q_HRSG_sol(i);
 W_tv_sol(i)=W_sol(i)*(h1_H2O-h2_H2O);
 W_b_sol(i)=W_sol(i)*(h4_H2O-h3_H2O);
 W_net_Rankine_sol(i)=W_tv_sol(i)-W_b_sol(i);
 eta_Rankine_sol(i)=W_net_Rankine_sol(i)/Q_Rankine_sol(i);


eta_COMBINADO_sol(i)=eta_Brayton+eta_Rankine_sol(i)*eta_HRSG_sol
(i)*(1-eta_Brayton);
 end