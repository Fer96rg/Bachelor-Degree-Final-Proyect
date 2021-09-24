%CICLO RANKINE Y HRSG DE DOS NIVELES DE PRESIÓN
%Variables de entrada
T1a_H2O=565+273;
T1b_H2O=400+273;
p1a_H2O=10000000;
p1b_H2O=800000;
%DATOS
PP_alta_min=10;
PP_baja_min=10;
T10_gases_min=90+273;
eta_tv_iso=0.9;
eta_b1_iso=0.85;
eta_b2_iso=0.85;
% (1) A la entrada de la turbina de vapor de cada nivel de presión.
h1a_H2O=CoolProp.PropsSI ('H', 'P', p1a_H2O, 'T', T1a_H2O, 'Water');
s1a_H2O=CoolProp.PropsSI ('S', 'P', p1a_H2O, 'T', T1a_H2O, 'Water');

h1b_H2O=CoolProp.PropsSI ('H', 'P', p1b_H2O, 'T', T1b_H2O, 'Water');
s1b_H2O=CoolProp.PropsSI ('S', 'P', p1b_H2O, 'T', T1b_H2O, 'Water');

% (2) Entre la turbina de vapor y el condensador.
T2_H2O=T_amb+15;
p2_H2O=CoolProp.PropsSI ('P', 'T', T2_H2O, 'Q', 0, 'Water');

% (1a) Vapor procedente del nivel de baja presión.
s2ai_H2O=s1a_H2O;
h2ai_H2O=CoolProp.PropsSI ('H', 'P', p2_H2O, 'S', s2ai_H2O, 'Water');
h2a_H2O=h1a_H2O-eta_tv_iso*(h1a_H2O-h2ai_H2O);
% (1b) Vapor procedente del nivel de alta presión.
s2bi_H2O=s1b_H2O;
h2bi_H2O=CoolProp.PropsSI ('H', 'P', p2_H2O, 'S', s2bi_H2O, 'Water');
h2b_H2O=h1b_H2O-eta_tv_iso*(h1b_H2O-h2bi_H2O);

% (3) Entre el condensador y la bomba de alimentación.
p3_H2O=p2_H2O;
T3_H2O=CoolProp.PropsSI ('T', 'P', p3_H2O, 'Q', 0, 'Water');
h3_H2O=CoolProp.PropsSI ('H', 'P', p3_H2O, 'Q', 0, 'Water');
s3_H2O=CoolProp.PropsSI ('S', 'P', p3_H2O, 'Q', 0, 'Water');

% (4) Entre la bomba y la entrada de la caldera de recuperación de
calor (HRSG)
p4_H2O=p1a_H2O;
s4i_H2O=s3_H2O;
h4i_H2O=CoolProp.PropsSI ('H', 'S', s4i_H2O, 'P', p4_H2O, 'Water');
h4_H2O=h3_H2O+(h4i_H2O-h3_H2O)/eta_b1_iso;
T4_H2O=CoolProp.PropsSI ('T', 'P', p4_H2O, 'H', h4_H2O, 'Water');

%% HRSG
p4_H2O=p1b_H2O;
p5_H2O=p1b_H2O;
p6_H2O=p1b_H2O;

p7_H2O=p1a_H2O;
p8_H2O=p1a_H2O;
p9_H2O=p1a_H2O;
% (5) Entre el economizador y el calderín del nivel de baja presión.
T5_H2O=CoolProp.PropsSI ('T', 'P', p5_H2O, 'Q', 0, 'Water');
h5_H2O=CoolProp.PropsSI ('H', 'P', p5_H2O, 'Q', 0, 'Water');
s5_H2O=CoolProp.PropsSI ('S', 'P', p5_H2O, 'Q', 0, 'Water');
% (6) Entre el calderín y el sobrecalentador del nivel de baja presión.
T6_H2O=CoolProp.PropsSI ('T', 'P', p6_H2O, 'Q', 1, 'Water');
h6_H2O=CoolProp.PropsSI ('H', 'P', p6_H2O, 'Q', 1, 'Water');

% (7) Entre la bomba de alimentación y el economizador del nivel de
alta presión.
s7i_H2O=s5_H2O;
h7i_H2O=CoolProp.PropsSI ('H', 'S', s7i_H2O, 'P', p7_H2O, 'Water');
h7_H2O=h5_H2O+(h7i_H2O-h5_H2O)/eta_b2_iso;
T7_H2O=CoolProp.PropsSI ('T', 'P', p7_H2O, 'H', h7_H2O, 'Water');
% (8) Entre el economizador y el calderín del nivel de alta presión.

T8_H2O=CoolProp.PropsSI ('T', 'P', p8_H2O, 'Q', 0, 'Water');
h8_H2O=CoolProp.PropsSI ('H', 'P', p8_H2O, 'Q', 0, 'Water');
% (9) Entre el calderín y el sobrecalentador del nivel de alta presión.

T9_H2O=CoolProp.PropsSI ('T', 'P', p9_H2O, 'Q', 1, 'Water');
h9_H2O=CoolProp.PropsSI ('H', 'P', p9_H2O, 'Q', 1, 'Water');
% Se llevará a cabo un proceso iterativo para el cálculo de la
%relación de vapor de agua y aire (W). Para ello, se inicializará a un
%valor mínimo y evaluaremos para cada iteración si cumple los
%requisitos termodinámicos impuestos:
% -Temperatura mínima de Pinch Point de 10ºC para los niveles de
%alta y baja.
% -Temperatura mínima de los gases a la salida de la caldera de
% recuperación de 90ºC.

N=99;
[alpha,W]=ndgrid(0.01:0.01:0.99,0.01:0.01:0.99);
for i=1:N
 for j=1:N
T5_gases=T4_gases-W(i,j)*(1-
alpha(i,j))/((1+F)*Cp_gases)*(h1a_H2O-h9_H2O);
T6_gases=T5_gases-W(i,j)*(1-
alpha(i,j))/((1+F)*Cp_gases)*(h9_H2O-h8_H2O);

T7_gases=T6_gases-W(i,j)*(1-
alpha(i,j))/((1+F)*Cp_gases)*(h8_H2O-h7_H2O);

T8_gases=T7_gasesW(i,j)*alpha(i,j)/((1+F)*Cp_gases)*(h1b_H2O-h6_H2O);

T9_gases=T8_gases1/((1+F)*Cp_gases)*(W(i,j)*alpha(i,j)*h6_H2O+W(i,j)*(1-
alpha(i,j))*h5_H2O-W(i,j)*h5_H2O);

 T10_gases=T9_gases-W(i,j)/((1+F)*Cp_gases)*(h5_H2O-h4_H2O);
 PP_alta=T9_gases-T5_H2O;
 PP_baja=T6_gases-T8_H2O;

if (T10_gases>T10_gases_min) && (PP_alta>PP_alta_min) &&
(PP_baja>PP_baja_min) && (T4_gases>T5_gases) &&
(T5_gases>T6_gases) && (T6_gases>T7_gases) &&
(T7_gases>T8_gases) && (T8_gases>T9_gases) &&
(T9_gases>T10_gases)
 W_sol(i,j)=W(i,j);
 alpha_sol(i,j)=alpha(i,j);
 T10_gases_sol(i,j)=T10_gases;
 PP_alta_sol(i,j)=PP_alta;
 PP_baja_sol(i,j)=PP_baja;
 else
 W_sol(i,j)=NaN;
 alpha_sol(i,j)=NaN;
 T10_gases_sol(i,j)=NaN;
 PP_alta_sol(i,j)=NaN;
 PP_baja_sol(i,j)=NaN;
 end
 end

end
% Calores
for j=1:N
 for i=1:N
 Q_HRSG_sol(i,j)=F*LHV-(W_tg-W_c);

 Q_Rankine_sol(i,j)=(1+F)*Cp_gases*(T4_gases-T10_gases_sol(i,j));
 Q_perdidas_sol(i,j)=Q_HRSG_sol(i,j)-Q_Rankine_sol(i,j);
 eta_HRSG_sol(i,j)=Q_Rankine_sol(i,j)/Q_HRSG_sol(i,j);

 W_tv_sol(i,j)=alpha_sol(i,j)*W_sol(i,j)*(h1b_H2O-h2b_H2O)+(1-
alpha_sol(i,j))*W_sol(i,j)*(h1a_H2O-h2a_H2O);
 W_b1_sol(i,j)=W_sol(i,j)*(h4_H2O-h3_H2O);
 W_b2_sol(i,j)=(1-alpha_sol(i,j))*W_sol(i,j)*(h7_H2O-h5_H2O);
 W_net_Rankine_sol(i,j)=W_tv_sol(i,j)-W_b1_sol(i,j)-W_b2_sol(i,j);
 eta_Rankine_sol(i,j)=W_net_Rankine_sol(i,j)/Q_Rankine_sol(i,j);

eta_COMBINADO_sol(i,j)=eta_Brayton+eta_Rankine_sol(i,j)*eta_HRSG
sol(i,j)*(1-eta_Brayton);
 end
end