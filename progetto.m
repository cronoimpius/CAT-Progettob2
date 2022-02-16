%% inizializzazione
clear all; close all; clc
omega_plot_min=1e-1;
omega_plot_max=1e4;
s=tf('s');
b1=0.3;
b2=0.1;
m=1;
k=1.5;
x3=3*10^7;
K=6.67*10^(-11);
M=5.98*10^24;
x2=0;
x1=sqrt((K*M)/(x3^3));
u=b2*sqrt((K*M)/x3);
A=[ -2*x2/x3-b2/m, -2*x1/x3, (2*x2*x1)/x3^2 - u/(m*x3^2); -(k-1)*2*x1*x3, -b1/m, (k-1)*(-2*K*M/x3^3 -x1^2); 0,1,0 ];
B=[ 1/(m*x3);0;0];
D=0;
C=[1,0,0];
num=(s+b1/m)*s+(k-1)*(2*K*M/x3^3 +x1^2);
den=s*(s+b2/m)*(s+b1/m)-u/(m*x3)*(k-1)*2*x1-4*x1^2*s*(k-1)+(k-1)*(2*K*M/x3^3 +x1^2)*(s+b2/m);
G=(num/(m*x3))/den;
step(G);
%% diagramma di bode
figure(1)
h_G=bodeplot(G);
grid on; zoom on;
%% iniziamo usando un regolatore statico con guadagno unitario
mu_s=1; % cambio il valore per rispettare le specifiche
R_s=mu_s/s; 
G_e=R_s*G
%%
% specifiche su d

if 0
    figure(2)
    h_Ge = bodeplot(G_e);
    grid on, zoom on;
    return
 end

figure(2)

omega_d_min=0.0001; %per fare bene il disegno 
omega_d_max=0.08; 
Ad=45; 
bnd_d_x =[omega_d_min; omega_d_max; omega_d_max; omega_d_min];
bnd_d_y= [Ad;Ad;-700;-700];
patch(bnd_d_x,bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

if 0
    h_Ge=bodeplot(G_e);
    grid on, zoom on;
    %return
end

% specifiche su n
omega_n_min=(5)*10^4;
omega_n_max=(7.5)*10^7;
An=-85;
bnd_n_x=[omega_n_min; omega_n_max; omega_n_max; omega_n_min]; 
bnd_n_y=[An; An; 1000; 1000];
patch(bnd_n_x,bnd_n_y,'b','FaceAlpha',0.2,'EdgeAlpha',0); 
hold on;

if 0
    h_Ge=bodeplot(G_e);
    grid on, zoom on;
    %return
end

% sovrelongazione
S_100_spec = 1/100;
xi = sqrt(log(S_100_spec)^2/(pi^2+log(S_100_spec)^2));
Mf=xi*100
%S_100 = 100*exp(-pi*xi/sqrt(1-xi^2))

% specifiche su Ta

Ta5_spec=0.15;
omega_Ta_low = 1e-4; % lower bound just for the plot
omega_Ta_MAX = 300/(Mf*Ta5_spec);

Bnd_Ta_x = [omega_Ta_low; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_low];
Bnd_Ta_y = [0; 0; -350; -350];
patch(Bnd_Ta_x, Bnd_Ta_y,'y','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

h_Ge = bodeplot(G_e);
grid on, zoom on;

% STOP qui per disegnare solo le specifiche sul GUADAGNO
if 0
  return;
end

%specifiche su Mf

omega_c_min=omega_Ta_MAX;
omega_c_max=omega_n_min;

phi_spec=Mf-180;
phi_low=-270; % per disegnare meglio

bnd_Mf_x=[omega_c_min; omega_c_max; omega_c_max; omega_c_min];
bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(bnd_Mf_x, bnd_Mf_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

%STOP qui per le specifiche statiche
if 0
  return 
end
%%
% Da quello che abbiamo visto sopra ci troviamo in uno scenario B, per cui:
% Progetto regolatore dinamico (Rete Anticipatrice)

Mf_star=Mf;

omega_c_star=omega_c_min+5;
%scegliamo una frequenza di taglio nell'intervallo [omega_c_min;
%omega_c_max]

[mag_omega_c_star, arg_omega_c_star, omega_c_star]=bode(G_e,omega_c_star);

mag_omega_c_star_dB=20*log10(mag_omega_c_star)

M_star=10^(-mag_omega_c_star_dB/20)
%scegliamo un phi_star>82.6 per avere più margine
phi_star= Mf_star - 180 - arg_omega_c_star+7
cos(phi_star)
%formule di inversione
alpha_tau=(cos(phi_star*pi/180)-1/(M_star))/(omega_c_star*sin(phi_star*pi/180))
tau=(M_star- cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180))

alpha=alpha_tau/tau

check = cos(phi_star*pi/180)-inv(M_star)
if check<0
    disp('Alpha negativo')
    return;
end

%%
R_d = (1+tau*s)/(1+alpha*tau*s);
LL=R_d*G_e;

figure(3)
hold on;
%specifiche su d

omega_d_min=0.0001; %per fare bene il disegno 
omega_d_max=0.08; 
Ad=45; 
bnd_d_x =[omega_d_min; omega_d_max; omega_d_max; omega_d_min];
bnd_d_y= [Ad;Ad;-350;-350];
patch(bnd_d_x,bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

%specifiche su n

omega_n_min=(5)*10^4;
omega_n_max=(7.5)*10^7;
An=-85;
bnd_n_x=[omega_n_min; omega_n_max; omega_n_max; omega_n_min]; 
bnd_n_y=[An; An; 100; 100];
patch(bnd_n_x,bnd_n_y,'b','FaceAlpha',0.2,'EdgeAlpha',0); 
hold on;

%specifiche su T_a

Bnd_Ta_x = [omega_Ta_low; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_low];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'y','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

h_Ge=bodeplot(G_e,{ 1e-1,1e7},'r');
h_LL= bodeplot(LL,{ 1e-1,1e7},'b');
grid on, zoom on;

margin(LL,{ omega_plot_min,omega_plot_max});

%specifiche su M_f

bnd_Mf_x=[omega_c_min; omega_c_max; omega_c_max; omega_c_min];
bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(bnd_Mf_x, bnd_Mf_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

%stop per vedere il diagramma
if 0
    return;
end

%%
%check delle caratteristiche in anello chiuso del sistema

F=LL/(1+LL);
figure(4);
%step(F)
W=1; %ampiezza del gradino
T_sim=5;

[y_step,t_step]=step(W*F,T_sim);
S_F=stepinfo(F,'SettlingTimeThreshold',0.05);
st_F=S_F.SettlingTime



%S_F=stepinfo(F);

plot(t_step,y_step)
grid on,zoom on, hold on;
%sovraelongazione

patch([0,T_sim,T_sim,0],[W*(1+S_100_spec),W*(1+S_100_spec),W+1,W+1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
ylim([0,W+1]);
%aggiungo il vincolo sul tempo di assestamento


patch([Ta5_spec,T_sim,T_sim,Ta5_spec],[W*(1-0.05),W*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([Ta5_spec,T_sim,T_sim,Ta5_spec],[W*(1+0.05),W*(1+0.05),W+1,W+1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

figure(5);
bode(F)
grid on, zoom on

if 0 
    return;
end
%% utilities per simulink
%open('progetto.slx');
x=[1,2,3,4];
WW=8*10^(-5);
DD=3*10^(-5);
NN=2*10^(-4);
%parte sistema lineare
R=R_d*R_s;
[nuR,deR]=tfdata(R);
numR=nuR{1};
denR=deR{1};
[nuG,deG]=tfdata(G);
numG=nuG{1};
denG=deG{1};
%parte sistema non lineare
x01=sqrt((K*M)/(x3^3));
x02=0;
x03=3*10^7;
ye=x01;
ue=b2*sqrt((K*M)/x03);