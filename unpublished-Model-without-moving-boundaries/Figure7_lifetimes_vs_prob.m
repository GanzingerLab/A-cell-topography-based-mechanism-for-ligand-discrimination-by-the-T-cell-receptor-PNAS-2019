% These simulatiions calculates AUC to get probability that a TCR stay
% longer than 2 sec where we can change koff, kon, pMHC etc.

function Figure7_lifetimes_vs_prob
%close all; %Do used because we want to plot several curves in the same plot

global tao kon pMHC koff % Parameters used in ODE function
tmin=2; % Min to stay in CCZ to trigger
r=0.215; % Radius of CCZ
D=0.05; % Diffusion rate
tao=r^2/(8*D) % Mean residence time (In the ODE I use 1/tao to get the rate)
area=r.^2.*pi; % Area of ccz
pMHC=30; % Density of pMHC (number of pMHC/µm2)
TCR=100*0.37;  % Density of TCR (number of TCRs/µm2)
numb_TCRs=TCR.*area; % Total number of TCR that an exist in a ccz
numb_TCRs=numb_TCRs%*10; % Assum that 10 ccz exist in a "real" cell-to-cell contact


konV=[0.01 0.1 1 10 100]; % = Kon = 2D association rate (µm2/s) 
param=logspace(0,4); % Different values for a parameter to be changed (0.01 - 1000) 

for ii=1:length(konV)
    kon = konV(ii);
for i=1:length(param) % Loop over different Koffs
koff=param(i); % Change Koff during the loop (of another parameter if you like)

K(i)=1/koff; % Convert to Life time to be used as x-axis in the figure

n0=zeros(3,1); % Initial values for ODEs
n0(1)=TCR; % Initial number of TCRs in ccz

tmax=1e6; % Run until all TCR are depleated from the CCZ (1 million s should be enough)
[t,n]=ode15s(@func,[0,tmax],n0); %
N1=n(:,1); % Density of Unbound TCR over time 
N2=n(:,2); % Density of bound TCR over time
AUC_tot=n(end,3); % AUC for total area

[t1,n]=ode15s(@func,[0,tmin],n0); % Run to tmin s
N11=n(:,1); % Density of Unbound TCR over time 
N22=n(:,2); % Density of bound TCR over time
AUC_tmin=n(end,3); % AUC for total area until 2 s

% Calculate fraction of the area that is on the right hand side of tmin
probability=(AUC_tot-AUC_tmin)/AUC_tot;

pp(i)=probability; % Probability to stay longer than 2 s
trails=round((60./tao).*numb_TCRs); % Trails assuming CCZ lasts for 1 minute
pr2(i)=1 - binopdf(0,trails,pp(i)); % Prob that at least one TCR triggers

end

% figure(1)
% semilogx(K,pp)
% hold on
% ylabel('Probability that a TCR stay longer than 2 seconds')
% xlabel('Life-time (s) of TCR/pMHC compex')
%axis([0.01 100 0 1])

figure(2)
semilogx(K,pr2)
hold on
ylabel('Probability that at least one TCR stay longer than 2 seconds')
xlabel('Life-time (s) of TCR/pMHC compex')
%axis([0.01 100 0 1])
end
legend('k_{on} = 0.01','k_{on} = 0.1','k_{on} = 1','k_{on} = 10','k_{on} = 100')


% Function for the ordinary differential equations (eq 9 and 10 in supp info by Andreas in January)
function dndt=func(t,n);
global tao kon pMHC koff
dndt=zeros(3,1); %
dndt(1)=-kon*n(1)*(pMHC-n(2))+koff*n(2)-(1/tao)*n(1); 
dndt(2)=kon*n(1)*(pMHC-n(2))-koff*n(2);
dndt(3)=n(1); % Area under curve


