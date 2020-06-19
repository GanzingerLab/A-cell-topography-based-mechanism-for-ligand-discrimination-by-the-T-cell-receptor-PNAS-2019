% These simulatiions show fraction of TCR that are inside the CCZ after 2s
% in the prcense of pMHC. Note that the simulations are done under the
% assumtion that CCZ exist for only 2 seconds.

function [rV p_pMHC pLIT] = Figure7_radVsContactTimes

global tao kon pMHC koff % Parameters used in ODE function
rV = 0.15:0.01:2;
for rad = 1:length(rV)
    
    r = rV(rad);
    tmin=2;
    D=0.05; % Diffusion rate (measured)
    area=r.^2.*pi; % Area of ccz
    tao=r^2/(8*D); % Mean residence time (In the ODE I use 1/tao to get the rate)
    kon=0.01; % = Kon = 2D association rate (µm2/s) -- value you suggested - reference?
    koff=0.5; % value you suggested - reference?
    pMHC=30; % Density of pMHC (number of pMHC/µm2) - value commonly used in LIT according to Ricardo
   
    %% Density of TCR (number of TCRs/µm2), taking segregation into account
    TCR=100*0.37;  % (0.37 of all TCRs found in contact)
    numb_TCRs=TCR.*area; % Total number of TCR that an exist in a ccz
    numb_TCRs=numb_TCRs;%*10; % Assum that 10 ccz exist in a "real" cell-to-cell contact
    tot_timeR=logspace(-5,15); % Different time points - wide range because of comparison to LIT
    
    
    n0=zeros(3,1); % Initial values for ODEs
    n0(1)=TCR; % Initial number of TCRs in ccz
    
    tmax=1e6; % Run until all TCR are depleated from the CCZ
    [t,n]=ode15s(@func,[0,tmax],n0); %
    N1=n(:,1); % Density of Unbound TCR over time
    N2=n(:,2); % Density of bound TCR over time
    AUC_tot=n(end,3); % AUC for total area
    
    [t1,n]=ode15s(@func,[0,tmin],n0); % Run to tmin s
    N11=n(:,1); % Density of Unbound TCR over time
    N22=n(:,2); % Density of bound TCR over time
    AUC_tmin=n(end,3); % AUC for total area until 2 s
    
    % Calculate fraction of the area that is on the right hand side if tmin
    p_pMHC(rad)=(AUC_tot-AUC_tmin)/AUC_tot;
    pLIT(rad)=exp(-tmin./tao); % Probability that a TCR stay longer than 2s in the ccz in absence of pMHC
    
    
    if r<0.16
        for i=1:length(tot_timeR)
            
            trailsR=round((tot_timeR(i)./tao).*numb_TCRs); % Total number of TCRs that will exist in the cczs during the cell-to-cell contact
            prR(i)=1 - binopdf(0,trailsR,p_pMHC(rad)); % % Prob that at least one TCR triggers (stay longer than 2s)
            
        end
    end
    
    
    %% The idea behind the following lines is that this should give us a means
    %% to compare curves of contact times vs probability of at least a single TCR staying in the CCZ > 2s
    %% for multiple radii: how long would a contact need to last for (at a given radius) in order for the
    %% probability of at least a single TCR staying in the CCZ > 2s to be 0.5?
    
    numtrialsreq(rad) = log(0.5)./log(1-p_pMHC(rad)); %I found this by calculating P = 1-p(no successes) = 1-(1-p(success))^n 
                                                      %and rearranging for n: n trials required for P(at least one success) =0.5 given 
                                                      %the probability of a success is p (= p_pMHC,pLIT -- your AUC) 
    timereq(rad) = (numtrialsreq(rad)*tao)./numb_TCRs;
    
    numtrialsreqLIT(rad) = log(0.5)./log(1-pLIT(rad));
    timereqLIT(rad) = (numtrialsreqLIT(rad)*tao)./numb_TCRs;
end
figure

semilogy(rV,timereqLIT,'k')
hold on
semilogy(rV,timereq,'r')
xlim([0 2])
ylabel('Time of cell-to-cell contact (s) required for p(singleTCR>2s) = 0.5')
xlabel('radius (um)')
legend('LIT','pMHC')


% Function for the ordinary differential equations (eq 9 and 10 in supp info by Andreas in January)
function dndt=func(t,n);
global tao kon pMHC koff
dndt=zeros(3,1); %
dndt(1)=-kon*n(1)*(pMHC-n(2))+koff*n(2)-(1/tao)*n(1);
dndt(2)=kon*n(1)*(pMHC-n(2))-koff*n(2);
dndt(3)=n(1); % Area under curve