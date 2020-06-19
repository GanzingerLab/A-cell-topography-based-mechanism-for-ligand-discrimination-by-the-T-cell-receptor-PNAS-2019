% These simulatiions show fraction of TCR that are inside the CCZ after 2s
% in the prcense of pMHC. Note that the simulations are done under the
% assumtion that CCZ exist for only 2 seconds.

function [xy4plot,trigg] = Figure7_konVsHalfMaxTrig_different_koffs

global tao kon pMHC koff % Parameters used in ODE function

tmin=2; % Min to stay in CCZ to trigger
r=0.215; % Radius of CCZ
D=0.05; % Diffusion rate
tao=r^2/(8*D) % Mean residence time (In the ODE I use 1/tao to get the rate)
area=r.^2.*pi; % Area of ccz
pMHC=30; % Density of pMHC (number of pMHC/µm2)
TCR=100*0.37;  % Density of TCR (number of TCRs/µm2)
numb_TCRs=TCR.*area; % Total number of TCR that an exist in a ccz
numb_TCRs=numb_TCRs*10; % Assum that 10 ccz exist in a "real" cell-to-cell contact



koffV=1:1:10; % = Kon = 2D association rate (µm2/s)
konV=logspace(-6,-2)
count = 0;

for i=1:length(koffV)
    
    for ii = 1:length(konV)
        count = count + 1;
        kon=konV(ii); % = Kon = 2D association rate (µm2/s) -- value you suggested - reference?
        koff=koffV(i); % value you suggested - reference?
        
        xy4plot(count,1:2) = [kon koff];
        
        pMHC=30; % Density of pMHC (number of pMHC/µm2) - value commonly used in LIT according to Ricardo
        
        
        n0=zeros(3,1); % Initial values for ODEs
        n0(1)=numb_TCRs; % Initial number of TCRs in ccz
        
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
        p_pMHC(ii)=(AUC_tot-AUC_tmin)/AUC_tot;
        % pLIT(ii)=exp(-tmin./tao); % Probability that a TCR stay longer than 2s in the ccz in absence of pMHC
        
        
           
        
        %% The idea behind the following lines is that this should give us a means
        %% to compare curves of contact times vs probability of at least a single TCR staying in the CCZ > 2s
        %% for multiple radii: how long would a contact need to last for (at a given radius) in order for the
        %% probability of at least a single TCR staying in the CCZ > 2s to be 0.5?
        
        numtrialsreq(count) = log(0.5)./log(1-p_pMHC(ii)); %I found this by calculating P = 1-p(no successes) = 1-(1-p(success))^n
        %and rearranging for n: n trials required for P(at least one success) =0.5 given
        %the probability of a success is p (= p_pMHC,pLIT -- your AUC)
        timereq(count) = (numtrialsreq(ii)*tao)./numb_TCRs;
        
        %numtrialsreqLIT(ii) = log(0.5)./log(1-pLIT(ii));
        %timereqLIT(ii) = (numtrialsreqLIT(ii)*tao)./numb_TCRs;
        
        
    end
    
    
end

trigg = timereq<120;

    figure(1)
        
        semilogy(konV,timereq(1:50),'k')
        
       % scatter3(xy4plot(1,:),xy4plot(2,:),trigg')
        %xlim([0 2])
        ylabel('Time of cell-to-cell contact (s) required for p(singleTCR>2s) = 0.5')
        xlabel('k_{on} (s^{-1}) of TCR/pMHC complex')
    
end

% Function for the ordinary differential equations (eq 9 and 10 in supp info by Andreas in January)
function dndt=func(t,n);
global tao kon pMHC koff
dndt=zeros(3,1); %
dndt(1)=-kon*n(1)*(pMHC-n(2))+koff*n(2)-(1/tao)*n(1);
dndt(2)=kon*n(1)*(pMHC-n(2))-koff*n(2);
dndt(3)=n(1); % Area under curve
end