% These simulatiions show fraction of TCR that are inside the CCZ after 2s
% in the prcense of pMHC. Note that the simulations are done under the
% assumtion that CCZ exist for only 2 seconds.

function Figure7_ContactTimesVsProbTrig(panel)



global tao kon pMHC koff % Parameters used in ODE function

tmin=2;
D=0.05; % Diffusion rate (measured)
konV=[0 0.01];%.01; % = Kon = 2D association rate (µm2/s) -- value you suggested - reference?
koff=0.5;%10; %0.5 agonist, 10 non-agonist value you suggested - reference?
pMHC=30;%300; %30 for koff = 0.5, Density of pMHC (number of pMHC/µm2) - value commonly used in LIT according to Ricardo
TCR=(15000/415)*0.37;  % Density of TCR (number of TCRs/µm2), taking segregation into account
% (0.37 of all TCRs found in contact)
tot_timeR=logspace(0.5,8); % Different time points - wide range because of comparison to LIT
tot_timeR=[2 tot_timeR];
switch panel
    
    case 'a'
        %% Figure panel a: compare the probability of triggering in absence of pMHC (LIT only) for two radii as explained below
        
        
        %% diameter of real cell-cell contact measured to be 0.43 (+/- 0.034) um (Springer AH et al. (2012), J Immunol 188: 3686–3699.)
        %% radius we measured on glass is 0.385um (and there we saw LIT)- the idea is to compare the situation at those two values
        kon = konV(1);
        r=[0.215 0.3225];
        area=r.^2.*pi; % Area of ccz
        numb_TCRs=TCR.*area; % Total number of TCR that an exist in a ccz
        numb_TCRs=numb_TCRs;%*10; % Assum that 10 ccz exist in a "real" cell-to-cell contact
        taoV=r.^2/(8*D); % Mean residence time (In the ODE I use 1/tao to get the rate)
        
        
        for i=1:length(tot_timeR) % Loop over different Koffs
            for ii=1:length(r)
                tao = taoV(ii);
                n0=zeros(3,1); % Initial values for ODEs
                n0(1)=TCR;%  TCR; % Initial number of TCRs in ccz
                
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
                p_LIT(ii)=(AUC_tot-AUC_tmin)/AUC_tot;% Probability that a TCR stay longer than 2s in the ccz in presence of pMHC
                
                trailsR=round((tot_timeR(i)./tao).*numb_TCRs(ii)); % Total number of TCRs that will exist in the cczs during the cell-to-cell contact
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                P(i,ii)=1 - binopdf(0,trailsR,p_LIT(ii)); % % Prob that at least one TCR triggers (stay longer than 2s) (pMHC)
                
                
            end
        end
        
        %% calculating the x value (contact time) for y (prob) = 0.5:
        numtrialsreq = log(0.5)./log(1-p_LIT);
        timereq = (numtrialsreq.*taoV)./numb_TCRs;
        
        
        
        figure
        semilogx(tot_timeR,P(:,1),'b')
        hold on
        semilogx(tot_timeR,P(:,2),'k')
        %plot([timereq(1) timereq(1)],[0 1],'r')
        plot([2 2],[0 1],'r')
        ylabel('Probability that at least one TCR stay > 2s in a ccz during a cell-to-cell contact')
        xlabel('Time of cell-to-cell contact (s)')
        axis([0.1 10^10 0 1])
        legend(strcat('LIT 0.215, t = ', num2str(timereq(1)/(60*60)), ' h'),strcat('LIT 0.32, t = ', num2str(timereq(2)),' s'));
        title('Figure 7 a')
        
        
        
    case 'b'
        %% Figure panel b: compare the probability of triggering in the presence of pMHC to LIT for contacts of the size expected in vivo
        
        r=0.215; %diameter measured to be 0.43 +/- 0.034 um (Springer AH et al. (2012), J Immunol 188: 3686–3699.)
        area=r.^2.*pi; % Area of ccz
        tao=r^2/(8*D); % Mean residence time (In the ODE I use 1/tao to get the rate
        numb_TCRs=TCR.*area; % Total number of TCR that an exist in a ccz
        numb_TCRs=numb_TCRs;%*10; % Assum that 10 ccz exist in a "real" cell-to-cell contact
        
        % Loop over different Kons
        for ii=1:length(konV)
            kon = konV(ii);
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
            p(ii)=(AUC_tot-AUC_tmin)/AUC_tot;% Probability that a TCR stay longer than 2s in the ccz in presence of pMHC
            
            for i=1:length(tot_timeR)
                trailsR=round((tot_timeR(i)./tao).*numb_TCRs); % Total number of TCRs that will exist in the cczs during the cell-to-cell contact
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                P(i,ii)=1 - binopdf(0,trailsR,p(ii)); % % Prob that at least one TCR triggers (stay longer than 2s)
            end
            
        end
        
        %% calculating the x value (contact time) for y (prob) = 0.5:
        numtrialsreq = log(0.5)./log(1-p);
        timereq = (numtrialsreq.*tao)./numb_TCRs;
        
        figure
        semilogx(tot_timeR,P(:,2),'b')
        hold on
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %LIT case
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        semilogx(tot_timeR,P(:,1),'k')
         plot([2 2],[0 1],'r')
        ylabel('Probability that at least one TCR stay > 2s in a ccz during a cell-to-cell contact')
        xlabel('Time of cell-to-cell contact (s)')
        
        legend(strcat('pMHC, t = ', num2str(timereq(2)), ' s'),strcat('LIT, t = ', num2str(timereq(1)/(60*60)),' h'));
        axis([1 10^10 0 1])
        title('Figure 7 b')
end

% Function for the ordinary differential equations (eq 9 and 10 in supp info by Andreas in January)
function dndt=func(t,n);
global tao kon pMHC koff
dndt=zeros(3,1); %
dndt(1)=-kon*n(1)*(pMHC-n(2))+koff*n(2)-(1/tao)*n(1);
dndt(2)=kon*n(1)*(pMHC-n(2))-koff*n(2);
dndt(3)=n(1); % Area under curve