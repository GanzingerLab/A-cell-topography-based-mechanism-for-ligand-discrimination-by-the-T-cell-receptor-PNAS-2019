close all;

%%plot the cummulative probability distribution

load CaptureDataSelfPeptideR0.mat

figure(5); set(gcf,'color','w'); box on; hold on;
plot(t,C,'linewidth',2); 



%load CaptureDataSelfPeptide2R0.mat

%plot(t,C,'linewidth',2); 

load CaptureDataAgonistR0.mat
plot(t,C,'linewidth',2); 


%load CaptureDataAgonist2R0.mat

%plot(t,C,'linewidth',2);

load CaptureDataNopMHCR0.mat

plot(t,C,'linewidth',2); 



%load CaptureDataNopMHC2R0.mat

%plot(t,C,'linewidth',2);

legend('Self Peptide R_0 = 0.22','Agonist R_0 = 0.22','No pMHC, R_0 = 0.22');
title('Occupation time distribution'); xlabel('time [s]'); ylabel('Probability Density');
xlim([01 3.0]); hold off;