SNRdB=[10 20 30]; 
R_thresholds = [1,1.5,2:0.1:6,6.5,7]; 
Feasibility_RS = zeros(length(SNRdB),length(R_thresholds));
Feasibility_MULP = zeros(length(SNRdB),length(R_thresholds));

save('Feasibility_Jun1.mat','Feasibility_RS','Feasibility_MULP')


