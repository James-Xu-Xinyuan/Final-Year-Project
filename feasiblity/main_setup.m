SNRdB=[10 20 30]; 
R_thresholds = [1,1.5,2:0.1:6,6.5,7]; 
Feasibility_RS = zeros(length(SNRdB),length(R_thresholds));
Feasibility_MULP = zeros(length(SNRdB),length(R_thresholds));
save('Feasibility_Jun1.mat','Feasibility_RS','Feasibility_MULP')

R_thresholds = [1,1.5:0.1:5.5,6];  % 43 points
u2 = [10,3,1,1/3,0.1]; % 5 rate pairs
Feasibility_RS = zeros(length(u2),length(R_thresholds));
Feasibility_MULP = zeros(length(u2),length(R_thresholds));
save('Feasibility_Jun15_Weights.mat','Feasibility_RS','Feasibility_MULP')


