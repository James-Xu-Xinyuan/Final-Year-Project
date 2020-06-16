% load Feasibilities_EW.mat
R_thresholds = [1,1.5,2:0.1:6,6.5,7];
figure
plot(R_thresholds, Feasibility_MULP)
legend('MULP 10dB', 'MULP 20dB', 'MULP 30dB')
xlabel('R_th (bits/s/Hz)')
ylabel('feasibility')
title('Feasibility -  MULP - Equal Weight')

% load Feasibilities_EW.mat
R_thresholds = [1,1.5,2:0.1:6,6.5,7];
figure
plot(R_thresholds, Feasibility_RS)
legend('RS 10dB', 'RS 20dB', 'RS 30dB')
xlabel('R_th (bits/s/Hz)')
ylabel('feasibility')
title('Feasibility -  RS - Equal Weight')

% load Feasibility_3WeightPairs
R_thresholds = [1,1.5,2:0.1:6,6.5,7];
figure
plot(R_thresholds, Feasibility_MULP)
legend('MULP 3:1', 'MULP 1:1', 'MULP 1:3')
xlabel('R_th (bits/s/Hz)')
ylabel('feasibility')
title('Feasibility -  MULP - Variable Weight')
xlim([1 4])

% load Feasibility_3WeightPairs
R_thresholds = [1,1.5,2:0.1:6,6.5,7];
figure
plot(R_thresholds, Feasibility_RS)
legend('RS 3:1', 'RS 1:1', 'RS 1:3')
xlabel('R_th (bits/s/Hz)')
ylabel('feasibility')
title('Feasibility -  RS - Variable Weight')
xlim([1 4])

% load Feasibility_5WeightPairs
R_thresholds = [1,1.5:0.1:5.5,6];
figure
plot(R_thresholds, Feasibility_RS)
legend('RS 10:1','RS 3:1', 'RS 1:1', 'RS 1:3', 'RS 1:10')
xlabel('R_th (bits/s/Hz)')
ylabel('feasibility')
title('Feasibility -  RS - Variable Weight')

figure
plot(R_thresholds, Feasibility_MULP)
legend('MULP 10:1','MULP 3:1', 'MULP 1:1', 'MULP 1:3', 'MULP 1:10')
xlabel('R_th (bits/s/Hz)')
ylabel('feasibility')
title('Feasibility -  MULP - Variable Weight')

