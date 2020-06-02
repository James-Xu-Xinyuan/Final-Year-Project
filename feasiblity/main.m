% Xinyuan Xu
% May 2020
% feasibility of RS with partial CSIT with extra QoS contraints
% main, the basic version:
% fixed SNR, fixed weight pair
% the only variable - Rth in array R_thresholds

clear all; 

tolerance = 10^-2;
Nt=4; 
Nr=1; 
SNRdB=[10 20 30]; 
% SNRdB = 30;
weights = [1,1];
M = 10;% Realization of Sample Average Functions
% M = 1000;
    
% R_thresholds = [0.01, 0.1,0.5, 1:0.05:2,2.5,3];
% R_thresholds = [1:1:12]; 
% test input, to find roughly the region of interst
R_thresholds = [1,1.5,2:0.1:6,6.5,7]; 
Feasibility_RS = zeros(length(SNRdB),length(R_thresholds));
Feasibility_MULP = zeros(length(SNRdB),length(R_thresholds));
NoTest = 10; % local
% NoTest = 100; HPC

for i_SNR = 1:3    
    SNR = 10^(SNRdB(i_SNR)/10);
    P_t=SNR; %total transmission power, unit norm variance
    P_e=P_t^(-0.6); %error variance

for i_th = 1:length(R_thresholds) % loop throught different threshold values
%     i_th
    Feasibility_rs = 0; % counter/register
    Feasibility_mulp = 0;
    Rth = R_thresholds(i_th);
    for i_test = 1:NoTest
        % test feasiblity with a lot of random channels
        H_est(:,:,1)= (randn(1,Nt)+j*randn(1,Nt));
        H_est(:,:,1) = H_est(:,:,1)/norm(H_est(:,:,1))*sqrt(Nt);
        H_est(:,:,2)= (randn(1,Nt)+j*randn(1,Nt));
        H_est(:,:,2) = H_est(:,:,2)/norm(H_est(:,:,2))*sqrt(Nt);

        for i=1:M
            H_error_1(:,:,i)=((randn(Nr,Nt)+j*randn(Nr,Nt))/sqrt(2))*sqrt(P_e);
            H_error_2(:,:,i)=((randn(Nr,Nt)+j*randn(Nr,Nt))/sqrt(2))*sqrt(P_e);
        end

        [Feasibile_order1_RS] = RS_Rate_Region_feasibility(Rth,M,weights,H_est,H_error_1,H_error_2,SNR,tolerance);
        % change the order
        H_est1(:,:,1)=H_est(:,:,2);
        H_est1(:,:,2)=H_est(:,:,1);
        weights1(1)=weights(2);
        weights1(2)=weights(1);
        % different decoding order
        [Feasibile_order2_RS] = RS_Rate_Region_feasibility(Rth,M,weights1,H_est1,H_error_2,H_error_1,SNR,tolerance);

        if Feasibile_order1_RS == 0 && Feasibile_order2_RS == 0
%             Feasibility_rs = Feasibility_rs;
        else
            Feasibility_rs = Feasibility_rs + 1;
        end

        Feasible_MULP = MULP_Rate_Region_feasibility(Rth, M,weights,H_est,H_error_1,H_error_2,SNR,tolerance);
        if Feasible_MULP == 0 
%             Feasibility_mulp = Feasibility_mulp;
        else
            Feasibility_mulp = Feasibility_mulp + 1;
        end

    end % for i_test
    Feasibility_rs = Feasibility_rs / NoTest; 
    Feasibility_RS(iSNR,i_th) = Feasibility_rs;
    % approximate real percentage feasibility (probability) as percentage of feasbile channels in tests 
    Feasibility_mulp = Feasibility_mulp/ NoTest; 
    Feasibility_MULP(iSNR,i_th) = Feasibility_mulp;
end     % for i_th
end     % for i_SNR
save('Feasibilities_30dB', 'R_thresholds','Feasibility_RS','Feasibility_MULP');

% plot(R_thresholds, Feasibility_RS)
% title('Feasibility, RS, 10dB, equal weights, 1000random test channels')
% hold on
% plot(R_thresholds, Feasibility_MULP)
% title('Feasibility, NoRS, 10dB, equal weights, 1000random test channels')
