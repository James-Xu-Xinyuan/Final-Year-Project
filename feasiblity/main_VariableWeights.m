clear all; 
cvx_setup;
LASTN=maxNumCompThreads(1);
array_index=getenv('PBS_ARRAY_INDEX');
i_th=str2num(array_index);
i1 = str2num(array_index);
rng(i1); 

% Xinyuan Xu
% June 2020
% feasibility of RS with partial CSIT with extra QoS contraints
% main, the basic version:
% fixed SNR, variable weight pairs
% key variable - Rth in array R_thresholds

tolerance = 10^-3;
Nt=4; 
Nr=1; 
SNRdB = 10; % rate region should intersect with axis at around 3
SNR = 10^(SNRdB/10);
P_t=SNR; %total transmission power, unit norm variance
P_e=P_t^(-0.6); %error variance

u2 = [10,3,1,1/3,0.1];
u1 = ones(1,length(u2));
M = 1000;
R_thresholds = [1,1.5:0.05:3]; 
NoTest = 100; %HPC

for i_weight = 1:5    
    weights = [u1(i_weight), u1(i_weight)];

% for i_th = 1:length(R_thresholds) % loop throught different threshold values
%     i_threshold
    Feasibility_rs = 0; % counter/register
    Feasibility_mulp = 0;
%     i_th
    Rth = R_thresholds(i_th);
    for i_test = 1:NoTest
% use batch processing in HPC to replace this loop
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

        Feasible_MULP = MULP_Rate_Region_feasibility(Rth, M,weights,H_est,H_error_1,H_error_2,SNR,tolerance);

        % caution! save data operation! read write!
%         load Feasibility_Jun1.mat
        if Feasibile_order1_RS == 0 && Feasibile_order2_RS == 0
%             Feasibility_rs = Feasibility_rs;
        else
            Feasibility_rs = Feasibility_rs  + 1;
        end

        if Feasible_MULP == 0 
%             Feasibility_mulp = Feasibility_mulp;
        else
            Feasibility_mulp = Feasibility_mulp + 1;
        end        


    end % for i_test
    
    % caution! save data operation! read write!
    % File might be corrupt.
    load Feasibility_Jun2_Weights.mat
    Feasibility_rs = Feasibility_rs / NoTest; 
    Feasibility_RS(i_weight,i_th) = Feasibility_rs;
    % approximate real percentage feasibility (probability) as percentage of feasbile channels in tests 
    Feasibility_mulp = Feasibility_mulp/ NoTest; 
    Feasibility_MULP(i_weight,i_th) = Feasibility_mulp;
    save('Feasibility_Jun2_Weights.mat','Feasibility_RS','Feasibility_MULP');
    clear Feasibility_MULP
    clear Feasibility_RS
% end     % for i_th
end     % for i_SNR