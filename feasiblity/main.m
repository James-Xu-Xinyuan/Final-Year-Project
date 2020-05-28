% Xinyuan Xu
% May 2020

clear all; 

tolerance = 10^-3;
Nt=4; 
Nr=1; 
% SNRdB=[10 20 30]; 
SNRdB = 10;
weights = [1,1];
theta=pi/3;
M=100;% Realization of Sample Average Functions

% weight=[-3 -1:0.05:1 3];
% u2=10.^weight;
% u1=ones(1,length(u2));
% H_est(:,:,1)=[1 1 1 1];
% H_est(:,:,2)=[1 exp(1i*theta) exp(1i*2*theta) exp(1i*3*theta)]*bias;

% for i_SNRdB=1:length(SNRdB)
    
    %Channel error
%     SNR = 10^(SNRdB(i_SNRdB)/10);
SNR = 10^(SNRdB/10);
P_t=SNR; %total transmission power, unit norm variance
P_e=P_t^(-0.6); %error variance
    
R_thresholds = [0.01, 0.1,0.5, 1:0.05:2,2.5,3];
% R_thresholds = [1, 1.3, 2];
% in a test: 10dB, equal weight, 1->100%, 1.3->78%
Feasibility_RS = zeros(1,length(R_thresholds));
NoTest = 100;
    
for i_th = 1:length(R_thresholds); % loop throught different threshold values
    Feasibility_rs = 0; % counter/register
    Rth = R_thresholds(i_th);
    for i_test = 1:NoTest
            % test feasiblity with 100 random channels
            H_est(:,:,1)= (randn(1,Nt)+j*randn(1,Nt));
            H_est(:,:,1) = H_est(:,:,1)/norm(H_est(:,:,1))*sqrt(Nt);
            H_est(:,:,2)= (randn(1,Nt)+j*randn(1,Nt));
            H_est(:,:,2) = H_est(:,:,2)/norm(H_est(:,:,2))*sqrt(Nt);
            
            for i=1:M
                H_error_1(:,:,i)=((randn(Nr,Nt)+j*randn(Nr,Nt))/sqrt(2))*sqrt(P_e);
                H_error_2(:,:,i)=((randn(Nr,Nt)+j*randn(Nr,Nt))/sqrt(2))*sqrt(P_e);
            end

%     for i_weight=1:length(u1)
%         SNRdB(i_SNRdB)
%         i_weight
%         weights=[u1(i_weight),u2(i_weight)];

%         [Capacity_RS_order1, Capacity_RS_order2, P_common1, P_common2] = RS_rateRegion(M,weights,H_BC_estimate,H_BC_error_1,H_BC_error_2,SNRdB(i_SNRdB),tolerance)
            [Feasibile_order1_RS] = RS_Rate_Region(Rth,M,weights,H_est,H_error_1,H_error_2,SNR,tolerance);
            % change the order
            H_est1(:,:,1)=H_est(:,:,2);
            H_est1(:,:,2)=H_est(:,:,1);
            weights1(1)=weights(2);
            weights1(2)=weights(1);
            % different decoding order
            [Feasibile_order2_RS] = RS_Rate_Region(Rth,M,weights1,H_est1,H_error_2,H_error_1,SNR,tolerance);

            if Feasibile_order1_RS == 0 && Feasibile_order2_RS == 0
                Feasibility_rs = Feasibility_rs;
            else
                Feasibility_rs = Feasibility_rs + 1;
            end
        
        end % for i_test
        Feasibility_rs = Feasibility_rs / NoTest; 
        Feasibility_RS(i_th) = Feasibility_rs;
        % approximate real percentage feasibility (probability) as percentage of feasbile channels in tests 
end     % for i_th
%         C_RS_order2(2)=C_order2(1);
%         C_RS_order2(1)=C_order2(2);
%         
%         C_RS_order1_user1(i_SNRdB,i_weight) = C_RS_order1(1);
%         C_RS_order1_user2(i_SNRdB,i_weight) = C_RS_order1(2);
%         C_RS_order2_user1(i_SNRdB,i_weight) = C_RS_order2(1);
%         C_RS_order2_user2(i_SNRdB,i_weight) = C_RS_order2(2);
% 
%         C_NoRS = MULP_rateRegion(M,weights,H_est,H_error_1,H_error_2,SNRdB(i_SNRdB),tolerance)
%         C_NoRS_user1(i_SNRdB,i_weight)=C_NoRS(1);
%         C_NoRS_user2(i_SNRdB,i_weight)=C_NoRS(2);

%     end %end looping user weights
    
% end %end looping SNR
 
% x=[C_RS_order1_user1(1,:) C_RS_order2_user1(1,:)];
% y=[C_RS_order1_user2(1,:) C_RS_order2_user2(1,:)];
% k=convhull(x,y);
% x1 = x(k);
% y1 = y(k);
% xx=floor(x1.*10^(5))./(10^(5));
% ind=find(xx==0);
% [~,ind_ini]=max(xx);
% plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
% hold on,grid on

% x=C_NoRs_user1(1,:);
% y=C_NoRs_user2(1,:);
% k=convhull(x,y);
% x1 = x(k);
% y1 = y(k);
% xx=floor(x1.*10^(5))./(10^(5));
% ind=find(xx==0);
% [~,ind_ini]=max(x1);
% plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
% hold on,grid on

% x=[C_RS_order1_user1(2,:) C_RS_order2_user1(2,:)];
% y=[C_RS_order1_user2(2,:) C_RS_order2_user2(2,:)];
% k=convhull(x,y);
% x1 = x(k);
% y1 = y(k);
% xx=floor(x1.*10^(5))./(10^(5));
% ind=find(xx==0);
% [~,ind_ini]=max(xx);
% plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
% hold on,grid on

% x=C_NoRs_user1(2,:);
% y=C_NoRs_user2(2,:);
% k=convhull(x,y);
% x1 = x(k);
% y1 = y(k);
% xx=floor(x1.*10^(5))./(10^(5));
% ind=find(xx==0);
% [~,ind_ini]=max(x1);
% plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
% hold on,grid on

% x=[C_RS_order1_user1(3,:) C_RS_order2_user1(3,:)];
% y=[C_RS_order1_user2(3,:) C_RS_order2_user2(3,:)];
% k=convhull(x,y);
% x1 = x(k);
% y1 = y(k);
% xx=floor(x1.*10^(5))./(10^(5));
% ind=find(xx==0);
% [~,ind_ini]=max(xx);
% plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
% hold on,grid on

% x=C_NoRs_user1(3,:);
% y=C_NoRs_user2(3,:);
% k=convhull(x,y);
% x1 = x(k);
% y1 = y(k);
% xx=floor(x1.*10^(5))./(10^(5));
% ind=find(xx==0);
% [~,ind_ini]=max(x1);
% plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
% hold on,grid on

% title('RS with partial CSIT')
% % legend('10dB RS','10dB MU-LP','20dB RS','20dB MU-LP','30dB RS','30dB MU-LP')
% xlabel('R1 (bit/s/Hz)')
% ylabel('R2 (bit/s/Hz)')

% save('Partial_CSIT_Nt4','C_NoRS_user1','C_NoRS_user2','C_RS_order1_user1','C_RS_order2_user1','C_RS_order1_user2','C_RS_order2_user2')
save('Feasibilities', 'R_thresholds','Feasibility_RS');

