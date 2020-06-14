% This code is to reproduced Figure 8 of 
% "Sum-Rate Maximization for Linearly Precoded Downlink Multiuser MISO Systems with Partial CSIT: A Rate-Splitting Approach"
% by Hamdi Joudeh and Bruno Clerks

% Xinyuan Xu
% 1st written: April 2020
% revised: 2020 May 29 

% code templete for main
% different job variation in main{x} for different submission to HPC.

clear all; 
tolerance = 10^-4;
Nt = 4; 
Nr = 1; 
SNRdB = [10 20 30]; 
theta = pi/3;
M = 1000;
% Realization of Sample Average Functions
bias = 1;
weight = [-3 -1:0.05:1 3];
u2 = 10.^weight;
u1 = ones(1,length(u2));

H_est(:,:,1) = [1 1 1 1];
H_est(:,:,2) = [1 exp(1i*theta) exp(1i*2*theta) exp(1i*3*theta)]*bias;

C_NoRs_user1  = zeros(length(SNRdB), length(weight));
C_NoRs_user2  = zeros(length(SNRdB), length(weight));

for i_SNRdB=1:length(SNRdB)
    
    %Channel error
    SNR = 10^(SNRdB(i_SNRdB)/10);
    P_t=SNR; %total transmission power, unit norm variance
    P_e=P_t^(-0.6); %error variance
    
    for i=1:M
        H_error_1(:,:,i)=((randn(Nr,Nt)+j*randn(Nr,Nt))/sqrt(2))*sqrt(P_e);
        H_error_2(:,:,i)=((randn(Nr,Nt)+j*randn(Nr,Nt))/sqrt(2))*sqrt(P_e)*bias;
    end

    for i_weight=1:length(u1)
%         SNRdB(i_SNRdB)
%         i_weight
        weights=[u1(i_weight),u2(i_weight)];

        [C_RS_order1,Pc1] = RS_Rate_Region(M,weights,H_est,H_error_1,H_error_2,SNR,tolerance)
        % change the order
        H_est1(:,:,1)=H_est(:,:,2);
        H_est1(:,:,2)=H_est(:,:,1);
        weights1(1)=weights(2);
        weights1(2)=weights(1);
        % different decoding order
        [C_order2,Pc2] = RS_Rate_Region(M,weights1,H_est1,H_error_2,H_error_1,SNR,tolerance)
        C_RS_order2(2)=C_order2(1);
        C_RS_order2(1)=C_order2(2);
        
        C_RS_order1_user1(i_SNRdB,i_weight) = C_RS_order1(1);
        C_RS_order1_user2(i_SNRdB,i_weight) = C_RS_order1(2);
        C_RS_order2_user1(i_SNRdB,i_weight) = C_RS_order2(1);
        C_RS_order2_user2(i_SNRdB,i_weight) = C_RS_order2(2);

        C_NoRS = MULP_Rate_Region(M,weights,H_est,H_error_1,H_error_2,SNRdB(i_SNRdB),tolerance)
%         C_NoRS = MULP_rateRegion(M,weights,H_est,H_error_1,H_error_2,SNRdB(i_SNRdB),tolerance)
        C_NoRs_user1(i_SNRdB,i_weight)=C_NoRS(1);
        C_NoRs_user2(i_SNRdB,i_weight)=C_NoRS(2);

    end %end looping user weights
    
end %end looping SNR
 
x=[C_RS_order1_user1(1,:) C_RS_order2_user1(1,:)];
y=[C_RS_order1_user2(1,:) C_RS_order2_user2(1,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
hold on,grid on

x=C_NoRs_user1(1,:);
y=C_NoRs_user2(1,:);
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
hold on,grid on

x=[C_RS_order1_user1(2,:) C_RS_order2_user1(2,:)];
y=[C_RS_order1_user2(2,:) C_RS_order2_user2(2,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
hold on,grid on

x=C_NoRs_user1(2,:);
y=C_NoRs_user2(2,:);
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
hold on,grid on

x=[C_RS_order1_user1(3,:) C_RS_order2_user1(3,:)];
y=[C_RS_order1_user2(3,:) C_RS_order2_user2(3,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
hold on,grid on

x=C_NoRs_user1(3,:);
y=C_NoRs_user2(3,:);
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)))
hold on,grid on

title('RS with partial CSIT')
legend('10dB RS','10dB MU-LP','20dB RS','20dB MU-LP','30dB RS','30dB MU-LP')
xlabel('R1 (bit/s/Hz)')
ylabel('R2 (bit/s/Hz)')

% save('Partial_CSIT_Nt4','C_NoRS_user1','C_NoRS_user2','C_RS_order1_user1','C_RS_order2_user1','C_RS_order1_user2','C_RS_order2_user2')
legend('RS','MULP')

