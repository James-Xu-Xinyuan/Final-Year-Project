% This code is to reproduced Figure 8 of 
% "Sum-Rate Maximization for Linearly Precoded Downlink Multiuser MISO Systems with Partial CSIT: A Rate-Splitting Approach"
% by Hamdi Joudeh and Bruno Clerks

% Xinyuan Xu
% April 2020

function [Feasible] = RS_Rate_Region(Rth,M,weights,H_est,H_error_1,H_error_2,SNR,tolerance)  
% function [Feasible] = RS_Rate_Region(M,weights,H_est,H_error_1,H_error_2,SNR,tolerance)  
[Nr,Nt,K] = size(H_est);
% introduce Rth, threshold rate / QoS requirement
% introduce Boolean Variable 'feasible'

% % test only
% Rth = 0; % 0, should be all feasible, 100 low SNR should not be feasible
% Rth = 0,0.5, 1, 10, 100; 
% in an experiement, 1.3 is a kind of 50% 50% for SNR around 20dB, theta =
% pi/3, weight = [1,1]
% weights = [1,1];
% SNR = 20;
% test variables -> to fix

% Feasible: the outcome of trying all 4 initialization methods
feasible = zeros(1,4);

for i=1:M
    H_1(:,:,i) = H_est(:,:,1)+H_error_1(:,:,i);
    H_2(:,:,i) = H_est(:,:,2)+H_error_2(:,:,i);
end

for method = 1:4 % loop through 4 different initializaiton methods
    [p1,p2,pc] = initialize_precoders(H_est, SNR, method);

%     loop=1;
%     AWSMSE_old=0.1;
%     count=0;
%     while (loop)
        % reset realizations
        tc1 = 0;        tc2 = 0;
        t1 = 0;         t2 = 0;

        Uc1 = 0;        Uc2 = 0;
        U1 = 0;         U2 = 0;

        psi_c1 = 0;     psi_c2 = 0;
        psi_1 = 0;      psi_2 = 0;

        fc1 = zeros(Nt,Nr);     fc2 = zeros(Nt,Nr);
        f1 = zeros(Nt,Nr);      f2 = zeros(Nt,Nr);

        vc1 = 0;        vc2 = 0;
        v1 = 0;         v2 = 0;

        for i=1:M
            h1 = H_1(:,:,i);        h2 = H_2(:,:,i);

            % average receive power for a given channel state. equation (3)
            I1 = abs(h1*p2)^2+1;        I2 = abs(h2*p1)^2+1;
            T1 = abs(h1*p1)^2+I1;       T2 = abs(h2*p1)^2+I2;
            Tc1 = abs(h1*pc)^2+T1;      Tc2 = abs(h2*pc)^2+T2;

            % Optimum Minimum MSE equalizers. equation (20)
            gc1 = pc'*h1'*inv(Tc1);     gc2 = pc'*h2'*inv(Tc2);
            g1 = p1'*h1'*inv(T1);       g2 = p2'*h2'*inv(T2);

            % MMSE, equation (21)   % Ick = Tk
            MMSE_c1 = inv(Tc1)*T1;      MMSE_c2 = inv(Tc2)*T2;
            MMSE1 = inv(T1)*I1;         MMSE2 = inv(T2)*I2;

            % optimum MMSE weights, between equation 24 and 25
            uc1 = inv(MMSE_c1);         uc2 = inv(MMSE_c2);
            u1 = inv(MMSE1);            u2 = inv(MMSE2);

            Uc1 = Uc1+uc1;              Uc2 = Uc2+uc2;
            U1 = U1+u1;                 U2 = U2+u2;

            % D, 1) Updating the Equalizers and Weights
            % sum of realizations
            tc1 = tc1+uc1*abs(gc1)^2;   tc2 = tc2+uc2*abs(gc2)^2;
            t1 = t1+u1*abs(g1)^2;       t2 = t2+u2*abs(g2)^2;

            psi_c1 = psi_c1+(uc1*abs(gc1)^2)*h1'*h1;
            psi_c2 = psi_c2+(uc2*abs(gc2)^2)*h2'*h2;
            psi_1 = psi_1+(u1*abs(g1)^2)*h1'*h1;
            psi_2 = psi_2+(u2*abs(g2)^2)*h2'*h2;

            fc1 = fc1+uc1*h1'*gc1';         fc2 = fc2+uc2*h2'*gc2';
            f1 = f1+u1*h1'*g1';             f2 = f2+u2*h2'*g2';

            vc1 = vc1+log2(uc1);            vc2 = vc2+log2(uc2);
            v1 = v1+log2(u1);               v2 = v2+log2(u2);
        end % for

        % averaging over the correpsonding realizaitons
        Uc1 = Uc1/M;        Uc2 = Uc2/M;
        U1 = U1/M;          U2 = U2/M;

        tc1 = tc1/M;        tc2 = tc2/M;
        t1 = t1/M;          t2 = T2/M;

        psi_c1 = psi_c1./M;     psi_c2 = psi_c2./M;
        psi_1 = psi_1./M;       psi_2 = psi_2./M;

        fc1 = fc1./M;       fc2 = fc2./M;
        f1 = f1./M;         f2 = f2./M;

        vc1 = vc1/M;        vc2 = vc2/M;
        v1 = v1/M;          v2 = v2/M;

        %Step 4: Update P
        [feasible(method), AWSMSE,p1,p2,pc]=RS_update_P(Rth, H_est,weights,SNR,Uc1,Uc2,U1,U2,tc1,tc2,t1,t2,...
                             psi_c1,psi_c2,psi_1,psi_2,fc1,fc2,f1,f2,vc1,vc2,v1,v2);
        if feasible(method) == 0 
%             loop = 0;
%             Capacity(1) = 0 ;   Capacity(2) = 0 ;
%             P_c = 0;
%             break; % while
        end

%         if abs((AWSMSE-AWSMSE_old)/AWSMSE_old)<=tolerance
%             loop=0;
%         else
%             AWSMSE_old = AWSMSE;
%             count = count+1;
%         end
% 
%         if count>=500
%             loop=0;
%             break;
%         end

%     end %while

    
end % for method

if feasible == zeros(1,4)
    Feasible = 0;
else
    Feasible = 1;
end



% P_c = real(trace(pc*pc'));

% this experiment concerns feasiblity only, so capacity calculation is not
% needed
% for i=1:M
%     h1 = H_1(:,:,i);            h2 = H_2(:,:,i);
%     
%     % equation 3
%     I1 = abs(h1*p2)^2+1;        I2 = abs(h2*p1)^2+1;;
%     T1 = abs(h1*p1)^2+I1;       T2 = abs(h2*p2)^2+I2;
% %     Tc1 = abs(h1*pc)^2+T1;      Tc2 = abs(h2*pc)^2+T2;
%     % Ick = Tk      %SINR, gama on paper        % equation 4
%     yc1 = abs(h1*pc)^2*inv(T1);        yc2 = abs(h2*pc)^2*inv(T2);    
%     y1 = abs(h1*p1)^2*inv(I1);         y2 = abs(h2*p2)^2*inv(I2);     
%     
%     % equation 5
%     Rc1 = real(log2(1+yc1));           
%     Rc2 = real(log2(1+yc2));
%     Rc = min(Rc1,Rc2);
%     R1(i) = real(log2(1+y1))+Rc; %  R1(:,:,i)=
%     R2(i) = real(log2(1+y2));
% end
% 
% Capacity(1)=sum(R1)/M;
% Capacity(2)=sum(R2)/M;




          
          
          
 