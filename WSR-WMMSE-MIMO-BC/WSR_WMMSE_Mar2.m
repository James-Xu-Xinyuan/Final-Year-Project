clear all
P = 4; % transmit antennas
K = 2; % two users
Q =2; % each user has two antennas

SNRdB =[-10,0,10,20];

weight=[-3,-1:0.05:1,3];
u2=10.^weight;
u1=ones(1,length(u2));
rate1=zeros(length(u2),1);
rate2=zeros(length(u2),1);
bias  = 1;
% bias  = 1/sqrt(2);

theta = pi/3;

for s = 1:length(SNRdB)
    Etx = db2pow(SNRdB(s))
    for i=1:length(u2)
    R1_average = 0;
    R2_average = 0;
    
    iteration  = 1;
    for repetition=1:iteration
         H(:,:,1)=[1,1,1,1;1,1,1,1];
         H(:,:,2)=[1, exp(1i*1*theta),exp(1i*2*theta),exp(1i*3*theta);1, exp(1i*1*theta),exp(1i*2*theta),exp(1i*3*theta)]*bias;
%         H(:,:,1)=(1/sqrt(2))*(randn(Q,P)+1i*randn(Q,P));
%         H(:,:,2)=(1/sqrt(2))*(randn(Q,P)+1i*randn(Q,P));
        % one random pair of channel realization

       %initialize transmit filter for the 2 user
        B(:,:,1)=H(:,:,1)';
        B(:,:,2)=H(:,:,2)';

        % energy of current random transmit filter  
        Etf = trace(B(:,:,1)*B(:,:,1)')+trace(B(:,:,2)*B(:,:,2)');

        % imagining about this normalization
        B(:,:,1)= B(:,:,1)./sqrt( Etf ).*sqrt(Etx);
        B(:,:,2)=B(:,:,2)./sqrt( Etf ).*sqrt(Etx);
        Etf = trace(B(:,:,1)*B(:,:,1)')+trace(B(:,:,2)*B(:,:,2)');

        % this threshold probably slows things down
%         tolerence = 10^-6;
        tolerence = 10^-6;
        
%         R1 = 1;     R2 = 1;
        diff = 1;
        WSR = 0;
        counter = 0 ;

        while (diff > tolerence)
            
            % (6) effective noise covariance matrix
%             R = eye(Q,Q)+H(:,:,1)*B(:,:,2)*B(:,:,2)'*H(:,:,1)'+H(:,:,2)*B(:,:,1)*B(:,:,1)'*H(:,:,2)';
            R1 = eye(Q,Q)+H(:,:,1)*B(:,:,2)*B(:,:,2)'*H(:,:,1)';
            R2 = eye(Q,Q)+H(:,:,2)*B(:,:,1)*B(:,:,1)'*H(:,:,2)';


            % (7) MMSE receive filter
            A(:,:,1) = B(:,:,1)'*H(:,:,1)'*inv(H(:,:,1)*B(:,:,1)*B(:,:,1)'*H(:,:,1)'+R1);
            A(:,:,2) = B(:,:,2)'*H(:,:,2)'*inv(H(:,:,2)*B(:,:,2)*B(:,:,2)'*H(:,:,2)'+R2);

            % (8) MMSE-matrix
            E(:,:,1)= inv( eye(Q,Q)+B(:,:,1)'*H(:,:,1)'*inv(R1)*H(:,:,1)*B(:,:,1) );
            E(:,:,2)= inv( eye(Q,Q)+B(:,:,2)'*H(:,:,2)'*inv(R2)*H(:,:,2)*B(:,:,2) );

            % (21) weight matrix of user 
            W(:,:,1) = u1(i)*inv(E(:,:,1));
            W(:,:,2) = u2(i)*inv(E(:,:,2));

            % is there easier way to express these block-diagonal matrices?
            W_diag = [W(:,:,1),zeros(2,2);zeros(2,2),W(:,:,2)]; 
            A_diag = [A(:,:,1),zeros(2,2);zeros(2,2),A(:,:,2)]; 
            H_stacked = [H(:,:,1); H(:,:,2)];

            % (22)
%             Etf = trace(B(:,:,1)*B(:,:,1)')+trace(B(:,:,2)*B(:,:,2)');
            B_bar = H_stacked'*A_diag'*W_diag*A_diag*H_stacked;
            B_bar = B_bar+trace(W_diag*(A_diag*(A_diag')))/Etx*eye(size(B_bar));
            B_bar = inv(B_bar)*H_stacked'*A_diag'*W_diag;

            % (23)
            b = sqrt(Etx/trace(B_bar*B_bar'));
            B_WMMSE = b *B_bar;
            B(:,:,1)=B_WMMSE(:,1:2);
            B(:,:,2)=B_WMMSE(:,3:4);
            Etf = trace(B(:,:,1)*B(:,:,1)')+trace(B(:,:,2)*B(:,:,2)');

            R1 = abs(log2(det(inv(E(:,:,1)))));
            R2 = abs(log2(det(inv(E(:,:,2)))));
            
            WSR_new = u1(i)*R1+u2(i)*R2;

            diff = abs(WSR_new/WSR-1);
            WSR = WSR_new;
            
%             counter = counter +1; % sometimes the code do not converge
%             if counter >= 2000
%                 break
%             end
        end
%         if counter~=2000
%             R1_average = R1_average + R1/100;
%             R2_average = R2_average + R2/100;
%         end
        R1_average = R1_average + R1/iteration;
        R2_average = R2_average + R2/iteration;
    end
    [i, R1_average, R2_average]
    rate1(i) =  R1_average;
    rate2(i) =  R2_average;
    end
    
    figure
    % plot(rate1,rate2,'*')
    k = convhull(rate1, rate2);
    x1=rate1(k);
    y1=rate2(k);
    xx = floor(x1.*10^5)./(10^5);
    ind = find(xx==0);
    [l,ind_ini]= max(x1);
    plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'bo-')
    title("WSR-MMSE,"+SNRdB(s)+"dB")
    xlabel('R1')
    ylabel('R2')
end




