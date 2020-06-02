function [p1,p2,pc] = initialize_precoders(H_est, SNR, m)
% m: method number
warning('off')

h1 = H_est(:,:,1);
h2 = H_est(:,:,2);
[Nr,Nt,K] = size(H_est);
P_c = SNR*0.6;
P_k = (SNR-P_c)/2;

% private parts
% MRC 
p1_mrc = h1'/norm(h1)*sqrt(P_k);
p2_mrc = h2'/norm(h2)*sqrt(P_k);
% ZF % Es =  SNR
H = [h1; h2];
Gzf = sqrt(Nt/SNR)*inv(H'*H)*H';
p1_zf = Gzf(:,1)/norm(Gzf(:,1))*sqrt(P_k);
p2_zf = Gzf(:,2)/norm(Gzf(:,2))*sqrt(P_k);

% common part
% SVD
H = [h1' h2'];
[U,~,~] = svd(H);
p_c_hat = U(:,1);
pc_svd = p_c_hat*sqrt(P_c);   
pc_rnd = randn(Nt,1)+j*randn(Nt,1);
pc_rnd = pc_rnd/norm(pc_rnd)*sqrt(P_c);

if m == 1
    p1 = p1_mrc;
    p2 = p2_mrc;
    pc = pc_svd;
end

if m == 2
    p1 = p1_mrc;
    p2 = p2_mrc;
    pc = pc_rnd;
end

if m == 3
    p1 = p1_zf;
    p2 = p2_zf;
    pc = pc_svd;
end

if m == 4
    p1 = p1_zf;
    p2 = p2_zf;
    pc = pc_rnd;
end


end