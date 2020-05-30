function [WMSE,p1,p2]=MULP_update_P(M,weights,H_est,SNR,U1,U2,t1,t2,psi_1,psi_2,f1,f2,v1,v2);
% The algorithm used is from papers:
% Weighted Sum-Rate Maximization using Weighted MMSE for MIMO-BC 
% Beamforming Design, 2008 TWC, Christensen, Agarval, Carvalho, Cioffi
% "Sum-Rate Maximization for Linearly Precoded Downlink Multiuser MISO Systems with Partial CSIT: A Rate-Splitting Approach"
% by Hamdi Joudeh and Bruno Clerks

    u1 = weights(1);    u2 = weights(2);
    [Nr,Nt,K] = size(H_est);  
    
    cvx_begin quiet
    
    variable p1(Nt,Nr) complex
    variable p2(Nt,Nr) complex
 
    expression P_total;

    private_1 = quad_form(p1,psi_1)+quad_form(p2,psi_1)+t1-2*real(f1'*p1)+U1-v1;
    private_2 = quad_form(p1,psi_2)+quad_form(p2,psi_2)+t2-2*real(f2'*p2)+U2-v2;
    
    obj_func = u1*private_1+u2*private_2;
    minimize(obj_func)

    P_total = trace(p1'*p1)+trace(p2'*p2);

    subject to
        P_total <= SNR
        
    cvx_end

    WMSE=obj_func;
      
end