% This code is to reproduced Figure 8 of 
% "Sum-Rate Maximization for Linearly Precoded Downlink Multiuser MISO Systems with Partial CSIT: A Rate-Splitting Approach"
% by Hamdi Joudeh and Bruno Clerks

% Xinyuan Xu
% 1st written: April 2020
% revised: 2020 May 29

function [AWSMSE,p1,p2,pc ]=RS_update_P(H_est,weights,SNR,Uc1,Uc2,U1,U2,tc1,tc2,t1,t2,...
                         psi_c1,psi_c2,psi_1,psi_2,fc1,fc2,f1,f2,vc1,vc2,v1,v2)
   
    u1=weights(1);
    u2=weights(2);
    
    [Nr,Nt,K] = size(H_est) ;
    
    cvx_begin quiet
    
    % equations 29
    
    variable p1(Nt,Nr) complex
    variable p2(Nt,Nr) complex
    variable pc(Nt,Nr) complex
 
    expression constraints;

    % constrains on commmon part
    common_1=quad_form(pc,psi_c1)+quad_form(p1,psi_c1)+quad_form(p2,psi_c1)+tc1-2*real(fc1'*pc)+Uc1-vc1;
    common_2=quad_form(pc,psi_c2)+quad_form(p1,psi_c2)+quad_form(p2,psi_c2)+tc2-2*real(fc2'*pc)+Uc2-vc2;
    c=max(common_1,common_2);
  
    private_1=quad_form(p1,psi_1)+quad_form(p2,psi_1)+t1-2*real(f1'*p1)+U1-v1;
    private_2=quad_form(p1,psi_2)+quad_form(p2,psi_2)+t2-2*real(f2'*p2)+U2-v2;
        
    object_func=u1*(c+private_1)+u2*private_2;
    minimize(object_func)
    
    % constraints on total power
    constraints=trace(p1'*p1)+trace(p2'*p2)+trace(pc'*pc);

    subject to
        constraints<=SNR
           
    cvx_end

    AWSMSE=object_func;
    
end