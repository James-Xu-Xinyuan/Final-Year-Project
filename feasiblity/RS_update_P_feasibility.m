function [feasible,AWSMSE,p1,p2,pc]=RS_update_P_feasibility(Rth, H_est,weights,SNR,Uc1,Uc2,U1,U2,tc1,tc2,t1,t2,...
                         psi_c1,psi_c2,psi_1,psi_2,fc1,fc2,f1,f2,vc1,vc2,v1,v2)
    % introduce Rth, threshold rate / QoS requirement
    % introduce Boolean Variable 'feasible"
    feasible = 1; % default output
    u1=weights(1);
    u2=weights(2);
    h1 = H_est(:,:,1);
    h2 = H_est(:,:,2);
    [Nr,Nt,K] = size(H_est) ;
    try
        cvx_begin quiet

        % equations 29

        variable p1(Nt,Nr) complex
        variable p2(Nt,Nr) complex
        variable pc(Nt,Nr) complex

        expression P_total;
        expression R1;
        expression R2;

        % constrains on commmon part
        common_1=quad_form(pc,psi_c1)+quad_form(p1,psi_c1)+quad_form(p2,psi_c1)+tc1-2*real(fc1'*pc)+Uc1-vc1;
        common_2=quad_form(pc,psi_c2)+quad_form(p1,psi_c2)+quad_form(p2,psi_c2)+tc2-2*real(fc2'*pc)+Uc2-vc2;
        c=max(common_1,common_2);
        
        Rc1 = 1 - common_1;
        Rc2 = 1 - common_2;
        Rc = min(Rc1, Rc2);
        % Rc has problem defined above, not sure why
        Rc = 1-c;

        private_1=quad_form(p1,psi_1)+quad_form(p2,psi_1)+t1-2*real(f1'*p1)+U1-v1;
        private_2=quad_form(p1,psi_2)+quad_form(p2,psi_2)+t2-2*real(f2'*p2)+U2-v2;
        
        object_func=u1*(c+private_1)+u2*private_2;
        minimize(object_func)
       
        % constraints on total power
        P_total=trace(p1'*p1)+trace(p2'*p2)+trace(pc'*pc);
        R1 = 1 - private_1 + Rc;
        R2 = 1 - private_2;

        subject to
            P_total <= SNR;
            R1 >= Rth;
            R2 >= Rth;

        cvx_end
        
        AWSMSE=object_func;
        
        if cvx_status  == "Solved"
            feasible = 1;
        else
            feasible = 0;
        end
            
    catch ME
        ME.identifier;
        feasible = 0;   
        AWSMSE = 0;
        p1 = zeros(4,1); % return zeros
        p2 = zeros(4,1); % position occupier only
        pc = zeros(4,1); % values does not matter
    end % try catch 
end %function