function [F,delta_ss] = F_pf( Vx,Vy,Slack_H,Slack_L,V_max,V_min,l,tnr,Bus_fstage,a_bus_p,b_bus_q)
%pfF Calculate the residual of power flow Jacobian
%   Detailed explanation goes here

    [pv, pq]  = deal(tnr.pv, tnr.pq);

    V = Vx + 1j*Vy;
    
    delta_ss = zeros(max(size(V)),1);
    delta_ss(Bus_fstage) = a_bus_p*(real(tnr.dS(Bus_fstage))) + 1j*(b_bus_q*(imag(tnr.dS(Bus_fstage))));
    
    
    delta_sss =  tnr.dS + (l.*delta_ss);
    
    residual =  (V .* conj(tnr.Ybus * V)) - delta_sss;
    
    residual_PV = (Vx.^2) + (Vy.^2)- ((abs(tnr.V0)).^2);
    
    residual_max = (Vx.^2) + (Vy.^2)- (V_max).^2 + (Slack_H).^2;
    
    residual_min = (Vx.^2) + (Vy.^2)- (V_min).^2 - (Slack_L).^2;
    
    F = [ 
      real(residual(pv));
      real(residual(pq));      
      imag(residual(pq));
      residual_PV(pv);
      residual_max(pq);
      residual_min(pq);
    ];
end

