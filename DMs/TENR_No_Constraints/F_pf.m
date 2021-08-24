function [F, delta_injection_F] = F_pf(Vx,Vy,l,tnr,Bus_fstage,a_bus_p,b_bus_q)
%pfF Calculate the residual of power flow Jacobian
%   Detailed explanation goes here

    [pv, pq]  = deal(tnr.pv, tnr.pq);

    V = Vx + 1j*Vy;
    
    delta_injection = sparse(max(size(tnr.dS)),1);
    
    delta_injection(Bus_fstage) = a_bus_p.*real(tnr.dS(Bus_fstage,1)) + 1j*(b_bus_q.*(imag(tnr.dS(Bus_fstage))));
    
    delta_injection_F =  delta_injection;
        
%     residual_current = (tnr.Ybus*V);
    
    residual =  (V .* conj(tnr.Ybus * V)) -  ( tnr.dS + l.*(delta_injection) );
    
    residual_PV = (Vx.^2) + (Vy.^2)- ((abs(tnr.V0)).^2);
    
    F = [ 
      real(residual(pv));
      real(residual(pq));      
      imag(residual(pq));
      residual_PV(pv);
    ];
end

