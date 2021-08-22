function [F] = F_pf( Vx,Vy,l,tnr )
%pfF Calculate the residual of power flow Jacobian
%   Detailed explanation goes here

    [pv, pq]  = deal(tnr.pv, tnr.pq);

    V = Vx + 1j*Vy;
    
    residual =  (V .* conj(tnr.Ybus * V))- (l.*tnr.dS);
    
    residual_PV = (Vx.^2) + (Vy.^2)- ((abs(tnr.V0)).^2);
    
    F = [ 
      real(residual(pv));
      real(residual(pq));      
      imag(residual(pq));
      residual_PV(pv);
    ];
end

