function [ V ] = x2V_pf( x, tnr )
%pf_x2V Converts the x vector to complex voltagess
%   Detailed explanation goes here
    [pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);

    V = tnr.V0;
    V(pq) = x(npv+npq+1:npv+2*npq).*exp(1j*x(npv+1:npv+npq));
    V(pv) = abs(V(pv)).*exp(1j*x(1:npv));

end

