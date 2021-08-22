function [ x ] = V2x_pf( V, tnr )
%pf_V2x Convert voltage vector to x variable
%   Detailed explanation goes here
    [pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);

    Vm = abs(V);
    Va = angle(V);
    
    x = zeros(2*npq + npv,1);
    x(1:npv+npq) = Va([pv;pq]);
    x(npv+npq+1:end) = Vm(pq);

end

