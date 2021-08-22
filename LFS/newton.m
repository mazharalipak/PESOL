clc;
clearvars;

% function [Iter,Vx, Vy, Tole] = newton(data,step) 

step = 1;
maxIters = 20;                                                             % Max Iteration bound.................
data = case14;                                                         % Load case file ............................

%% loading case file .........................

mpc = ext2int( loadcase( data ) );

%% Initializing ........................................................

tnr = tnr_init(mpc);                                                       % Initial data format .......
[ pv, pq, npv, npq ]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);            % PV, PQ buses and ..........
nbus = 1 + npv + npq;                                                      % total number of buses.........

%% intial germ ..............................

Vx = real(tnr.V0);                                                         % real component of complex voltage phasor..............
Vy = imag(tnr.V0);                                                         % imag component of complex voltage phasor..............

% Vx = ones(size(tnr.V0));                                                 % real component of complex voltage phasor .....................
% Vy = 0.0*( ones(size(tnr.V0)) );                                         % imag component of complex voltage phasor..............

l = 1.0;                                                                   % Loadability parameter.................................

[ AAX, BBY ] = derivativexy(nbus,tnr);                                     % partial J/ partial x_{i} .........................

%% Newtpn Iterations .............................................

Tole = 1;  
Iter = 1;
counter = 0; 

while (Tole > 1e-4)  

F_B = tnr.F( Vx,Vy,l,tnr );                                                % Power blance mismatch vector (Rectangular formulation)........
J_F = tnr.J(Vx,Vy,tnr);                                                    % Power flow Jacobian .......
dv = - J_F \ F_B;                                                          % Newton correction step ....

%% optimal multiplier .................................

[ G_0, G_1, G_2, G_3 ] = optimal_multi(dv, J_F, AAX, BBY, F_B, nbus );
proot  = [ G_3  G_2  G_1  G_0 ];
r_x    = roots(proot);
[aff, bff] = min(abs( imag(r_x) ) );
step   = abs( real( r_x(bff) ) );

%% Updating state variables .....................

Vx([pv;pq])= Vx([pv;pq]) + step*dv(1:npv+npq);
Vy([pv;pq])= Vy([pv;pq]) + step*dv(npv+npq+1:end);


step_x (Iter) = step; 

%% Checking Tolerance ..................

Iter  = Iter + 1;
Tole = max( abs( F_B ) );
counter = counter + 1;
 
 if counter == maxIters
     break;
 end
 
end

%% check for non-convergence

if Iter > maxIters 
    disp('Newton did not converge, change step size'); 
 else
     disp ('Newton converged within 12 iterations')
end


%% complex nodal voltages .............................

 VL = Vx + 1i*Vy;

Tole

Iter

%end
