%% test_run script

clc;
clearvars;

maxIters  = 60;
%% Initialization Block ...........................

testcase   = case14;

%far_east_indexing;                                           

%% Initial seed / germ for starting TENR Algorithm.........................

mpc = ext2int( ( loadcase( testcase ) ) );

% mpc.bus (:,3) = 1.0*mpc.bus(:,3);

%% Newton step size alpha,   0 < alpha <= 1,  normally a small value is selected......

step    = 0.9;                                                 

%% Scenario for test case (Please "Bus_fstage" values shouldn't exceed from the number of total buses in the network.)

Bus_fstage = [ 1,2,3,4,5,6,7,8,9,10,11,12,13,14]; %,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 ];                                                  % Buses, where loadability \lambda is subjected .................     
a_bus_p    = 1.0;                                                          % On/off status of active power injection "1" ON and "0" is OFF ................
b_bus_q    = 0.0;                                                          % On/off status of reactive power injection "1" ON and "0" is OFF ................

%% running TENR algorithm ..........................

%s_tole = 0;                                                               
%count_l = 0;

%while 1

% random selection of optimal step size \alpha for global convergence ......................

%step = randsample(n_stepx, 1);

% Lambda is maximum Loadability, Tole is tolerance and Iter is the total iterations ........ 

tic;
[lambda, Tole, Iter, Vx, Vy, tnr, Stabi,ur,vr, step_iter, Jxx,delta_S]=tnr_test(mpc, step, Bus_fstage, a_bus_p, b_bus_q, maxIters);   
toc;


%count_l = count_l+1;

%s_tole = Tole;

%% closing while loop ..............................

% if s_tole > 1e-4
%     
%     continue
% else 
%     break
% end
% 
% end

%% Observing right eigenvector corresponding to zero eig value..........

%  lam  = vr( 1:(length(vr)/2) ) +  1i* vr( (length(vr)/2)+1 :end);
%  [ a_value , b_indx ] = min( imag( lam(tnr.pq) ) );
%   
%  bus_indx_pq = tnr.pq(b_indx)

 lambda
 Tole

 
% emergency control actions to find the load which we need to shed to have ,....
% largest impact on value of lambda (increasing margin to solvability
% boundary)
 
Lamda_F = diag( -delta_S);
HJ = ur'*(Lamda_F);
AL = HJ/norm(HJ);
BL = ( AL(tnr.npv + 1:end) ) ;
BL_P = BL(1:tnr.npq);
BL_Q = BL (tnr.npq + 1: end);

[ax_Q, bx_Q] = max( abs( BL_Q ) );

[ax_P, bx_P] = max( abs( BL_P ) );
tnr.pq(bx_P)



