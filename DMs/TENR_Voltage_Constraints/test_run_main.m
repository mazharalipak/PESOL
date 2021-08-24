%% test_run script

clc;
clearvars;

%% Initialization Block ........................................................

data       = case300; 
testcase   = ext2int( data );                                              % Change the test case ..........................
n_stepx    = 0.5;                                                         % Newton step size alpha,   0 < alpha <= 1,  normally a small value is selected......

%% Scenario for test case (Please "Bus_fstage" values shouldn't exceed from the number of total buses in the network.)

Bus_fstage = [ 1:1:300 ];                                                   % Buses, where loadability \lambda is subjected .............    
a_bus_p    = 1.0;                                                          % On/off status of active power injection "1" ON and "0" is OFF ................................
b_bus_q    = 1.0;                                                          % On/off status of reactive power injection "1" ON and "0" is OFF ................

Vmax       = 1.10;                                                         % maximum bound on voltages magnitudes on PQ buses ......................
Vmin       = 0.90;                                                         % minimum bound on voltage  magnitudes on PQ buses .............. 

%% running TENR algorithm ..........................

s_tole = 0;                                                             
count_l = 0;

while 1

step = randsample(n_stepx, 1);                                             % random selection of optimal step size \alpha for global convergence ......................

[Vx, Vy, Slack_H, Slack_L, l, Tole, Iter, tnr, Stabi, Jx ] = tnr_testvol(testcase,step,Bus_fstage,a_bus_p,b_bus_q,Vmax,Vmin);           % Lambda is maximum Loadability, Tole is tolerance and Iter is the total iterations ........ 

count_l = count_l+1;
s_tole = Tole;

%% closing while loop ..............................

if s_tole > 1e-3
    continue
else
    break
end
end

l
Tole
Iter


