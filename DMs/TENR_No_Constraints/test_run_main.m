%% test_run script

clc;
clearvars;

maxIters  = 60;            % maximum number of iterations ....................................

%% Initialization Block ...........................

testcase   = case14;

%% Initial seed / germ for starting TENR Algorithm.........................

mpc = ext2int( ( loadcase( testcase ) ) );

%% Newton step size alpha,   0 < alpha <= 1,  normally a small value is selected......

step    = 0.9;                                                 

%% Scenario for test case (Please "Bus_fstage" values shouldn't exceed from the number of total buses in the network.)

Bus_fstage = [ 1,2,3,4,5,6,7,8,9,10,11,12,13,14];                          % Buses, where loadability \lambda is subjected .................     
a_bus_p    = 1.0;                                                          % On/off status of active power injection "1" ON and "0" is OFF ................
b_bus_q    = 0.0;                                                          % On/off status of reactive power injection "1" ON and "0" is OFF ................

%% running TENR algorithm ..........................

[lambda, Tole, Iter, Vx, Vy, tnr, Stabi,ur,vr, step_iter, Jxx,delta_S]=tnr_test(mpc, step, Bus_fstage, a_bus_p, b_bus_q, maxIters);   

 lambda
 Iter
 Tole
 


