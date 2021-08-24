function [Vx, Vy, Slack_H, Slack_L, l, Tole, Iter, tnr,Stabi, J ] = tnr_testvol(testcase,step,Bus_fstage,a_bus_p,b_bus_q,Vmax,Vmin)                                                             % max Iteration bound..................................
maxIters  = 100;

%% Initial seed / germ for starting TENR Algorithm................................................................................
mpc = loadcase(testcase);
mpc = ext2int( runpf (mpc) );
%% Initializing ..................................................................................................................

tnr = tnr_init(mpc);                                                       % Initial data format
[pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);              % pv, pq buses and total no of pq and pv buses..........
nbus=1+npv+npq;

%% Dearivative of J w.r.t. Variables...............................................................................................

[AAX,BBY,CCX,CCY]=derivativexy(nbus,tnr);

[ AXY , BXY, CXY, DXY, CLY ] = deri_extJ ( AAX, BBY, CCX, CCY, tnr);

%%

Vx = real(tnr.V0);                                                           % real component of complex voltage phasor.............
Vy = imag(tnr.V0);                                                           % imag component of complex voltage phasor.............


%  Vx = (ones( size(tnr.V0)) );
%  Vy = 0.5*( ones( size(tnr.V0) ) );

% tslack = 1;

%% Voltage Feasibility constraints................................................................................................

V_max = Vmax*(ones(1+npv+npq,1));
V_min = Vmin*(ones(1+npv+npq,1));

Slack_H = ones(1+npv+npq,1);
Slack_L = 0.5*ones(1+npv+npq,1);

%% ...............................................................................................................................

l=1.0;                                                                       % Loadability parameter.................................
%% Newtpn Iterations ................................

Tole = 1;  
Iter = 1;
counter = 0;

while (Tole > 1e-3)  

[ F, delta_ss ] = F_pf( Vx,Vy,Slack_H,Slack_L,V_max,V_min,l,tnr,Bus_fstage,a_bus_p,b_bus_q);

[ J ] = tnr.J(Vx,Vy,Slack_H,Slack_L,tnr);

[ G, u, v] = G_svd(J);

[ Lastrow ] = djdx ( AAX, BBY, CCX, CCY, v, u, nbus, npq );

delta_S = [ real(delta_ss([pv;pq]));
imag(delta_ss(pq));
zeros(npv,1)
zeros(npq,1)
zeros(npq,1) ];

ExtJacobb = [ J -delta_S; Lastrow ];

DMG = [ F; G ];

%%
% TM = l - tslack^2 ;
% DM = [DMG; TM];
% ra =zeros(1,size(ExtJacobb,2)+1);
% ra(1,end-1) = 1;
% ra(1,end) = -2*tslack;
% 
% ExtJacob  = [ExtJacobb zeros( size(ExtJacobb,1), 1);ra ];

%% correction step .....................

dv = - ExtJacobb \ DMG;

%% optimal Newton step-size  ..................................

[ G_0, G_1, G_2, G_3 ] = optimal_multi( dv, ExtJacobb, AXY , BXY, CXY, DXY, CLY, DMG, nbus, npq  );

proot    = [ G_3  G_2  G_1  G_0 ];
r_x      = roots( proot );
[ ~ , bx] = min( abs( imag( r_x ) ) );
step     = abs( real( r_x( bx ) ) );

%%

Vx( [pv;pq] ) = Vx([pv;pq]) + step*dv(1:npv+npq);
Vy( [pv;pq] ) = Vy([pv;pq]) + step*dv( npv+npq+1:2*(npv+npq) );
Slack_H( pq ) = Slack_H(pq) + step*dv(2*(npv+npq)+1: 2*npv +3*npq);
Slack_L( pq ) = Slack_L(pq) + step*dv(2*npv+3*npq+1:end-1);
l             = l + step*dv(end,1);

%tslack = tslack + step*dvf;

%%

Stabi (Iter) = G;

%% Checking Tolerance.......

 Iter = Iter + 1;
 Tole=max(abs(DMG));
 counter = counter + 1;
 
 if counter == maxIters
     break;
 end
 
end


semilogy(1:1:Iter-1, (Stabi./Stabi(1,1)),'LineWidth',1.5);
xlim([-Inf Inf])
ylim([-Inf Inf])

title ('IEEE test case','Interpreter','Latex','fontsize',14);
xlabel('Newton Iterations','Interpreter',' Latex','fontsize',14);
ylabel('Stability Index','Interpreter',' Latex','fontsize',14);  

%%

if Iter >maxIters % check for non-convergence
    disp('TENR algorithm did not converge and change value of "step" '); 
else
     disp ('TENR algorithm converged')
end

end

