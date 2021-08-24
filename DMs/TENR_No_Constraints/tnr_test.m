%% Conventional Power Flow.........
function [l,Tole,Iter,Vx,Vy,tnr,Stabi,u,v, step_iter, J, delta_S]=tnr_test(mpc, step, Bus_fstage, a_bus_p, b_bus_q, maxIters)
%% Initializing ...........................................................

tnr = tnr_init( runpf(mpc) );                                                       % Initial data format

[ pv, pq, npv, npq ]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);              % pv, pq buses and total no of pq and pv buses..........

nbus = 1 + npv + npq;

%% Gradient of J ....................

[ AAX, BBY ] = derivativexy(nbus,tnr);

[ AXY , BXY, CXY, DXY ] = deri_extJ (AAX, BBY, nbus);

%% Initial guess ............................................................

% Vx = 0.5*ones(size(tnr.V0));                                               % real component of complex voltage phasor .....................
% Vy = 0.5*ones(size(tnr.V0));                                               % imag component of complex voltage phasor..............

 Vx = real( tnr.V0 );                                                          
 Vy = imag( tnr.V0 );

% tslack = 1;

l = 1.0;                                                                     % Loadability parameter.................................

%% Newtpn Iterations ................................

Tole = 1;  
Iter = 1;
counter = 0;


while (Tole > 1e-6)  

[ F, delta_injection_F ] = tnr.F( Vx, Vy, l, tnr, Bus_fstage, a_bus_p, b_bus_q);

[ J ] = tnr.J(Vx, Vy, tnr);

[ G, u, v ] = G_svd(J);

[Lastrow] = djdx ( AAX, BBY, v, u, nbus);

delta_S = [ real(delta_injection_F([pv;pq]));
imag(delta_injection_F(pq));
zeros(npv,1) ];

%% Extended Jacobain .........................................

ExtJacobb = [ J -delta_S ; Lastrow ];

DMG = [ F; G ];

%%

% TM = l - tslack^2 ;

% DM = [DMG; TM];

% ra =zeros(1,size(ExtJacobb,2)+1);
% ra(1,end-1) = 1;
% ra(1,end) = -2*tslack;

% ExtJacob  = [ExtJacobb zeros( size(ExtJacobb,1), 1); ra ];

%% correction step ..........................

dv = - ExtJacobb \ DMG;

%% optimal Newton step-size  ..................................

[ step ] = optimal_multi( dv, ExtJacobb, AXY, BXY, CXY, DXY , DMG, nbus );


%%

Vx([pv;pq]) = Vx([pv;pq]) + step*dv(1:npv+npq);
Vy([pv;pq]) = Vy([pv;pq]) + step*dv(npv+npq+1:end-1);
l = l + step*dv(end,1);

% tslack = tslack + step*dv(end,1);

%%

Stabi (Iter) = G;
step_iter (Iter) = step;

%% Checking Tolerance.......
 Iter = Iter + 1;
 Tole=max(abs(DMG));
 Toele_F (Iter) = Tole;
 counter = counter + 1;
 if counter == maxIters
     break;
 end
end


%% Convergence plot ............................

 semilogy(Stabi./Stabi(1,1),'LineWidth',1.5);
% hold on
% semilogy(Toele_F,'-.r*','LineWidth',1.5);

xlim([-Inf Inf])
ylim([-Inf Inf])

title ('IEEE test case','Interpreter','Latex','fontsize',14);
xlabel('Newton Iterations','Interpreter',' Latex','fontsize',14);
ylabel('Stability Index','Interpreter',' Latex','fontsize',14);  





% hold on
% 
% plot(step_iter, '+')
%%

if Iter > maxIters % check for non-convergence
    disp('TENR algorithm did not converge and change value of "step" '); 
else
     disp ('TENR algorithm converged')
end

end 
