function [J] = J_pf(Vx,Vy,tnr)
%J_pf Calculate the residual of power flow Jacobian
%   Detailed explanation goes here

    [pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);              % pv, pq buses and total no of pq and pv buses..........

    G_y = real(tnr.Ybus);
    B_y = imag(tnr.Ybus);
    
    nubs = 1+npv+npq;
    
  %% Partial P/ partial vx................... P(x)...........
    
    J11 = (Vx.*G_y) + (Vy.*B_y);                                           % Off diagonal..............
   
  % Jexp = diag(Vx)*G_y + diag(Vy)*B_y;
  
%    Jexp =  2*( diag(Vx)*diag(diag(G_y))) +   (  ( G_y*Vx- ( diag(Vx)*diag(diag(G_y)) )) - ( B_y*Vy- ( diag(Vy)*diag(diag(B_y))))  );
%    JexpL = full ( (sparse(1:nubs,1:nubs, (2*Vx.*diag(G_y)) +   (  (G_y*Vx- (Vx.*diag(G_y))) - (B_y*Vy- (Vy.*diag(B_y)))  ), nubs,nubs )) );
%      
    
    J11= (J11-diag(diag(J11))) + (sparse(1:nubs,1:nubs, (2*Vx.*diag(G_y)) +   (  (G_y*Vx- (Vx.*diag(G_y))) - (B_y*Vy- (Vy.*diag(B_y)))  ), nubs,nubs ));
    
    %% Partial P/ partial vy...................
    
    J12 =(-Vx.*B_y)+(Vy.*G_y);
    J12 =(J12- diag(diag(J12))) + sparse(1:nubs,1:nubs, (2*Vy.*diag(G_y))  + ( (B_y*Vx -(Vx.*diag(B_y))) + (G_y*Vy -(Vy.*diag(G_y)) ) ), nubs,nubs  );
    
   %% Partial Q/ partial vx...................Q(x)..............
   
    J21 = (-Vx.*B_y) +(Vy.*G_y); 
    J21 = (J21 - diag(diag(J21))) + sparse(1:nubs,1:nubs, (-2*Vx.*diag(B_y)) - ( (B_y*Vx -(Vx.*diag(B_y))) + (G_y*Vy  -(Vy.*diag(G_y)))  )   , nubs,nubs );
   
   %% Partial Q/ partial vy...................
   
   J22 = (-Vx.*G_y) - (Vy.*B_y);
   J22 = (J22 -diag(diag(J22))) + sparse(1:nubs,1:nubs, (- 2*Vy.*diag(B_y))  + ((G_y*Vx -(Vx.*diag(G_y))) - (B_y*Vy- (Vy.*diag(B_y)))), nubs,nubs  );
    
   %% Partial |V|/ partial vx.................|V|^2............... 
   
   J31=sparse(1:nubs,1:nubs,(2*Vx),nubs,nubs);
   
   % Partial |V|/ partial vx.....................................
   
   J32=sparse(1:nubs,1:nubs,(2*Vy),nubs,nubs);


   
   %% Jacobian Matrix....................................................
   

J=[ J11([pv;pq],[pv;pq]) J12([pv;pq],[pv;pq]); 
    J21(pq,[pv;pq]) J22(pq,[pv;pq]); 
    J31(pv,[pv;pq]) J32(pv,[pv;pq])];

end

