% Program for computing partial_Jacob/ Partial_vx and  partial_Jacob/Partial_vy.....................
function [AAX,BBY]=derivativexy(nbus,tnr)

g = real(tnr.Ybus);
b = imag(tnr.Ybus);

[pv, pq,npv,~]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);                % pv, pq buses and total no of pq and pv buses..........

bus_length = 1:nbus;

order = [pv;pq];                                                             % Buses order............................................

%% Derivative of The Jacobian matrix w.r.t. VX

AAX=cell(nbus-1,1);

for s = 1:nbus-1
    
    J1 = sparse(order(s),bus_length,g(order(s),bus_length),nbus,nbus);
    J1 = (J1-diag(diag(J1)))+ sparse(1:nbus,1:nbus,g(order(s),bus_length),nbus,nbus);
    J1(order(s),order(s)) = 2*g(order(s),order(s));    
        
    
%%  J2/////////////////////
   
   J2 = sparse(order(s),bus_length,-b(order(s),bus_length),nbus,nbus);
   J2 = (J2-diag(diag(J2)))  + sparse(1:nbus,1:nbus,b(order(s), bus_length),nbus,nbus);                  % diagonal...........................
   J2(order(s),order(s)) = 0;    


%%  J3//////////////////////
  
  J3 = sparse(order(s),bus_length,-b(order(s),bus_length),nbus,nbus);
  J3 = (J3-diag(diag(J3))) + sparse(1:nbus,1:nbus,-b(order(s), bus_length),nbus,nbus);                  % diagonal...........................
  J3(order(s),order(s)) = -2*b(order(s),order(s));

%%  J4//////////////////////
  
  J4 = sparse(order(s),bus_length,-g(order(s),bus_length),nbus,nbus);
  J4 = (J4-diag(diag(J4))) + sparse(1:nbus,1:nbus,g(order(s),bus_length),nbus,nbus);  % diagonal...........................
  J4(order(s),order(s)) = 0;
  
%%  J5///////////////////////
  
  J5 = sparse(order(s),order(s),2,nbus,nbus); 

%%  J6///////////////////////////

  J6 = sparse(npv,max(size(order)));

%%
  J_VX = sparse([J1([pv;pq],[pv;pq]) J2([pv;pq],[pv;pq]); 
      J3(pq,[pv;pq]) J4(pq,[pv;pq]); 
      J5(pv,[pv;pq]) J6]);

  AAX(s,1) = {J_VX};

end
 


%% Derivative of The Jacobian matrix W.R.T VY
 
BBY=cell(nbus-1,1);

for s=1:nbus-1
    
    
    J1 = sparse(order(s),bus_length,b(order(s),bus_length),nbus,nbus);
    J1 = (J1-diag(diag(J1))) + sparse(1:nbus,1:nbus,-b(order(s),bus_length),nbus,nbus);
    J1(order(s),order(s)) = 0;    
    
    
    
%%  J2/////////////////////

  
   J2 = sparse(order(s),bus_length,g(order(s),bus_length),nbus,nbus);
   J2 = (J2-diag(diag(J2))) + sparse(1:nbus,1:nbus,g(order(s),bus_length),nbus,nbus);             % diagonal...........................
   J2(order(s),order(s)) = 2*g(order(s),order(s));
   
    

%%  J3//////////////////////

  J3 = sparse(order(s),bus_length,g(order(s),bus_length),nbus,nbus);
  J3 = (J3-diag(diag(J3))) + sparse(1:nbus,1:nbus,-g(order(s),bus_length),nbus,nbus);  % diagonal...........................
  J3(order(s),order(s)) = 0;
  

%%  J4//////////////////////
  
  J4 = sparse(order(s),bus_length,-b(order(s),bus_length),nbus,nbus);
  J4 = (J4-diag(diag(J4)))+ sparse(1:nbus,1:nbus,-b(order(s), bus_length),nbus,nbus);  % diagonal...........................
  J4(order(s),order(s)) = -2*b(order(s),order(s));
  

%% J5/////////////////////

  J5 = sparse(npv,max(size(order)));

%%  J6//////////////////////
 
  J6 = sparse(order(s),order(s),2,nbus,nbus);

  %%
  J_VY = sparse([J1([pv;pq],[pv;pq]) J2([pv;pq],[pv;pq]); 
      J3(pq,[pv;pq]) J4(pq,[pv;pq]); 
      J5 J6(pv,[pv;pq])]);

  BBY(s,1) = {J_VY};

end   
end