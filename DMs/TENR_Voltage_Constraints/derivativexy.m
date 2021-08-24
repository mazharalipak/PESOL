% Program for computing partial_Jacob/ Partial_vx and  partial_Jacob/Partial_vy.....................
function [AAX,BBY,CCX,CCY]=derivativexy(nbus,tnr)

g=real(tnr.Ybus);
b=imag(tnr.Ybus);

[pv, pq,npv,npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);                % pv, pq buses and total no of pq and pv buses..........

%la=tnr.branch(:,1);                                                        % From bus i to j........................................
%lb=tnr.branch(:,2);                                                        % From bus j to i........................................
bus_length=1:nbus;

order=[pv;pq];                                                             % Buses order............................................
%% Derivative of The Jacobian matrix w.r.t. VX

AAX=cell(nbus-1,1);

for s=1:nbus-1
    
    J1=sparse(order(s),bus_length,g(order(s),bus_length),nbus,nbus);
    J1= (J1-diag(diag(J1)))+ sparse(1:nbus,1:nbus,g(order(s),bus_length),nbus,nbus);
    J1(order(s),order(s))=2*g(order(s),order(s));    
        
    
%%  J2/////////////////////
   
   J2=sparse(order(s),bus_length,-b(order(s),bus_length),nbus,nbus);
   J2=(J2-diag(diag(J2)))  + sparse(1:nbus,1:nbus,b(order(s), bus_length),nbus,nbus);                  % diagonal...........................
   J2(order(s),order(s))=0;    


   J13 = sparse(1+npv+npq, 1+npv+npq);
   
%%  J3//////////////////////
  
   J3=sparse(order(s),bus_length,-b(order(s),bus_length),nbus,nbus);
   J3= (J3-diag(diag(J3))) + sparse(1:nbus,1:nbus,-b(order(s), bus_length),nbus,nbus);                  % diagonal...........................
   J3(order(s),order(s))=-2*b(order(s),order(s));

%%  J4//////////////////////
  
   J4=sparse(order(s),bus_length,-g(order(s),bus_length),nbus,nbus);
   J4=(J4-diag(diag(J4))) + sparse(1:nbus,1:nbus,g(order(s),bus_length),nbus,nbus);  % diagonal...........................
   J4(order(s),order(s))=0;
  
%%  J5///////////////////////
  
    J5=sparse(order(s),order(s),2,nbus,nbus); 

%%  J6///////////////////////////

    J6=sparse(1+npv+npq, 1+npv+npq);

  
%%


  
%%
  J_VX=sparse([J1([pv;pq],[pv;pq]) J2([pv;pq],[pv;pq])  J13([pv;pq],pq)     J13([pv;pq],pq) ; 
      J3(pq,[pv;pq])               J4(pq,[pv;pq])       J13(pq,pq)          J13(pq,pq); 
      J5(pv,[pv;pq])               J6(pv,[pv;pq])       J13(pv,pq)          J13(pv,pq)
      J5(pq,[pv;pq])               J6(pq,[pv;pq])       J13(pq,pq)          J13(pq,pq)
      J5(pq,[pv;pq])               J6(pq,[pv;pq])       J13(pq,pq)          J13(pq,pq)]);

  AAX(s,1)={J_VX};

end
 


%% Derivative of The Jacobian matrix W.R.T VY
 
BBY=cell(nbus-1,1);

for s=1:nbus-1
    
    
    J1=sparse(order(s),bus_length,b(order(s),bus_length),nbus,nbus);
    J1=(J1-diag(diag(J1))) + sparse(1:nbus,1:nbus,-b(order(s),bus_length),nbus,nbus);
    J1(order(s),order(s))=0;    
    
    
    J13 = sparse(1+npv+npq, 1+npv+npq);

%%  J2/////////////////////

  
   J2=sparse(order(s),bus_length,g(order(s),bus_length),nbus,nbus);
   J2=(J2-diag(diag(J2))) + sparse(1:nbus,1:nbus,g(order(s),bus_length),nbus,nbus);             % diagonal...........................
   J2(order(s),order(s))=2*g(order(s),order(s));
   
    

%%  J3//////////////////////

  J3=sparse(order(s),bus_length,g(order(s),bus_length),nbus,nbus);
  J3=(J3-diag(diag(J3))) + sparse(1:nbus,1:nbus,-g(order(s),bus_length),nbus,nbus);  % diagonal...........................
  J3(order(s),order(s))=0;
  

%%  J4//////////////////////
  
  J4=sparse(order(s),bus_length,-b(order(s),bus_length),nbus,nbus);
  J4= (J4-diag(diag(J4)))+ sparse(1:nbus,1:nbus,-b(order(s), bus_length),nbus,nbus);  % diagonal...........................
  J4(order(s),order(s))=-2*b(order(s),order(s));
  

%% J5/////////////////////

  J5=sparse(1+npv+npq, 1+npv+npq);

%%  J6//////////////////////
 
  J6=sparse(order(s),order(s),2,nbus,nbus);

  %%
  J_VY=sparse([J1([pv;pq],[pv;pq]) J2([pv;pq],[pv;pq])  J13([pv;pq],pq)     J13([pv;pq],pq) ; 
      J3(pq,[pv;pq])               J4(pq,[pv;pq])       J13(pq,pq)          J13(pq,pq); 
      J5(pv,[pv;pq])               J6(pv,[pv;pq])       J13(pv,pq)          J13(pv,pq)
      J5(pq,[pv;pq])               J6(pq,[pv;pq])       J13(pq,pq)          J13(pq,pq)
      J5(pq,[pv;pq])               J6(pq,[pv;pq])       J13(pq,pq)          J13(pq,pq)]);

  BBY(s,1)={J_VY};

end   

%%

CCX=cell(npq,1);

for s=1:npq
    
    J1=sparse(1+npv+npq,1+npv+npq);
    
    
    Jxh=sparse(pq(s),pq(s),2,nbus,nbus);
    
    
    J_SH=sparse([J1([pv;pq],[pv;pq]) J1([pv;pq],[pv;pq])  J1([pv;pq],pq)     J1([pv;pq],pq) ; 
      J1(pq,[pv;pq])                 J1(pq,[pv;pq])       J1(pq,pq)          J1(pq,pq); 
      J1(pv,[pv;pq])                 J1(pv,[pv;pq])       J1(pv,pq)          J1(pv,pq)
      J1(pq,[pv;pq])                 J1(pq,[pv;pq])       Jxh(pq,pq)         J1(pq,pq)
      J1(pq,[pv;pq])                 J1(pq,[pv;pq])       J1(pq,pq)          J1(pq,pq)]);
   
  CCX(s,1)={J_SH};
    
end


%% 

CCY=cell(npq,1);

for s=1:npq
    
    J1 = sparse(1+npv+npq,1+npv+npq);
    
    Jxh = -sparse(pq(s),pq(s),2,nbus,nbus);
    
    
    J_SH = sparse( [J1([pv;pq],[pv;pq]) J1([pv;pq],[pv;pq])  J1([pv;pq],pq)     J1([pv;pq],pq) ; 
      J1(pq,[pv;pq])                 J1(pq,[pv;pq])       J1(pq,pq)             J1(pq,pq); 
      J1(pv,[pv;pq])                 J1(pv,[pv;pq])       J1(pv,pq)             J1(pv,pq)
      J1(pq,[pv;pq])                 J1(pq,[pv;pq])       J1(pq,pq)             J1(pq,pq)
      J1(pq,[pv;pq])                 J1(pq,[pv;pq])       J1(pq,pq)             Jxh(pq,pq)] );
   
  CCY(s,1)={J_SH}; 
    
end



end