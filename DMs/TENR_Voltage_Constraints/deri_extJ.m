function [ AXY , BXY, CXY, DXY, CLY ] = deri_extJ ( AAX, BBY, CCX, CCY, tnr)

% Derivative of the extended Jacobian ................................................

% First w.r.t. Vx ....................

[ ~ , ~ , npv , npq ]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);              % pv, pq buses and total no of pq and pv buses..........
nbus  = 1 + npv + npq;

set_eq = 4*npq + 2*npv;

AXY = cell(nbus-1,1);

for s = 1:nbus-1
    
JX = [ cell2mat(AAX(s)) zeros( 1,(set_eq) )' ; zeros( 1,(set_eq +1) )];  
    

    AXY(s,1) = {JX};
    
end

% Second w.r.t. Vy ......................

BXY = cell(nbus-1,1);


for s = 1:nbus-1    

JY = [ cell2mat(BBY(s)) zeros( 1,(set_eq) )' ; zeros( 1,(set_eq +1) )];  

  BXY(s,1) = {JY};     
    
end


% Second w.r.t. Slack_L ......................

CXY = cell(npq,1);


for s = 1:npq    

JYC = [ cell2mat(CCX(s)) zeros( 1,(set_eq) )' ; zeros( 1,(set_eq +1) )];  

  CXY(s,1) = {JYC};     
    
end

% Second w.r.t. Slack_H  ......................

DXY = cell(npq,1);

for s = 1:npq    

JYD = [ cell2mat(CCY(s)) zeros( 1,(set_eq) )' ; zeros( 1,(set_eq +1) )];  

  DXY(s,1) = {JYD};     
    
end

%%

% Third w.r.t. lambda ...................

CLY = zeros(set_eq +1 , set_eq +1 );
CLY (end,end-1) = 0;
CLY = {CLY};

% Fourth w.r.t. ''t'' .................

% DY = zeros(2*nbus, 2*nbus);
% DY (end,end) = -2;
% DXY = {DY};


end