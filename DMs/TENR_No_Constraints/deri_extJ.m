function [ AXY , BXY, CXY, DXY ] = deri_extJ (AAX, BBY, nbus)

%%%%%%%%%%% Derivative of the extended Jacobian ................................................

% First w.r.t. Vx ....................

AXY = cell(nbus-1,1);

for s = 1:nbus-1
    
    
%   JX = [ cell2mat(AAX(s)) zeros(2,(2*nbus-2))' ; zeros(2,(2*nbus))];  

JX = [ cell2mat(AAX(s)) zeros(1,(2*nbus-2))' ; zeros(1,(2*nbus -1))];  
    

    AXY(s,1) = {JX};
    
end


% Second w.r.t. Vy ......................

BXY = cell(nbus-1,1);


for s = 1:nbus-1
    
%   JY = [ cell2mat(BBY(s)) zeros(2,(2*nbus-2))' ; zeros(2,(2*nbus))];  
    

JY = [ cell2mat(BBY(s)) zeros(1,(2*nbus-2))' ; zeros(1,(2*nbus -1))];  

  BXY(s,1) = {JY};  
    
    
end

% Third w.r.t. lambda ...................

CY = zeros(2*nbus-1, 2*nbus-1);
CY (end,end-1) = 0;
CXY = {CY};

% Fourth w.r.t. ''t'' .................

DY = zeros(2*nbus, 2*nbus);
DY (end,end) = -2;
DXY = {DY};


end