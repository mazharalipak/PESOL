function [Lastrow]=djdx (AAX, BBY, v, u, nbus)

VV_R = cell(nbus-1,1);
VV_R (1:end,1) = {v};

AX = cellfun(@(x,y) x*y, AAX,VV_R, 'UniformOutput',false);
%% 

BY = cellfun(@(x,y) x*y, BBY,VV_R, 'UniformOutput',false);

%% Last row of extended Jacobain ....................................................

Lastrow = u'*sparse(cell2mat([AX' BY']));

Lastrow=[Lastrow 0];
end