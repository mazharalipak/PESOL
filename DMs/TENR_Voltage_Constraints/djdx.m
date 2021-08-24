function [ Lastrow ] = djdx ( AAX, BBY, CCX, CCY, v, u, nbus, npq )

%%

VV_R = cell(nbus-1,1);
VV_R (1:end,1) = {v};

VV_PQ = cell(npq,1);
VV_PQ (1:end,1) = {v};

%%

AX = cellfun(@(x,y) x*y, AAX, VV_R, 'UniformOutput',false);

BY = cellfun(@(x,y) x*y, BBY, VV_R, 'UniformOutput',false);

CX = cellfun(@(x,y) x*y, CCX, VV_PQ, 'UniformOutput',false);

CY = cellfun(@(x,y) x*y, CCY, VV_PQ, 'UniformOutput',false);

%%

Lastrow = u'*sparse(cell2mat([AX' BY' CX' CY']));

Lastrow=[Lastrow 0];

end