% optimal step size ...............

function [ G_0, G_1, G_2, G_3 ] = optimal_multi(delta_x, Jacob, AAX, BBY, F_x, nbus )
% delta_x               ..............
% Jacob                 ..............
% AAX                   .............
% BBY                   .............
% F_x                   ............

a_x = F_x;
b_x = Jacob*delta_x; ......................

VV_R = cell(nbus-1,1);
VV_R (1:end,1) = {delta_x};


%%

AX = cellfun(@(x,y) x*y, AAX, VV_R, 'UniformOutput',false);

BY = cellfun(@(x,y) x*y, BBY, VV_R, 'UniformOutput',false);

Lastrow = 0.5*( delta_x'*sparse(cell2mat([AX' BY'])) );

c_x = Lastrow'; 

%%

G_0 = -a_x'*b_x;
G_1 = b_x'*b_x + 2*a_x'*c_x;
G_2 = -3*b_x'*c_x;
G_3 = 2*( c_x'*c_x);


end
