function [Vx, Vy, Slack_H, Slack_L, Lamda_1, Lamda_2] = sortbuses(Vx, Vy, Slack_H, Slack_L, Lamda_1, Lamda_2, tnr, vector, step)

[pv, pq, npv, npq]  = deal(tnr.pv, tnr.pq, tnr.npv, tnr.npq);              % pv, pq buses and total no of pq and pv buses..........

Vx([pv;pq]) = Vx([pv;pq]) + step *vector(1:npv+npq);
 
Vy([pv;pq]) = Vy([pv;pq]) + step *vector(npv+npq+1:2*(npv+npq));

Slack_H(pq) = Slack_H(pq) + step *vector(2*(npv+npq)+1: 2*npv +3*npq);

Slack_L(pq) = Slack_L(pq) + step *vector(2*npv+3*npq+1:end-2);

Lamda_1 = Lamda_1 + step *vector(end-1);
 
Lamda_2 = Lamda_2 + step *vector(end,1);




end 