function [ RP_pred ] = Overpressure(a,b,c,Pves_pred,Pop,Buoy,RP_pred,RP0,RPmax )
%OVERPRESSURE Summary of this function goes here
%   Detailed explanation goes here

RP_op_beta=(8/10000)*(1-exp(-a/Pves_pred));
RP_op_max=RP0+((RPmax-RP0)*exp(-b/Pves_pred));
RP_op_min=RPmax-((RPmax-RP0)*exp(-c*Pves_pred));

RP_op=RP_op_max-((RP_op_max-RP_op_min)*exp(-RP_op_beta*Pves_pred));
alpha=RP_pred-RP_op;
RP_pred_op=RP_op_max-((RP_op_max-RP_op_min)*exp(-RP_op_beta*(Pves_pred-Pop-Buoy)));
RP_pred=RP_pred_op+alpha;

end

