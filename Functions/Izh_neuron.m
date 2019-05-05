function [ dV ] = Izh_neuron(t,V,a,b,I)
% [ dV ] = Izh_neuron(t,V,a,b,c,d,I)
% 
% Dynamics for Izhikevich neurons
% Used for calling from RK4.

% from 2003 paper
dV(1,:)=0.04*V(1,:).^2+5*V(1,:)+140-V(2,:)+I;
dV(2,:)=a.*(b.*V(1,:)-V(2,:));
