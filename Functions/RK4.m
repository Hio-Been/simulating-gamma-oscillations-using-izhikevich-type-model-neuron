function [yNew, tNew]=RK4(t,y,h,funcName,varargin)
% [yNew, tNew]=RK4(t,y,h,funcName,varargin)
% 
% where dy/dt=funcName(t,y,varargin)

k1=feval(funcName,t,y,varargin{:});
k2=feval(funcName,t+h/2,y+h/2*k1,varargin{:});
k3=feval(funcName,t+h/2,y+h/2*k2,varargin{:});
k4=feval(funcName,t+h,y+h*k3,varargin{:});

yNew=y+h/6*(k1+2*k2+2*k3+k4);

tNew=t+h;