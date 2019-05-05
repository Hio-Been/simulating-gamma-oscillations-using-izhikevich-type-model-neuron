function [kernel,tau] = get_dend_kernel(t)
% alpha = 1/10; f = @factorial; 
% tau = 1:length(t);
% kernel = [1/max([exp(-alpha*tau).*(  ((alpha*tau).^5)/f(5) - ((alpha*tau).^7)/f(7)  )] )] ...
%     * [exp(-alpha*tau).*(  ((alpha*tau).^5)/f(5) - ((alpha*tau).^7)/f(7)  )];
alpha = 1/5; f = @factorial; 
tau = 1:length(t);
p1 = 1; p2 = 3;
kernel = [1/max([exp(-alpha*tau).*(  ((alpha*tau).^p1)/f(p1) - ((alpha*tau).^p2)/f(p2)  )] )] ...
    * [exp(-alpha*tau).*(  ((alpha*tau).^p1)/f(p1) - ((alpha*tau).^p2)/f(p2)  )];
