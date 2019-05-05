function [ V, tvec, Spikes, FR ] = get_Izhik( param, plotOption )
if nargin < 2
    plotOption=false; end; if nargin < 1
    abcd = [ .02, .20, -55, 2 ]; param = [];
    param.abcd = abcd; dt = 1; tmax = 1000;
    param.time = dt:dt:tmax; param.inject = 0 * ones(size(param.time));
    param.dend_input = [zeros([1,250]), ones([1,10]), zeros([1,740]) ];
    param.ths = 30; param.v_init = -65; param.u_init = -10; plotOption= true;
end

% Temporal Kernel
if ~isempty( param.dend_input)
    alpha = 1/10; f = @factorial;
    tau = 1:length(param.time);
    Dt_raw = exp(-alpha*tau).*(  ((alpha*tau).^5)/f(5) - ((alpha*tau).^7)/f(7)  ) ;
    C = 1/max(Dt_raw );
    Dt = C * Dt_raw;
end

% plot( tau, Dt , 'b', 'LineWidth', 2);
% set(gca, 'XDir','reverse')
% xlim([0 300]);
% xlabel('\tau (ms)');
% ylabel('D(\tau)' );
% title('Kernel');
% enhance_fig_visibility( gca ) ;


% Set parameters
a = param.abcd(1);
b = param.abcd(2);
c = param.abcd(3);
d = param.abcd(4);

% Set initial values
v=param.v_init;
u=param.u_init;

% Handy variables
V = [];
U = [];
Spikes = logical(zeros([1,length(param.time)]));

for t2= 1:length(param.time)          % simulation of 1000 ms
    t = param.time(t2);
    % With dendritic Input
    if ~isempty( param.dend_input )
        D = param.dend_input(t:-1:1);
        D_with_kernel = sum(D .* Dt(1:t));
        I = param.inject(t) + D_with_kernel;
        
        % Without dendritic input
    else
        I = param.inject(t2);
    end
    
    v = v+ 0.5*(0.04*v.^2+5*v+140-u+I);
    u = u+a.*(b.*v-u);
    
    V(t2) = v;
    U(t2) = u;
    
    if v > param.ths
        v=c;
        u=u+d;
        Spikes(t2) = true;
    end
    
end

nSpikes = sum(Spikes);
FR = nSpikes / (max( param.time )/1000);
tvec = param.time;

if plotOption
    subplot(1,3,1); hold off;
    plot( tvec, V, 'k' );
    xlabel('Time (ms)'); ylabel('V(t) (mV)');
    hold on;
    plot( find(Spikes), ones([size(find(Spikes))])*200, 'rs' );
    axis tight;
    axis([0 max(param.time) -80 220] );
    title(['V(t), FR = ' num2str(FR) ' Hz']);
    enhance_fig_visibility( gca );
    
    subplot(1,3,2); hold off;
    plot( tvec, U, 'k' );
    xlabel('Time (ms)'); ylabel('U(t) (a. u.)');
    hold on;
    plot( find(Spikes), ones([size(find(Spikes))])*2, 'rs' );
    axis tight;
    axis([0 max(param.time) -9 3] );
    title(['U(t)']);
    enhance_fig_visibility( gca );
    
    subplot(1,3,3); hold off;
    histogram( diff(find(Spikes)), 15, 'FaceColor', 'k' );
    title(['ISI Histogram']);
    xlabel('ISI (ms)'); ylabel('#Observation');
    enhance_fig_visibility( gca );
end
% ylim([-100 200]);


