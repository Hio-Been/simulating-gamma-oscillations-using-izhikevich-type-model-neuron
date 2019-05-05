addpath Functions
%% Set Initial Parameters
dt = .5;
tmax = 1000;

param.time = dt:dt:tmax;
param.ths = 30;
param.v_init = -65;
param.u_init = -10;
param.inject = 5 * ones(size(param.time));
param.dend_input = []; % No dendritic input

for neuron = 1:2
    figure(neuron); clf;
    
    %% Set parameters a, b, c, d
    if neuron==1
        param.abcd = [ .2 .20 -65 2]; % Given typical values
    else
        param.abcd = [ .02, .2, -50, 2 ]; % Bursting
    end
    
    
    %% Get activity
    [V, t, Spikes, FR ] =get_Izhik( param );
    
    %% Plot figures
    subplot(1,2,1); hold off;
    plot( t, V, 'k' );
    xlabel('Time (ms)'); ylabel('V(t) (mV)');
    hold on;
    plot( dt*find(Spikes), ones([size(find(Spikes))])*200, 'r+' );
    axis tight;
    axis([0 max(param.time) -80 220] );
    title(['Result of [a,b,c,d]=[' num2str(param.abcd(1)) ',' num2str(param.abcd(2)) ','...
        num2str(param.abcd(3)) ',' num2str(param.abcd(4)) ']' ]);
    % enhance_fig_visibility( gca );
    
    subplot(1,2,2); hold off;
    histogram( diff(find(Spikes)), [0:2:70], 'FaceColor', 'k' );
    title(['ISI Histogram (mean FR = ' num2str(FR) ' Hz)']);
    xlabel('Inter-Spike Interval (ms)'); ylabel('#Observation');
    % enhance_fig_visibility( gca );
    disp(['Mean ISI = ' num2str(nanmean(diff(find(Spikes)))) ' (' num2str(nanstd(diff(find(Spikes)))) ') ms' ])
    
    % saveas(gcf, ['Fig_' num2str(mfilename) '.png'] );
    
end