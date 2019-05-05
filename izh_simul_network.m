clear all
warning off
input = [];
addpath Functions
m = 'bicubic';

F_input = .1;
F_noise = 500;
% gen_input = 0;
syn_factor = 1;
pause(5);

%% (1) Defining Neuronal Network
% Parameter setting
input.dt = .5; % 2000 Hz
input.tLim = [0,3000]; % Total 1000 ms
input.t = input.dt:input.dt:input.tLim(end);
dummy_list = 1:length(input.t);

N = 100; % Total Neuron Number

nNeuron = round( [ .8, .2*[.8 .1 .1] ] * N );
vidname = 'with_interneuron';

%nNeuron = round( [ .97, .01 .01 .01 ] * N );
%vidname = 'without_interneuron';

cellTypes = {'E_{PYR}', 'I_{PV}', 'I_{SST}', 'I_{VIP}'};
nType = length(cellTypes);
input.nNeuron = nNeuron;
input.cellTypes = cellTypes;

input.EI = [ ones([1,nNeuron(1)]), zeros([1,N-nNeuron(1)]) ];
EI_index = find(input.EI);
input.cell_type = [ ...
    1 * ones([1,nNeuron(1)]), ...
    2 * ones([1,nNeuron(2)]) ,...
    3 * ones([1,nNeuron(3)]) ,...
    4 * ones([1,nNeuron(4)]) ,...
    ];

s0 = .0001;% practically zero
s1 = .010; % weak connection
s2 = .230; % moderate
s3 = .330; % strong
s_jitter = [-1 1]*.0001; % Uniform

weights = [ ...
    % PYR PV SST VIP
    s1  s2  s0  s0; % PYR
    s2  s3  s0  s0; % PV
    s2  s1  s0  s0; % SST
    s0  s1  s2  s0  % VIP
    ];

synapse = zeros(N);

% Getting weight matrix template (exponentially-decreasing)
weight_template= zeros(max(nNeuron));
expF = 5; % The higher, the less local connection
weight_decreasing = exp(linspace(expF, 0, max(nNeuron))); weight_decreasing=weight_decreasing/max(weight_decreasing);
for x = 1:max(nNeuron), for y = 1:max(nNeuron)
        dist = abs(x-y)+1;
        weight_template(x,y)=weight_decreasing(dist);
    end, end

% Obtain synapse matrix - reflecting connection-specific factor (weights)
for pre = 1:nType, for post = 1:nType
        synapse(find(input.cell_type==pre), find(input.cell_type==post)) = ...
            (ones([nNeuron(pre), nNeuron(post)])*weights(pre,post))... Given weight
            .* imresize(weight_template,[nNeuron(pre), nNeuron(post)], m)... Distance
            + (rand([nNeuron(pre), nNeuron(post)])*range(s_jitter)+s_jitter(1)); ... Jitter (noise)
    end, end; synapse(find(synapse<0))=synapse(find(synapse<0))*0;
synapse = syn_factor * synapse;
input.S = synapse;
imagesc(synapse)

%% (2) Dendritic input generator (Locally connected)
try input_dendrite;
catch;
    fig = figure(5);
    kernels = [];
    % Input only for PYR and PV
    thalamic_input =find( [ ones([1, sum( nNeuron(1:2) )]), zeros([1, sum( nNeuron(3:4) )])]) ;
    input_target = [ ones([1, sum(nNeuron(1:2))]), zeros([1, sum(nNeuron(3:4))]) ];
    nTarget = sum(input_target);
    
    % Temporal Kernel
    [kernels.temporal,tau] = get_dend_kernel( input.t );
    
    % Get 2D Spatial Weight Kernel
    input_unit_size = [20];
    expF = 1; % The higher, the less local connection
    distance = 1:input_unit_size/2;
    kernel_func = get_exp_decreasing_kernel( distance, expF );
    [x, y] = meshgrid(1:input_unit_size, 1:input_unit_size);
    dist = sqrt((x-input_unit_size/2).^2 + (y-input_unit_size/2).^2);
    kernels.spatial_template = reshape(kernel_func(dist(:)+1), size(dist));
    
    nSkip = 1;
    left = 1:nSkip:nNeuron(1)*nSkip;
    right = left+input_unit_size-1;
    k = zeros([ input_unit_size, right(end) ]);
    
    kernels.spatial_mask = single([]);
    % PYR
    for cellIdx = 1:nNeuron(1)
        kernels.spatial_mask(:,:,cellIdx) = k;
        kernels.spatial_mask(:,left(cellIdx):right(cellIdx),cellIdx) = kernels.spatial_template;
    end
    % PV
    scale = round(linspace( 1, nNeuron(1), nNeuron(2) ));
    for cellIdx = nNeuron(1)+1:nNeuron(1)+nNeuron(2)
        kernels.spatial_mask(:,:,cellIdx) = kernels.spatial_mask(:,:,scale(cellIdx-nNeuron(1))) ;
    end
    input_size = [size( kernels.spatial_mask,1), size(kernels.spatial_mask,2)];
    if false
        figure(5);
        for cellIdx = 1:3:nNeuron(1)+nNeuron(2)
            imagesc( kernels.spatial_mask(:,:,cellIdx));
            colormap gray
            title(['Input kernel, Cell#' num2str(cellIdx)]);
            drawnow;
        end
    end
    
    % Getting  poisson spikes
    input_Hz = 100;
    poisson_rate=input_Hz*input.dt/1000;
    input_spike = logical(rand([input_size, length(input.t)]) < poisson_rate); % Poisson
    
    % Get dendrite input
    input_dendrite = single(zeros([ length(input.t) input_size(1) ]));
    kernel_effective_size = 80;
    kernel_short = reshape(kernels.temporal(1:kernel_effective_size), [1 1 kernel_effective_size]);
    kernel_repmat = repmat(kernel_short, [input_size(1), input_size(2), 1]);
    for tIdx = 1:length(input_spike)
        input_index = tIdx-kernel_effective_size+1:tIdx;
        input_index = input_index(find(input_index>0));
        dend_input = sum( input_spike(:,:,input_index) ...
            .* kernel_repmat(:,:,length(input_index):-1:1), 3);
        for cellIdx = 1:length(thalamic_input)
            conv_input = dend_input .* kernels.spatial_mask(:,:,cellIdx);
            input_dendrite( tIdx, cellIdx ) = sum(conv_input(:));
        end
        input_dendrite( tIdx, cellIdx+1:N) = zeros([1,N-cellIdx]);
        
        %% Visualization
        if mod(tIdx, 200) == 0
            subplot(2,2,[1 3]);
            plot_multi_ch( input_dendrite', input.t );
            xlabel('Time (ms)');
            ylabel('Cell #');
            title(['t = ' num2str(tIdx)])
            
            subplot(4,2,2);
            imagesc( input_spike(:,:,tIdx) );
            colormap gray
            title(['Spike input (t = ' num2str(tIdx) ')' ]);
            
            subplot(4,2,4);
            ex_cell = randi(length(thalamic_input));
            imagesc( kernels.spatial_mask(:,:,ex_cell));
            colormap gray
            title(['Spatial kernel (ex, #' num2str(ex_cell) ')' ]);
            
            subplot(4,2,6);
            hold off;
            plot(tau, kernels.temporal, 'k');
            xlim([0 60]);
            set(gca, 'XTick', [0:10:300]);
            xlabel('\tau (ms)');
            ylabel('weight');
            title('Temporal kernel');
            
            subplot(4,2,8);
            imagesc( corrcoef( input_dendrite ))
            title(['Similarity of input']);
            colormap gray;
            xlabel('V1 neuron #'); ylabel('V1 neuron #');
            cb=colorbar;
            hold off;
            ylabel(cb, 'Corr coeff (r)');
            
            suptitle(['Generating dendrite input (Poisson, ' num2str(input_Hz) 'Hz)' ])
            drawnow;
        end
        
    end
end
input.I = input_dendrite;



%% Etc.
noise_range = [ -1 1]*.05;
input.noise = rand([ length(input.t), N ]) * diff(noise_range) + noise_range(1);
input.tau_G=[2; 8];
input.TAHpars=[.5, 1];
input.delay = 1; delayFlag = 1;
E_idx = find(input.EI);
I_idx = find(~input.EI);


%% (3) Model visualization
fig=figure(1); clf;

sp1=subplot(2,2,1); hold off;
cellType_N = {}; for ij = 1:4, cellType_N{ ij } = [ cellTypes{ij} ' (n=' num2str(nNeuron(ij)) ')' ]; end
pie( nNeuron, cellType_N )
% bar(1:4, nNeuron, 'k');
ylabel('#Cell');
set(gca , 'XTick', 1:4, 'XTIckLabel', cellTypes);
set(gca, 'LineWidth', 2, 'Box', 'off', 'FontSize', 13)

sp2=subplot(2,2,2); hold off;
imagesc( input.S ); %axis xy;
set(gca, 'YTick', [ nNeuron(1)*.5, sum(nNeuron(1))+nNeuron(2)*.5, sum(nNeuron(1:2))+nNeuron(3)*.5, sum(nNeuron(1:3))+nNeuron(4)*.5], ...
    'YTickLabel', cellTypes);
set(gca, 'XTick', [ nNeuron(1)*.5, sum(nNeuron(1))+nNeuron(2)*.5, sum(nNeuron(1:2))+nNeuron(3)*.5, sum(nNeuron(1:3))+nNeuron(4)*.5], ...
    'XTickLabel', cellTypes);
colormap gray; cb=colorbar;
caxis([0 max(synapse(:))]);
ylabel(cb,['weights']);
xlabel('Post-synapse');
ylabel('Pre-synapse');
set(gca, 'LineWidth', 2, 'Box', 'off', 'FontSize', 13)


subplot(2,2,3); hold off;
% histogram(input.I(find(input.I>0)))
xlabel('Input intensity');
set(gca, 'LineWidth', 2, 'Box', 'off', 'FontSize', 13)

subplot(2,2,4); hold off;
histogram(input.noise(find(input.noise>0)))
xlabel('Noise intensity');
set(gca, 'LineWidth', 2, 'Box', 'off', 'FontSize', 13)


%% (3) Izhikevich Parameters
param_E = [ .02, .20, -65, 8 ]'; % E-cell
param_I = [ .02, .28, -65, 2 ]'; % I-cell

V_init=repmat([param_E(3); param_E(2)*param_E(3)],1,N);
V_init=V_init+bsxfun(@times,randn(1,N),[1; .2]);


Izhparam = [ ...
    .02, .20, -65, 8; % E
    .02, .28, -65, 2; % PV
    .02, .21, -65, 5; % SST
    .02, .22, -65, 24; % VIP
    ];

V_init=repmat([Izhparam(1,3); Izhparam(1,2)*Izhparam(1,3)],1,N);
V_init=V_init+bsxfun(@times,randn(1,N),[1; .2]);
input.V_init = V_init;

param_var = [...
    .001; % a
    .05; % b
    .50; % c
    .50]; % d

input.a =[]; input.b=[]; input.c=[]; input.d=[];
param_idx = 1;
for type=1:4 % a
    cell_idx = find( input.cell_type == type ); param_noise = 2*(rand([1,length(cell_idx)])*param_var(param_idx) - param_var(param_idx)*.5);
    input.a(cell_idx) = param_noise+Izhparam(type,param_idx); end
param_idx = 2;
for type=1:4 % b
    cell_idx = find( input.cell_type == type ); param_noise = 2*(rand([1,length(cell_idx)])*param_var(param_idx) - param_var(param_idx)*.5);
    input.b(cell_idx) = param_noise+Izhparam(type,param_idx); end
param_idx = 3;
for type=1:4 % c
    cell_idx = find( input.cell_type == type ); param_noise = 2*(rand([1,length(cell_idx)])*param_var(param_idx) - param_var(param_idx)*.5);
    input.c(cell_idx) = param_noise+Izhparam(type,param_idx); end
param_idx = 4;
for type=1:4 % d
    cell_idx = find( input.cell_type == type ); param_noise = 2*(rand([1,length(cell_idx)])*param_var(param_idx) - param_var(param_idx)*.5);
    input.d(cell_idx) = param_noise+Izhparam(type,param_idx); end

% Input scaling
input.I = input.I * F_input;
input.noise = input.noise * F_noise;

% Initializing variables
tic
tLim=input.tLim;
dt=input.dt;
t=dt:dt:tLim(2);
noise=input.noise;
numIt=numel(t);
I_inp=zeros(numIt,N);
tau_G=input.tau_G;


% parsing rest of input
a=input.a;
b=input.b;
c=input.c;
d=input.d;
V_init=input.V_init;

V=V_init;
Synapse=input.S;

%% STDP-related things
STDPflag=0; %STDPflag=input.STDP;
if STDPflag
    try
        tau_STDP=input.tau_STDP(:);
    catch
        tau_STDP=[20; 20];
        input.tau_STDP=tau_STDP;
    end
    try
        A_STDP=input.A_STDP(:);
    catch
        A_STDP=[.001; .001];
        input.A_STDP=A_STDP;
    end
    E_syn=false(size(Synapse));
    E_syn(1:sum(input.EI),1:sum(input.EI))=true;
    deltaS=zeros(numIt,2);
end
% STDP memory
if STDPflag, X=zeros(2,N); end

%% ETC
TAHpars=input.TAHpars;
maxSynVal=max(max(Synapse(E_idx,E_idx)));
input=orderfields(input);
output.input=input;
delay = input.delay;

main_input=input.I;
G=zeros(1,sum(input.nNeuron));

Smax=mean(sum(Synapse(E_idx,E_idx),2));
spiks=sparse(N,numIt);

V_AMPA = 50;
V_GABA = -90;

spikes=[];
spikes.label=cell(1,N);
spikes.label(E_idx)={'E'};
spikes.label(I_idx)={'I'};
spikes.timestamp=cell(1,N);

input.EIind=E_idx;
reverseStr=[];

I_orig=main_input;
noise_orig=noise;
interpIval=1e2; % added 1 such that simulations that start at 0 and end in multiples of 1e3 don't need one complete extra interpolation step for the final timepoint

output.S_orig=Synapse;
output.EI=input.EI;

V_list=nan(2,N,numIt);
G_list=nan(1,N,numIt);
G_list(:,:,1)=0;
if STDPflag
    X_list=nan(2,N,numIt);
    X_list(:,:,1)=0;
end

% if LFPout
try
    LFPkernel=input.LFPkernel;
    if size(LFPkernel,1) ~= N
        error('LFPkernel does not match number of neurons (should be numNeur x numLFPChan)')
    end
catch
    LFPkernel=input.EI;
end
LFPkernel=bsxfun(@rdivide,LFPkernel,sum(LFPkernel));
LFP = [];
LFP=nan(1,numIt);

% interpolate I and noise to desired resolution
% Note: this is done in steps of interpIval to save memory space
main_input=interp1(dummy_list,I_orig,t([1:interpIval]));
noise=interp1(dummy_list,noise_orig,t([1:interpIval]));
main_input(isnan(main_input))=0;
noise(isnan(noise))=0;
interpCnt=1;

if delayFlag, firSelHist=false(delay/dt+1,N); end

%% (4) Main integration loop
vidSave = False;
try,vid.close(); end
if vidSave
    vid = VideoWriter(vidname);
    vid.open();
end

draw_interval = 15;
show_win = 1001;

fig=figure(2);

clf;
set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position',  [173 159 745.6000 857]);
sample_neurons = round( linspace(N,1,20));

for n=2:numIt
    interpCnt=interpCnt+1;
    if interpCnt==(interpIval+1)
        % interpolate I and noise to desired resolution
        % Note: this is done in steps of interpIval to save memory space
        interpCnt=1;
        tVec=[1:interpIval]+n-1;
        tVec=tVec(tVec<=numel(t));
        main_input=interp1(dummy_list,I_orig,t(tVec));
        noise=interp1(dummy_list,noise_orig,t(tVec));
        main_input(isnan(main_input))=0;
        noise(isnan(noise))=0;
    end
    
    I_E= G(1,E_idx)*Synapse(:,E_idx).'.*(V_AMPA-V(1,:));
    I_I =G(1,I_idx)*Synapse(:,I_idx).'.*(V_GABA-V(1,:));
    I_tot=main_input(interpCnt,:)+I_E+I_I;
    I_tot=I_tot+noise(interpCnt,:).*randn(1,N);
    
    % membrane potential
    V(:,:)=RK4(t(n-1),V(:,:),dt,'Izh_neuron',a,b,I_tot);
    
    % synaptic conductance
    G(:,E_idx)=RK4(t(n-1),G(:,E_idx),dt,'exp_decay',tau_G(1));
    G(:,I_idx)=RK4(t(n-1),G(:,I_idx),dt,'exp_decay',tau_G(2));
    
    if STDPflag
        % STDP memory
        X(:,:)=RK4(t(n-1),X(:,:),dt,'exp_decay',tau_STDP);
        X_list(:,:,n)=X(:,:);
    end
    
    
    firSel=squeeze(V(1,:)>30);
    if delayFlag, firSelHist(mod(n,delay/dt+1)+1,:)=firSel; end
    
    if any(firSel)
        % reset membrane potential
        if numel(c)>1
            V(1,firSel)=c(firSel);
            V(2,firSel)=V(2,firSel)+d(firSel);
        else
            V(1,firSel)=c;
            V(2,firSel)=V(2,firSel)+d;
        end
        % update synaptic channels to fully open
        if delayFlag
            G(1,firSelHist(mod(n+1,delay/dt+1)+1,:))=1;
        else
            G(1,firSel)=1;
        end
        if STDPflag
            if any(firSel)% update synapses using STDP
                dsyn=TAH(t,Synapse/maxSynVal,X(:,:),A_STDP,TAHpars,firSel,Synapse>0 & E_syn);
                dsyn=dsyn*maxSynVal; % synapses are bounded/normalized by maxSynVal
                dsyn(dsyn<0)=max(dsyn(dsyn<0),-Synapse(dsyn<0)); % clip, because negative values will yield synapses with imaginary components
                Synapse=Synapse+dsyn;
                deltaS(n,1)=sum(dsyn(:));
                deltaS(n,2)=sum(abs(dsyn(:)));
                
                % update STDP memory (Only E-E and E-I interactions)
                % but only after doing STDP, cells firing at exactly the same time
                % have no influence
                %         X(:,firSel)=bsxfun(@plus,X(:,firSel),A_STDP*0+1);
                X(:,firSel)=1;
                X_list(:,:,n)=X(:,:);
            end
        end
        spiks(firSel,n)=true;
        firSel=find(firSel);
        for firIdx=1:numel(firSel)
            spikes.timestamp{firSel(firIdx)}=[spikes.timestamp{firSel(firIdx)} t(n)];
        end
    end
    LFP(:, n) = nanmean(V(1,E_idx), 2 );

    % save
    output.S=Synapse;
    output.spiks=spiks(:,1:n);
    output.t=t(1:n);
    output.simulationTime=toc;
    if STDPflag, output.deltaS=deltaS(1:n,:); output.X=X_list; end
    
    G_list(:,E_idx,n)=G(:,E_idx);
    G_list(:,I_idx,n)=G(:,I_idx);
    V_list(:,:,n)=V;
    I_inp(n,:)=I_tot;
    output.G=G_list(:,:,1:n);
    output.I_tot=I_inp(1:n,:);
    output.V=V_list(:,:,1:n);
    
    %% Visualization
    if mod(n,draw_interval)==1
        
        %% All cells' V(t)
        subplot(2,2,1); hold off;
        tickinterval = 250;
        imagesc( [1:size(V_list,3)]*input.dt, 1:size(V_list,2), ...
            squeeze(V_list(1,:,:) )  );
        hold on; plot( xlim, [0 0]+nNeuron(1)+1, 'r-', 'LineWidth', 2 ); axis xy;
        caxis([-80 -30]);
        xlabel('Time (ms)');
        set(gca, 'YTick', [ nNeuron(1)*.5, nNeuron(1)+nNeuron(2)*.5 ], ...
            'YTickLabel', {'E', 'I'} );
        set(gca, 'FontSize', 14, 'LineWidth', 2, 'Box', 'off') ;
        xlim([n-show_win n+10]*input.dt);
        title('Membrane voltage, V(t)');
        xlims=xlim;
        set(gca, 'XTick', round([xlims(1)+1 mean(xlims) xlims(2)-1]))
        
        %% Sample spike trajectory
        subplot(2,2,3); hold off;
        sample_volt = squeeze( V_list(1,sample_neurons ,:));
        ys = length(sample_neurons)*25 ;
        erps=plot_multi_ch( sample_volt, [1:size(V_list,3)]*input.dt, ys);
        xlim([n-show_win n+10]*input.dt);
        xlabel('Time (ms)');
        set(gca, 'FontSize', 14, 'LineWidth', 2, 'Box', 'off') ;
        ylabel('Sample neurons V(t)')
        ylim([-1.1 1.1]*ys );
        xlims=xlim;
        set(gca, 'XTick', round([xlims(1)+1 mean(xlims) xlims(2)-1]))
        hold on;
        
        %% Mean field (LFP)
        subplot(2,2,2); hold off;
        plot((1:size(V_list,3))*input.dt, LFP, 'k' )
        
        xlims=[n-show_win n+10]*input.dt;
        hold on; % Gamma indicator
        sshift= -100;
        plot([xlims(end)-25 xlims(end)]+sshift, [0 0]-80, 'k-' ,'LineWidth', 2 );
        text( xlims(end)+sshift*.85, -82, '40 Hz' );
        xlim(xlims);
        ylim([ -85, -60]);
        set(gca, 'YTick', []);
        title('LFP, sum of V(t)_E');
        xlabel('Time (ms)');
        ylabel('Voltage');
        set(gca, 'FontSize', 14, 'LineWidth', 2, 'Box', 'off') ;
        xlims=xlim;
        set(gca, 'XTick', round([xlims(1)+1 mean(xlims) xlims(2)-1]))
        hold on;
        
        try hann_win;
        catch;
            hann_win=hanning(length(n-show_win:n));
        end
        
        try
            subplot(2,2,4); hold off;
            dat_for_fft = LFP(n-show_win:n) .* hann_win';
            [X,f]=positiveFFT(dat_for_fft, 1000/dt, 0);
            rz = 3; smt=10;
            plot(imresize(f, [1, length(f)*rz], m), ...
                smooth( imresize(abs(X), [1, length(X)*rz],m), smt).^2, 'k-', 'LineWidth', 2 );
            
            xlim([7 100]);
            ylim([0 .6]);
            set(gca, 'YTick', [0 .6]);
            xlabel('Freq (Hz)'); ylabel('|FFT|');
            title('LFP power spectrum');
            set(gca, 'FontSize', 14, 'LineWidth', 2, 'Box', 'off') ;
        end
        colormap gray;
        drawnow;
        
        if vidSave
            frr = getframe( gcf );
            vid.writeVideo(frr.cdata);
        end
    end
    msg=sprintf(['Iteration %d/%d\n'], [n numIt]);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
try; vid.close(); end

%% (6) Collect simulation output
output.t=t;
output.S=Synapse;
output.spiks=spiks;
output.simulationTime=toc;
if STDPflag, output.deltaS=deltaS;
    output.X=X_list;end
output.LFP=LFP;

G_list(:,E_idx,n)=G(:,E_idx);
G_list(:,I_idx,n)=G(:,I_idx);
V_list(:,:,n)=V;
I_inp(n,:)=I_tot;
output.G=G_list;
output.I_tot=I_inp;
output.V=V_list;



