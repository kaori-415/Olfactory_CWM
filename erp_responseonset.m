function ERPres = erp_responseonset(pflag, rmth)

% ERP analysis for response phase
% baseline correction : use a common baseline (pre 1st. response)
% output ERPresponse structure and plot figures
%  input: pflag = plot flag. if pflag == 1, output and save figs.
%                 if pflag == 0, no plotting, just output ERP structure
%         rmth = removing threshold for EOG, ex. rmth=80, remove trials over
%                80uV EOG signals
%                default: 80

conditions = {'control', 'odorant'};

pretime = 500;
limittime = 1000;
baseline = 200;

colornum = 4;

% rmth = 80;

Behavior = load('behavior/result.mat');
load('EEG_source/PreEEG.mat', 'Response');
load(fullfile('EEG_source', 'PreEEG.mat'), 'Ch_order');

% check plot flag
if nargin < 1
    pflag = 1;
end

if nargin < 2
    rmth = 80;
end

% epoch by colors
for ncond = 1:size(Response,2)
    
    % get color idx
    tmpresponse = Behavior.Response(ncond,:);
    coloridx = arrayfun(@(x) x.answerColorIdx(:,1), tmpresponse, 'UniformOutput', false);
    coloridx = horzcat(coloridx{:});
    
    trials = [Response(ncond).order.epoch];
    trials = reshape(trials,numel(trials)/2, 2)';
    
    % eog artifact to reject trials
    EOG = cellfun(@(x) x(:,11:12), trials, 'UniformOutput', false);
    rmtrial_cell= cellfun(@(x) find(x>rmth), EOG, 'UniformOutput', false); % index containing artifacts
    
    % get condition labels
    condlabel = Behavior.Conditions.labels{ncond};
    
    for ncolor = 1:colornum
        
        ERP = zeros(pretime+limittime+1, 12, 2);
        ntrial = zeros(1,2); 
        
        for norder = 1:2
            
            % get ERP trials
            % from one color group and less containing EOG artifacts
            erptrials = trials(norder, coloridx(norder,:) == ncolor &  cellfun(@isempty, rmtrial_cell(norder,:)));
            
            % reject too much shor trials
            trlen = cellfun(@(x) size(x,1), erptrials);
            erptrials(trlen<[pretime+limittime+1]) = [];
            
            erptrials_fix = cellfun(@(x) x(1:pretime+limittime+1, :), erptrials, 'UniformOutput', false);
            erpmat = cat(3, erptrials_fix{:});
            
            erp = squeeze(mean(erpmat,3)); % average
            
            % baseline correction
            d_base = pretime-baseline;
            base = mean(erp(d_base+1 : pretime+1, :));

            
            erp = erp - repmat(base, size(erp,1), 1);
            
            ERP(:,:,norder) = erp;
            ntrial(:,norder) = size(erpmat,3);
        end
        clear base
        
        ERPres.(condlabel).color(ncolor).erp = ERP;
        ERPres.(condlabel).color(ncolor).ntrial = ntrial; 
    end
end

% plot ERPs
switch pflag
    case 1
        
        linergb = [1, 0.51, 0; 0, 0.81, 0.28; 0, 0.78, 1; 1, 0.23, 1];
        
        for ncolor = 1:colornum
            h(ncolor) = figure('Position', [1,1,800,1200], 'Visible', 'off');
            
            for nplot = 1:8 %for channels
                
                for ncond = 1:length(conditions)
                    % set linestyle
                    condlabel = Behavior.Conditions.labels{ncond};
                    switch condlabel
                        case conditions{1}
                            style = '--';
                        case conditions{2}
                            style = '-';
                    end
                    
                    plotdata = ERPres.(condlabel).color(ncolor).erp;
                    
                    ax1 = subplot(8, 2, nplot*2-1);
                    hold(ax1,'on')
                    plot(-pretime:limittime, plotdata(:,nplot,1), 'LineStyle', style, 'Color', linergb(ncolor,:));
                    set(gca, 'YDir', 'reverse')
                    title([condlabel, ' ', Ch_order{nplot}, ' 1st']);
                    xlabel('Time');
                    ylabel('Amplitude');
                    axis tight
                    hold(ax1,'off')
                    
                    ax2 = subplot(8, 2, nplot*2);
                    hold(ax2, 'on')
                    plot(-pretime:limittime, plotdata(:,nplot,2), 'LineStyle', style, 'Color', linergb(ncolor,:));
                    set(gca, 'YDir', 'reverse')
                    title([condlabel, ' ', Ch_order{nplot}, ' 2nd']);
                    xlabel('Time');
                    ylabel('Amplitude');
                    axis tight
                    hold(ax1,'off')
                    
                end
                
            end
            
            legend(Behavior.Conditions.labels{:}, 'Location','southwest','Orientation','horizontal')
            legend('boxoff')
        end
        
        % save figs
        if ~exist('ERP', 'dir')
            mkdir('ERP')
        end
        
        for ncolor = 1:colornum
            set(h(ncolor), 'PaperPositionMode', 'auto')
            figname = fullfile('ERP', ['ERP_responseonset_color', num2str(ncolor)]);
            print(h(ncolor), figname,'-dpng', '-r0')
        end
        
end