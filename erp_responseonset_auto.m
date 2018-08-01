function erp_responseonset_auto(ch, diridx)

% grand average of erp waveforms
% ch = channel, diridx = subjidx

% ch = 'Cz';
% diridx = 1:20; 

if nargin < 2
    diridx = [1:4, 6:10, 12:16, 18:20]; % remove 5, 11, 17
end

root =cd;
addpath(root);

resultsRoot = 'Data';
DirHeadStr = 's';
condition = {'control', 'odorant'};

colornum = 4;

Ch_order = {'Fz', 'Cz', 'Pz', 'Oz', 'F3', 'F4', 'P3', 'P4','A1', 'A2', 'EOG1', 'EOG2', 'trig'};

pretime = 500;
posttime = 1000;

% get data
cd(resultsRoot);

% make dirnames
diridxCell = num2cell(diridx);
diridxStr = cellfun(@num2str, diridxCell, 'UniformOutput', false);
zeropadIdx = find(cellfun(@length, diridxStr)<2);
diridxStr(zeropadIdx) = strcat('0', diridxStr(zeropadIdx)); % padding zero as character like '01'
dirnames = strcat(DirHeadStr, diridxStr); % make dirnames including results
dirnamesLen = length(dirnames);


for n = 1:dirnamesLen
    cd(dirnames{n});
    RESPONSE(n,:) = orderfields(erp_responseonset(0));
    
    cd('../');
    
end

close all

linergb = [1, 0.51, 0; 0, 0.81, 0.28; 0, 0.78, 1; 1, 0.23, 1];

for ncolor = 1:colornum

    h = figure;
    
    for ncond = 1:length(condition)
        GrandAverage = cell(1, dirnamesLen);
        [GrandAverage{:}] = RESPONSE.(condition{ncond});
        GrandAverage_col = cellfun(@(x) x.color(ncolor).erp, GrandAverage, 'UniformOutput', false);
        
        %for nplot = 1:8
            
            nplot = find(strcmp(ch, Ch_order));
            condlabel = condition{ncond};
            
            subplot(1,2,ncond)
            hold on; 
            
            maxY = 25; minY = -15; 
            rectangle('Position', [300, minY, 200, maxY-minY], 'FaceColor', [0.7,0.7,0.7], 'EdgeColor', [0.7,0.7,0.7]); % latency range
            
            % plot individuals
            for nsubj = 1:dirnamesLen
                plot(-pretime:posttime, GrandAverage_col{nsubj}(:,nplot,1),...
                   'Color', [0.5,0.5,0.5]);
%                   plot(-pretime:posttime, GrandAverage_col{nsubj}(:,nplot,2),...
%                      'Color', [0.5,0.5,0.5]);
            end
            
            % plot grandaverage
            allerp = cellfun(@(x) x(:,nplot,1), GrandAverage_col, 'UniformOutput', false);
            % allerp = cellfun(@(x) x(:,nplot,2), GrandAverage_col, 'UniformOutput', false);
            granderp = mean(horzcat(allerp{:}),2);
            
            plot(-pretime:posttime, granderp, 'Color', linergb(ncolor, :), 'LineWidth', 2);
            title([condlabel, ' ', Ch_order{nplot}, ' 1st'], 'FontSize', 14);
            % title([condlabel, ' ', Ch_order{nplot}, ' 2nd'], 'FontSize', 14);
            xlabel('Time');
            ylabel('Amplitude');
            % axis tight
            xlim([-200, 1000]);
            ylim([minY, maxY]);
            set(gca, 'YDir', 'reverse')

            clear granderp

            fontname_helvetica; % font setting
        
    end
    
    savedir = fullfile('DAll', '002_ERP_fig', '002_response_GA');
    figname = ['ERP_color', num2str(ncolor), '_', Ch_order{nplot}, '_1st'];
    saveas(h, fullfile(savedir,figname), 'png');
end



cd(root) 
close all
            