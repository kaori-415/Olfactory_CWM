function FourColorExp_hue

% Multiple color presentation for color working memory experiment 
% This script refers to codes developed by Brady, T. F., Konkle, T.F., Gill, J., Oliva, A. and Alvarez, G.A. (2013). 
%  "Visual long-term memory has the same limit on fidelity as visual working memory". Psychological Science, 24(6), 981-990. 
%  https://bradylab.ucsd.edu/stimuli.html

%------------------
% config
%------------------

%Screen('Preference', 'SkipSyncTests', 1);

AssertOpenGL;
GetSecs;
WaitSecs(1);
rng('shuffle')


%------------------
% User definition
%------------------

dpath = uigetdir('../003_results');
cd(dpath)

% standard color from Zhang and Luck 2008, Nature
Lab = [70,20,38];
radius = 60;

% number of trials for each color
ntrial = 12*2; % 6‚Ì”{”
nsession = 2; 

% preferences for colors in an experiment
stepdegrees = 0:90:360;
groupdegrees = -20:10:20;

% rectangle size
nX = 100;
nY = 100;

% Make a base Rect
baseRect = [0, 0, nX, nY];
nrectAns = 2; 
nrectFill = 2;  
numsquares = nrectAns + nrectFill; 

% rectangle location
%pos = 0.25;
pos = 200;

% inter stimulus interval
dispdur = 0.2;
waitdur = 1.5;
fixsdur = 0.5; % duration of the fixation cross
iti = 0.2:0.1:0.5;  % inter stimulus interval - fixsdur

%----------------
% make colors
%----------------

% make cforms to convert colors
lab2lchform = makecform('lab2lch');
lch2labform = makecform('lch2lab');
lab2rgbform = makecform('lab2srgb');

% convert standard lab color to cielch
LchCenter = applycform(Lab, lab2lchform);
LchStandard = [LchCenter(1), LchCenter(2)+radius, LchCenter(3)];

% make gray color with the same luminance level
Lchmono = LchStandard; 
Lchmono(2:3) = 0;
gray = applycform(Lchmono, lab2rgbform);

% make colors
for ngroup = 1:length(stepdegrees)-1
    
    for ndeg = 1:length(groupdegrees)
        
        tmpLch = LchStandard; 
        tmpLch(3) = tmpLch(3) + stepdegrees(ngroup) + groupdegrees(ndeg); 
        Lchs(ndeg,:) = tmpLch;
        
        Labs = applycform(tmpLch, lch2labform);
        RGBs(ndeg, :) = applycform(Labs, lab2rgbform);
        
    end
    
    Lchcell{ngroup} = Lchs;
    RGBcell{ngroup} = RGBs; 
    
end

clear RGBs Labs

% randomize color selection
ColorIdx = repmat([1:length(groupdegrees)],[length(stepdegrees)-1,ntrial,nsession]); % ColorIdx(nrow, ntrial) == Nc, -> rgb = RGBcell{nrow}(Nc,:)
FillColorIdx = ColorIdx; 

% randomize ITI
ITI = repmat(iti', [ceil(size(ColorIdx,2)/length(iti)),1]);
ITIcell = repmat({ITI}, 1, nsession); 
rnITIcell = cellfun(@Shuffle, ITIcell, 'UniformOutput', false);
rnITI = cat(3, rnITIcell{:}); 

% trig rect 
trig_baseRect = [0, 0, 50, 50];


%---------------------
% Start Experiment
%---------------------
PsychDefaultSetup(2);

%try
% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, (white+black)/2);
%[window, windowRect] = Screen('OpenWindow',screenNumber,[],[10 20 1000 600]); % debug mode

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);
% Sync us and get a time stamp
vbl = Screen('Flip', window);
waitframes = 1;

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% set locations for rectangles
squareXpos = [xCenter, xCenter, xCenter+pos, xCenter-pos];
squareYpos = [yCenter+pos, yCenter-pos, yCenter, yCenter];

% trigger rectangle setting
trigRect = CenterRectOnPointd(trig_baseRect, 50, screenYpixels-50);

% gamma correction
load mygamma.mat % mygamma.mat is a gamma table for display gamma correction, please see 'Screen LoadNormalizedGammaTable?'
BackupCluts(screenNumber);
Screen('LoadNormalizedGammaTable', screenNumber, mygamma);

 for m = 1:nsession
    
    % draw an instruction
    Screen('TextSize', window, 50);
    DrawFormattedText(window, 'Please wait until start', 'center', 'center', black);
    % draw a trigger rectangle
    Screen('FillRect', window, black, trigRect); 
    
    vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    
    KbWait;
    %GetClicks;
    
    Screen('TextSize', window, 50);
    DrawFormattedText(window, 'Start', 'center', 'center', black);
    % draw a trigger rectangle
    Screen('FillRect', window, black, trigRect); 
    vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    
    WaitSecs(1);
    
    % draw a fixation cross
    Screen('TextSize', window, 50);
    DrawFormattedText(window, '+', 'center', 'center', black);
    % draw a trigger rectangle
    Screen('FillRect', window, black, trigRect); 
    
    % Flip the screen
    vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    
    WaitSecs(fixsdur);
    
    % draw an empty screen
    Screen('TextSize', window, 100);
    DrawFormattedText(window, '', 'center', 'center', black);
    % draw a trigger rectangle
    Screen('FillRect', window, black, trigRect); 
    
    % Flip the screen
    vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    WaitSecs(fixsdur*2);
    
    % make a randomized color idx for stepdegrees
    color_comb = combnk([1:length(stepdegrees)-1],2)';
    TargetC = repmat(color_comb,1, size(ColorIdx,2)/size(color_comb,2));
    TargetC_idx = randperm(size(TargetC,2)); 
    
    % make a randomized color for ColorIdx cols
    for ndeg = 1:length(stepdegrees)-1
        TargetC_deg(ndeg, :) = randperm(ntrial*(length(groupdegrees)));
    end
    
    % make a randomized filler colors for ColorIdx cols
    for ndeg = 1:length(stepdegrees)-1
        FillerC_deg(ndeg, :) = randperm(ntrial*(length(groupdegrees)));
    end
    
    % color counters 
    counters = ones(1, nrectAns + nrectFill); 
    
    % roop for trials within a session
    for n = 1:size(TargetC,2)
        
        HideCursor;
       
        % get colors and make rectangles
        targetidx = Shuffle(TargetC(:,TargetC_idx(n))); % rectangle indexes for target rects
        tmpidx = 1:numsquares;
        fillidx = tmpidx(~ismember(tmpidx, targetidx)); % filler rectangle indexes
        rectidx = [targetidx', fillidx]; % rectidx make random locations of rects
        
        Allrects = nan(4, numsquares);
        rgb = zeros(numsquares,3);
        lch = zeros(numsquares,3);
        rectlocation = [squareXpos(rectidx); squareYpos(rectidx)];
        
        for nr = 1:length(rectidx)
           % check counters
           count = counters(rectidx(nr));
           
           % get color_idx 
           if nr < 3
           color_idx = TargetC_deg(rectidx(nr), count);
           degidx = ColorIdx((rectidx(nr)), color_idx, m);
           else 
               color_idx = FillerC_deg(rectidx(nr),count);
               degidx = FillColorIdx((rectidx(nr)), color_idx, m);
           end
           
           % add counters
           counters(rectidx(nr)) = counters(rectidx(nr))+1;
           
           % make rgbs and lchs 
           rgb(nr,:) = RGBcell{rectidx(nr)}(degidx,:);
           lch(nr,:) = Lchcell{rectidx(nr)}(degidx,:);
           
           % make rectangles
           Allrects(:,nr) = CenterRectOnPointd(baseRect, rectlocation(1,rectidx(nr)), rectlocation(2,rectidx(nr)));
           % draw rectangles
           Screen('FillRect', window, rgb(nr,:), Allrects(:,nr));
        end
        
        % draw a fixation cross
        Screen('TextSize', window, 50);
        DrawFormattedText(window, '+', 'center', 'center', black);
        
        % draw a trigger rectangle
        Screen('FillRect', window, white, trigRect); 
        
        % Flip the screen
        vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        WaitSecs(dispdur);
        
        % mask the rects with gray
        for nr2 = 1:numsquares
            Screen('FillRect', window, gray, Allrects(:,nr2));
        end
        
        % draw a fixation cross
        Screen('TextSize', window, 50);
        DrawFormattedText(window, '+', 'center', 'center', black);
        
        % draw a trigger rectangle
        Screen('FillRect', window, black, trigRect); 
        
        % Flip the screen
        vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        WaitSecs(waitdur);
        

        
        for nc = 1:nrectAns
            ncount = nc; 

            % select the target rectangle            
            nrect = rectidx(ncount); 
            
            % set mouse location at the center
            SetMouse(rectlocation(1,nrect),rectlocation(2,nrect), window);
            ShowCursor; 
            
            % colorwheel settings
            colorWheel.radius = nX;
            colorWheel.rect = CenterRectOnPointd([0,0, colorWheel.radius*2, colorWheel.radius*2], ...
                                                 rectlocation(1,nrect),rectlocation(2,nrect));
            
            % draw a frame on an answer_target rectangle
            Screen('FrameOval', window, black, colorWheel.rect);
                    
            % mask the rects with gray
            for nr2 = 1:numsquares
                Screen('FillRect', window, gray, Allrects(:,nr2));
            end
            
            % draw a fixation cross
            Screen('TextSize', window, 50);
            DrawFormattedText(window, '+', 'center', 'center', black);
            % draw a trigger rectangle 
            Screen('FillRect', window, white, trigRect); 
            
            % Flip the screen
            vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                    
            % location of target rectangle
            NewRectlocate = [rectlocation(1,nrect)-(nX/2),rectlocation(2,nrect)-(nY/2), ...
                rectlocation(1,nrect)+(nX/2),  rectlocation(2,nrect)+(nY/2)];
            
            % make a random rgb color for the target rectangle
            d_lch = randperm(359,1);
            d_LchS = LchStandard;
            d_LchS(3) = d_lch; 
            tmprgb = applycform((applycform(d_LchS, lch2labform) ), lab2rgbform);
            
            % location of other rectangles
            nstore_rect = find([1:numsquares]~=nrect);
            nstoreRectlocate = [rectlocation(1,nstore_rect)'-(nX/2), rectlocation(2,nstore_rect)'-(nY/2), ...
                rectlocation(1,nstore_rect)'+(nX/2), rectlocation(2,nstore_rect)'+(nY/2)];
            % mask an alternative rectangle as gray 
            for nstore = 1:numsquares -1
                nstoreRect{nstore} = NewMakeTexture(gray, nX, nY, window);  
            end
            
            oldAngle = 0;
            buttons = [0,0,0];
            
            tic; 

            % animation untile get clicks
            while ~ismember(1,buttons)
                
                % show an oval
                Screen('FrameOval', window, black, colorWheel.rect);
                % show a fixation cross
                DrawFormattedText(window, '+', 'center', 'center', black);
                
                % get cursor
                [curX2,curY2, buttons] = GetMouse(window);
                    
                curAngle = GetPolarCoordinates(curX2,curY2,rectlocation(1,nrect),rectlocation(2,nrect)); 
                
                % change color
                if (round(curAngle) ~= round(oldAngle)) && round(curAngle) ~= 0
                    [newrgb,newlch] = RotateColors(round(curAngle), d_LchS, lch2labform,lab2rgbform);
                    tmprgb = newrgb;
                end
                
                NewRect = NewMakeTexture(tmprgb, nX, nY, window);
                
                % Draw New Rectangle
                Screen('DrawTexture', window, NewRect, [], NewRectlocate);
                
                % Draw an alternative Rectangle 
                for nstore = 1:numsquares -1
                    Screen('DrawTexture', window, nstoreRect{nstore}, [], nstoreRectlocate(nstore,:));
                end
                
                % draw a trigger rectangle
                Screen('FillRect', window, white, trigRect); 
                
                % Flip the screen
                vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                oldAngle = curAngle;
                
            end
            
            % RECORD response
            t = toc;
            
            % insert trig
            Screen('FillRect', window, black, trigRect); 
            % keep other shapes
            Screen('DrawTexture', window, NewRect, [], NewRectlocate);
            DrawFormattedText(window, '+', 'center', 'center', black);
            for nstore = 1:numsquares -1
                Screen('DrawTexture', window, nstoreRect{nstore}, [], nstoreRectlocate(nstore,:));
            end
            % Flip the screen
            vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            
            % save responses
            Response(m,n).rectID(nc) = nrect;
            Response(m,n).rectLocation(nc,:) = rectlocation(:,nrect); % [xLocation,yLocation]
            Response(m,n).answerColorIdx(nc,:) = [nrect,degidx];%answerColorIdx(1) = colorgroup, answerColorIdx(2) = groupdegree
            Response(m,n).answerRGB(nc,:) = rgb(nc,:);
            Response(m,n).responseRGB(nc,:) = newrgb;
            Response(m,n).answerLch(nc,:) = lch(nc,:);
            Response(m,n).responseLch(nc,:) = newlch;
            Response(m,n).responseTime(nc,:) = t;
            
            
            %clear newrgb nstore_rect_rgb
            assignin('base', 'Response', Response);

            WaitSecs(0.1);
            
        end
        % save randomized idx 
        Idx(m).TargetC = TargetC; 
        Idx(m).TargetC_idx = TargetC_idx;
        Idx(m).TargetC_deg = TargetC_deg; 

        % draw a fixation cross
        Screen('TextSize', window, 50);
        DrawFormattedText(window, '+', 'center', 'center', black);
        
        % insert trig
        Screen('FillRect', window, black, trigRect); 
        
        % Flip the screen
        vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        
        WaitSecs(fixsdur);
        
        % draw an empty screen
        Screen('TextSize', window, 100);
        DrawFormattedText(window, '', 'center', 'center', black);
        % insert trig
        Screen('FillRect', window, black, trigRect); 
        
        % Flip the screen
        vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        WaitSecs(rnITI(n,:,m));
    end
end

Screen('TextSize', window, 100);
DrawFormattedText(window, 'Finish', 'center', 'center', black);
% draw a trigger rectangle
Screen('FillRect', window, black, trigRect); 

% Flip the screen
vbl  = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
WaitSecs(2);

% Clear the screen
sca;
close all;
ShowCursor; 


% save response
if ~exist('behavior', 'dir')
    mkdir('behavior')
end
save behavior/result.mat Response Idx RGBcell Lchcell
        
end

% ----------------------------------------------------------
function [angle, radius] = GetPolarCoordinates(h,v,centerH,centerV)
% get polar coordinates
hdist   = h-centerH;
vdist   = v-centerV;
radius     = sqrt(hdist.*hdist + vdist.*vdist)+eps;

% determine angle using cosine (hyp will never be zero)
angle = acos(hdist./radius)./pi*180;

% correct angle depending on quadrant
angle(hdist == 0 & vdist > 0) = 90;
angle(hdist == 0 & vdist < 0) = 270;
angle(vdist == 0 & hdist > 0) = 0;
angle(vdist == 0 & hdist < 0) = 180;
angle(hdist < 0 & vdist < 0)=360-angle(hdist < 0 & vdist < 0);
angle(hdist > 0 & vdist < 0)=360-angle(hdist > 0 & vdist < 0);
end

% ----------------------------------------------------------
function [newrgb,lch] = RotateColors(Angle, oldlchs, lch2labform, lab2rgbform)

lch = oldlchs;
lch(3) = lch(3) + Angle;
lab = applycform(lch, lch2labform);
newrgb = applycform(lab, lab2rgbform);

end

function newrect = NewMakeTexture(rgb, nX, nY, window)

theRGBCalFormat = rgb'*ones(1,nX*nY);
theRGBImage = CalFormatToImage(theRGBCalFormat,nX,nY);
newrect = Screen('MakeTexture', window, theRGBImage);

end


