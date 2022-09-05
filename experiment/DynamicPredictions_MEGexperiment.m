% script to run MEG experiment which was used for scientific article:
% Predictive neural representations of naturalistic dynamic input
% Ingmar E.J. de Vries and Moritz F. Wurm (2022) 
% experiment can run as 1) behaviour only, 2) behaviour + eyetracker, 3) behaviour + eyetracker + MEG
% note that you need to have Psychophysics Toolbox Version 3 (PTB-3) downloaded and installed in order to run this experiment

% Clear the workspace
close all;
clear all;
clc;
sca;

% For debugging purposes this makes experiment screen transparent.
debug = false;
if debug  
    PsychDebugWindowConfiguration([],.4);
    Screen('Preference', 'SkipSyncTests', 1);%don't do screen sync + size tests
end
Screen('Preference', 'SkipSyncTests', 0);%do screen sync + size tests

%----------------------------------------------------------------------
%                        INITIALIZE EXPERIMENT
%----------------------------------------------------------------------

% SOME GENERAL SETTINGS
rootdir = 'type_rootdir_here\';
addpath(genpath(rootdir));

% turn eyetracker and MEG on/off
% if MEG is selected, eyetracker will be automatically included as well
trackeye = false;
MEG = true;
if MEG
    trackeye = true;
end
if trackeye
    Screen('Preference', 'SkipSyncTests', 0);
end

% Number of repetitions of unique trials/sequences per block (out of 14 unique sequences)
trialsPerCombination = 5;

% determine amount of practice trials.
trialAmountPractice = 60;

% seed random number generator
rng('default');
rng('shuffle'); 

% select subject number and whether this run is a practice run or not, then create file name
experiment_name = 'ActionPrediction';
subject_num = input('Please indicate subject number: ');
practice = input('Is this a practice run? [1/0]');
if practice
    experiment_name = [experiment_name '_practice'];
    block = 0;
else
    block = input('Please indicate number of real block: ');
end

if subject_num < 10
    filename = [experiment_name '_subject0' num2str(subject_num) '_block' num2str(block)];
else
    filename = [experiment_name '_subject' num2str(subject_num) '_block' num2str(block)];
end

savedir = [rootdir 'data\'];
addpath([rootdir 'experiment']);

if ~practice && exist([savedir filename '.mat'],'file')
    warning('Filename already exists, please make sure you have indicated the correct subject and block number!');
    return
end

% Monitor refresh rate, in MEG lab 120 Hz, in behavioural or eyetracking-only lab 60 Hz
% refrate = X; % currently we determine the refrate below with the PTB function: Screen('GetFlipInterval',window);, but sometimes when connected via
% remote desktop this doesn't work, so here you can force a certain refresh rate instad
display.dist = 58;% distance from screen (cm)
display.width = 41; %width of screen (cm)
display.resolution = 1280; % number of pixels of display in horizontal direction

% PSYCHTOOLBOX
PsychDefaultSetup(2);
screens = Screen('Screens');
screenNumber = 1;% max(screens); This determines which screen will be used for running the experiment. Depending on your monitor setup you might want to change this

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'UseVirtualFramebuffer');

% define colors
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white / 2;
bg = grey*0.8;
colors = [27 158 119;217 95 2;117 112 179;231,41,138;102,166,30]/255;%RGB values of colors that are colorblind and grayscale friendly
colorPlaceHolder = [0 0 0];% will be placed instead of video frames during inter-trial-interval

% Open the experiment window
[window, windowRect] = PsychImaging('OpenWindow',  screenNumber, bg, [], 32, 2);
refrate = 1/Screen('GetFlipInterval', window);
refrate = round(refrate);

% Set PTB to top priority so no other running processes on this PC interfer
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);

% Flip to clear
Screen('Flip', window);

% Get screen info
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[screenXpixels, screenYpixels] = Screen('WindowSize', window);%get screen size
[xCenter, yCenter] = RectCenter(windowRect);% get screen center

% For real experiment hide cursor
if ~debug
    HideCursor;
end

% STIMULI & CONDITION MATRIX
% define stimulus size, screen locations, and load videos
Screen('TextSize', window, 28);
fixCrossDimPix = angle2pix(display,.3);
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
lineWidthPix = 2;

%sizes in pixels
stimHeight= 400;
stimWidth = 376;
rect4stim = [0 0 stimWidth stimHeight];
%stimPos 1 = left, 2 = right, 3 = center
stimPos = CenterRectOnPointd(rect4stim, xCenter, yCenter);

%load all videos here, and transform into textures
stimdir = [rootdir 'experiment\stimuli'];
videos = dir(fullfile(stimdir, '*vid.mp4'));
stimnum = length(videos);
framenum = 5*refrate;% 5-sec videos
stimuli = cell(stimnum,framenum);% cell with size: [sequences x frames]

DrawFormattedText(window,...
    'Loading videos...',...
    'center', yCenter/3, black);
Screen('Flip', window);

for istim = 1:stimnum
    
    %load videos and change frames to textures
    videoName = fullfile(stimdir,videos(istim).name);
    videoHeader = VideoReader(videoName);
    video = read(videoHeader);
    video = im2single(video);
    
    %resample video to screen refresh rate using the 'nearest' method for refrate used in MEG scanner,
    %i.e. at new sample points it picks RGB values from nearest frame for each pixel
    fsOld = videoHeader.FrameRate;
    fsNew = refrate;
    tOld = 0:1/fsOld:5-1/fsOld;
    tNew = 0:1/fsNew:5-1/fsNew;
    
    video = permute(video,[4 1 2 3]);
    videoNew = interp1(tOld,video,tNew,'nearest');
    clear video% save WM storage
    videoNew = im2uint8(videoNew);
    videoNew = permute(videoNew,[2 3 4 1]);
    
    for iframe = 1:framenum
        stimuli{istim,iframe} = Screen('MakeTexture', window, squeeze(videoNew(:,:,:,iframe)));
    end
    clear videoNew
end

% All variables
% video sequence
sequence = 1:stimnum;

% catchtrial: 1 = catch trial, 2 = normal trial
if practice
    catchprob = 0.3;%for practice block set at 30% so subjects get enough practice with task
else
    catchprob = 0.2;%for normal block set at 20%, but because of division by 3 for 3 types of catch trials below, it ends up being 22%
end
% Make the matrix which will determine our condition combinations
% Row 1 = video sequence, row 2 = catch or no catch trial
condMatrixBase = sequence;

% Duplicate the condition matrix to get the full number of trials
condMatrix = repmat(condMatrixBase, 1, trialsPerCombination);
condMatrix(2,:) = 0;% no catch trials yet

%new total trial amount
trialAmountTotal = size(condMatrix,2);

% Randomise the conditions
shuffler = Shuffle(1:trialAmountTotal);
condMatrixShuffled = condMatrix(:, shuffler);

% Add catch or no catch trial variable
% of the total amount of catch trials in this block: 1/3 = kinematics task, 1/3 = semantics task, 1/3 = fixation cross task
numcatch = ceil((size(condMatrix,2)*catchprob)/3);

% catchtrialpool contains manually created catch trials (i.e., manually created because the statement regarding body motion or ballet figure needs to be possible to answer)
[CatchTrialsKinematics, CatchTrialsSemantics, CatchTrialsKeep] = CreateCatchTrials();
% or if previously created, load here:
% load([rootdir '\experiment\CatchTrialPool'],'CatchTrialsKinematics','CatchTrialsSemantics','CatchTrialsKeep');

% randomly pick right amount of catch trials from each category, and add the fixation cross catch trials of which the timing is set truly randomly:
catchTrialPool = [CatchTrialsKinematics(:,randperm(stimnum,numcatch)) CatchTrialsSemantics(:,randperm(stimnum,numcatch)) CatchTrialsKeep(:,randperm(stimnum,numcatch))];
for i=numcatch*3+1:numcatch*4
    catchTrialPool{1,i} = randperm(stimnum,1);
    catchTrialPool{2,i} = 0.4 + 4.4 * rand(1,1);
    catchTrialPool{3,i} = 0;
    catchTrialPool{4,i} = 1;% correct response is always left (red button) for fixation cross task
end
numcatch = size(catchTrialPool,2);

%randomize order of catch trials within a block
catchTrialPool = catchTrialPool(:,randperm(numcatch,numcatch));

% Make sure the same sequence never directly repeats (i.e. at least 1 other sequence in between)
for i=1:trialAmountTotal-2
    if condMatrixShuffled(1,i) == condMatrixShuffled(1,i+1)
        
        % if sequence i == sequence i+1, swap sequence i+1 with sequence i+2
        temp = condMatrixShuffled(:,i+1);
        condMatrixShuffled(:,i+1) = condMatrixShuffled(:,i+2);
        condMatrixShuffled(:,i+2) = temp;
    end
end

% Insert catch trials with a uniform distribution with the constraint:
% not too many normal trials in between two catch trials
% accomplished by dividing (trialAmountTotal + numcatch) by numcatch, and having 1 catch trial in each of these equal chunks of trials 
for icatch = 1:numcatch
    
    spot2insert = ceil(ceil((icatch-1)*((trialAmountTotal+numcatch)/numcatch)) + floor((trialAmountTotal+numcatch)/numcatch)*rand(1,1));
    
    condMatrixShuffled(:,spot2insert+1:trialAmountTotal+icatch) = condMatrixShuffled(:,spot2insert:trialAmountTotal+(icatch-1));
    condMatrixShuffled(1,spot2insert) = catchTrialPool{1,icatch};
    condMatrixShuffled(2,spot2insert) = icatch;
end
trialAmountTotal = trialAmountTotal + numcatch;

% MAKE AN OUTPUT MATRIX FOR STORING ALL VARIABLES PLUS RESPONSE DATA
% Matrix for storing behavioural data:
% 1st column is subject number
% 2nd column is block number
% 3nd column will record the trial number
% 4rd through 5th column trial information from CondMatrix
% 6th column correct versus incorrect for catch trials (i.e. 1 = correct, 0 = incorrect, 2 = not a catch trial)
% 7th column RT
% 8th column timing for comparing measured with intended timing
% 9th column measure actual ITI
% 10th column trigger value for that trial
output = struct('subject_num',num2cell(repmat(subject_num,1,trialAmountTotal)),...
    'block_num',num2cell(repmat(block,1,trialAmountTotal)),...
    'trial_num',num2cell(1:trialAmountTotal),...
    'sequence',num2cell(condMatrixShuffled(1,:)),...
    'catch',num2cell(condMatrixShuffled(2,:)),...
    'correct_response',num2cell(NaN(1,trialAmountTotal)),...
    'RT',num2cell(NaN(1,trialAmountTotal)),...
    'testTiming',num2cell(NaN(1,trialAmountTotal)),...
    'testITI',num2cell(NaN(1,trialAmountTotal)),...
    'trigger',num2cell(NaN(1,trialAmountTotal)),...
    'time_catch',num2cell(NaN(1,trialAmountTotal))...
    );

% timings in seconds that are constant over trials
VidOnset = 0; % in seconds from trial start
VidOffset = 5;
maxTestTime = 3;% extra time to respond after video stops on catch trials
timeCrossColor = 0.2;

% timings that are random over trials
% inter-trial-interval (ITI) jittered with random uniform distribution between 1.8 and 2.2 seconds
ITI = 1.8 + .4 * rand(1,trialAmountTotal);
for itrial = 1:trialAmountTotal
    output(itrial).ITI = ITI(itrial);
end

% create trigger values between 0 and 255 for MEG and eyetracker
% sequence onset trigger:
%       - 1:14 for sequence number
%       - + 100 for catch trial, except for the ones that have the test at
%       t = 5 sec, they can still be used as normal trials (there are 5 per run)
for itrial = 1:trialAmountTotal
    output(itrial).trigger = output(itrial).sequence + logical(output(itrial).catch) * 100;
    if output(itrial).time_catch == 5
       output(itrial).trigger = output(itrial).trigger - 100;
    end
end

%if not catch trial, response can already be set to 2
for itrial = 1:trialAmountTotal
    if ~output(itrial).catch
        output(itrial).correct_response = 2;
    end
end

%determine time of event in catch trials (i.e., extract from catch trial pool)
catchcount = 0;
for itrial = 1:trialAmountTotal
    if ~output(itrial).catch
        output(itrial).time_catch = 0;
    else
        catchcount = catchcount + 1;
        output(itrial).time_catch = catchTrialPool{2,catchcount};
    end
end

if trackeye
    % INITIALIZE EYELINK
    % obviously specific to eyelink equipment in the MEG lab at CIMeC (in the autumn of 2020)
    
    % It is better not to send too many Eyelink commands to the eye-tracker in a row. For this reason, between them, we wait for a short time, here defined.
    elk.wait = 0.01;
    
    % This code initializes the connection with the eyelink: if something fails, it exit program with error
    if EyelinkInit()~= 1
        error('Eyelink disconnected !!!');
    end
    
    % We need to provide Eyelink with details about the graphics environment and perform some initializations. The initialization information is returned in a
    % structure that also contains useful defaults and control codes (e.g. tracker state bit and Eyelink key values). The structure, moreover, acts as an handle
    % for subsequent commands, like "windowHandle" for Psychtoolbox.
    elk.el = EyelinkInitDefaults(window);
    
    % Here we create the name for the eyelink datafile. 
    % Specific for CIMeC MEG lab:
    % Data gathered from the eye tracker are saved on the eye-tracking PC in a file. Data from all users are
    % saved in the same folder and the folder is routinely cleaned up without any advice. So, be sure to copy your data after the experiment and choose an
    % unique name for the datafile (containing date/time, subject number etc...). It has to be less than 8 characters long.
    % here ivAP0b1 indicates: Ingmar de Vries, Action Prediction, Subject 0, block/run 1
    if subject_num < 10
        Eyefilename = ['ivAP0' num2str(subject_num) 'b' num2str(block)];
    else
        Eyefilename = ['ivAP' num2str(subject_num) 'b' num2str(block)];
    end
    Eyefilename = [Eyefilename '.edf'];
    elk.edfFile = sprintf(Eyefilename);		% Create file name
    Eyelink('Openfile', elk.edfFile);									% Open the file to the eye-tracker
    
    % Writing a short preamble to the file helps if the name became not that informative ;-)
    Eyelink('command', sprintf('add_file_preamble_text ''Ingmar de Vries Project Action Prediction: subject %d ; block %d ; practice %d ; time %s''', subject_num, block, practice, datestr(now, 'YYYYmmddhhMM')));
    
    % Setting the eye-tracker so as to record GAZE of  LEFT and RIGHT eyes, together with pupil AREA
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
    
    % Setting the proper recording resolution, proper calibration type, as well as the data file content
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, screenXpixels - 1, screenYpixels - 1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, screenXpixels - 1, screenYpixels - 1);
    
    % Setting the proper calibration type. Usually we use 9 points calibration. For a long range mount also 13 points (HV13) is a good (longer) calibration.
    % Note that for some subjects the eyetracker had trouble with capturing the corners of the screen, and thus I manually selected the 5 point cross
    % calibration, which excludes the corners. This should be totally fine here, as the video frame is much smaller than the total screen size and
    % restricted to the screen centre anyway. 
    Eyelink('command', 'calibration_type = HV9');
    
    % Setting the proper data file content
    Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
    Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS');
    
    % Setting link data (used for gaze cursor, optional)
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
    Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS');
    
    % Saccade detection thresholds (optional)
    Eyelink('command', 'saccade_velocity_threshold = 35');
    Eyelink('command', 'saccade_acceleration_threshold = 9500');
    
    % Now make sure that we are still connected to the Eyelink ... otherwise throw error
    if Eyelink('IsConnected')~=1
        error('Eyelink disconnected !!!');
    end
    
    % EYELINK CALIBRATION
    % This code allow the EyeLink software to take control of your psychtoolbox screen. This means that at this point you will see participant eyes as recorded
    % by the Eye-tracker camera on the projector in the MSR, a condition essential for setting up the camera. After setting up the camera you can perform calibration
    % and validation at this step.
    
    % Some calibration parameters
    elk.el.foregroundcolour = 0;
    elk.el.backgroundcolour = [bg bg bg] * 255;
    
    % Give eye-tracker control of the screen for camera setup and calibration, until you exit back to psychtoolbox by pressing ESC
    EyelinkDoTrackerSetup(elk.el);
    
end

if MEG
    % INITIALIZE DATAPIXX
    % obviously specific to datapixx equipment in the MEG lab at CIMeC (in the autumn of 2020)
    
    Datapixx('Open');					% Open DataPixx
    
    Datapixx('SetVideoMode', 0);		% This set video mode to normal passthrought, no stereo mode. C24, Straight passthrough from DVI 8-bit RGB to VGA RGB.
    % In this configuration luminance is linear with RGB (see our wiki).
    
    Datapixx('StopAllSchedules');		% Stop all schedules (audio waveforms, triggers etc...)
    
    Datapixx('SetDoutValues', 0);		% Set digital output to zero, as required to prepare for triggering
    
    Datapixx('EnableDinDebounce');		% Enable response debouncing. This is required to prune out spurious button presses after a real response
    
    Datapixx('SetDinLog');				% Clear digital input logger, i.e: clear old responses in the register
    Datapixx('StopDinLog');				% Stop running response logger
    
    Datapixx('RegWrRd');				% So far, no real changes occurred on the physical blue box devices. This command synchronize local and remote registers
    % in a read/write mode and immediately. Only now, the blue box status is as determined by the above initializations.
    
    responseButtonsMask = 2^0 + 2^1 + 2^2 + 2^3;	% Values of response buttons are stored in a cumbersome binary way. This is a binary mask useful to
    % transform them in decimal human-readable values. In particular, red = 1, blue = X, geen = X and yellow =
    % X. It works. Just believe it. I do, I am a true believer. Neo is the one. Follow the rabbit. 
end


% KEYBOARD (ONLY FOR BEHAVIOURAL VERSION)
KbName('UnifyKeyNames');
spaceKey = KbName('space');
escapeKey = KbName('ESCAPE');
upKey = KbName('UpArrow');
downKey = KbName('DownArrow');
RestrictKeysForKbCheck([escapeKey spaceKey upKey downKey]);

%----------------------------------------------------------------------
%                       START WITH INSTRUCTIONS
%----------------------------------------------------------------------

% Instructions
% Note that instructions on the screen were useful to show examples to the subjects, and for making sure instructions were the same between subjects,
% but instructions were actually given verbally by me (Ingmar) for each indivudal subject
if practice
    
    Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
    
    % Reset and fire up the response logger
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');                        % Commit changes to/from DP
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        DrawFormattedText(window,...
            'Thank you for participating in this experiment',...
            'center', yCenter/3, black);
        DrawFormattedText(window,'- Press a button to start the instructions -','center', yCenter+yCenter/2, black);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        DrawFormattedText(window,...
            'Each trial starts with a fixation cross \n\n Whenever you see a fixation cross, please fixate it',...
            'center', yCenter/3, black);
        DrawFormattedText(window,'- Press a button to continue -','center', yCenter+yCenter/2, black);
        Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        DrawFormattedText(window,...
            'Next you will see a video of a ballet dancer \n\n Pay attention to the dancer \n\n Keep fixating the cross while doing so',...
            'center', yCenter/3, black);
        DrawFormattedText(window,'- Press a button to continue -','center', yCenter+yCenter/2, black);
        Screen('DrawTexture', window, stimuli{1,1}, [], stimPos);
        Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        DrawFormattedText(window,...
            'On some trials the fixation cross shortly becomes purple \n\n You have to indicate the color change \n\n by pressing the red button',...
            'center', yCenter/3, black);
        DrawFormattedText(window,'- Press a button to continue -','center', yCenter+yCenter/2, black);
        Screen('DrawTexture', window, stimuli{1,50}, [], stimPos);
        Screen('DrawLines', window, allCoords,lineWidthPix, colors(3,:), [xCenter yCenter], 2);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        DrawFormattedText(window,...
            'On other trials the video is interruped, you will be asked \n\n  whether the dancer is performing a certain move. \n\n Pess red button for YES \n\n Press blue button for NO \n\n Examples of each move follow:',...
            'center', yCenter/3, black);
        DrawFormattedText(window,'- Press a button to continue -','center', yCenter+yCenter/2, black);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        Screen('DrawTexture', window, stimuli{1,(5/6)*refrate}, [], stimPos);
        % Fixation cross
        Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
        % Draw question
        DrawFormattedText(window,'bow?','center', yCenter+yCenter/2, black);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        Screen('DrawTexture', window, stimuli{1,(115/60)*refrate}, [], stimPos);
        % Fixation cross
        Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
        % Draw question
        DrawFormattedText(window,'swan?','center', yCenter+yCenter/2, black);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        Screen('DrawTexture', window, stimuli{1,(195/60)*refrate}, [], stimPos);
        % Fixation cross
        Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
        % Draw question
        DrawFormattedText(window,'pirouette?','center', yCenter+yCenter/2, black);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        Screen('DrawTexture', window, stimuli{3,(190/60)*refrate}, [], stimPos);
        % Fixation cross
        Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
        % Draw question
        DrawFormattedText(window,'passe?','center', yCenter+yCenter/2, black);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        Screen('DrawTexture', window, stimuli{2,(165/60)*refrate}, [], stimPos);
        % Fixation cross
        Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
        % Draw question
        DrawFormattedText(window,'jump?','center', yCenter+yCenter/2, black);
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    while status.newLogFrames == 0
        
        DrawFormattedText(window,...
            'If you need to respond \n\n you will receive direct feedback',...
            'center', yCenter/3, black);
        DrawFormattedText(window,...
            'Correct',...
            'center', yCenter*0.9, colors(1,:));
        DrawFormattedText(window,...
            'OR',...
            'center', yCenter, black);
        DrawFormattedText(window,...
            'Error',...
            'center', yCenter*1.1, colors(2,:));
        DrawFormattedText(window,'- Press a button to continue -','center', yCenter+yCenter/2, black)
        Screen('Flip', window);
        
        Datapixx('RegWrRd');                        % Commit changes to/from DP
        status = Datapixx('GetDinStatus');			% Get response logger status
        
    end
    
    WaitSecs(.5);
    Datapixx('StopDinLog');
    Datapixx('EnableDinDebounce');
    Datapixx('SetDinLog');
    Datapixx('StartDinLog');
    Datapixx('RegWrRd');
    status = Datapixx('GetDinStatus');
    
    DrawFormattedText(window,...
        'Please remember to relax and fixate the cross',...
        'center', 'center', black);
    DrawFormattedText(window,'- When ready, inform the experimenter -','center', yCenter+yCenter/2, black);
    Screen('Flip', window);
    KbStrokeWait;
else% if not a practice run:
    DrawFormattedText(window,...
        'Please remember to relax and fixate the cross',...
        'center', 'center', black);
    DrawFormattedText(window,'- When ready, inform the experimenter -','center', yCenter+yCenter/2, black);
    Screen('Flip', window);
    KbStrokeWait;
end

%countdown at start each block
DrawFormattedText(window, ['The experiment starts in \n\n ' num2str(3)], 'center', 'center', black);
Screen('Flip', window);
WaitSecs(1);
DrawFormattedText(window, ['The experiment starts in \n\n ' num2str(2)], 'center', 'center', black);
Screen('Flip', window);
WaitSecs(1);
DrawFormattedText(window, ['The experiment starts in \n\n ' num2str(1)], 'center', 'center', black);
Screen('Flip', window);
WaitSecs(1);

%----------------------------------------------------------------------
%                       START ACTUAL EXPERIMENT
%----------------------------------------------------------------------

% TRIAL LOOP
correctCount = 0;
catchcount = 0;
if practice
    trialAmountTotal = trialAmountPractice;% stop run before going through all trials on the practice run, but after finishing pre-set trial amount for pratice run
end

for trial = 1:trialAmountTotal
    
    % make sure there's a fixation cross and video frame - placeholder in between trials while current
    % trial settings are being generated and eyelink recording is started
    Screen('FillRect', window, bg, [xCenter-(screenXpixels/2) yCenter-(screenYpixels/2) xCenter+(screenXpixels/2) yCenter+(screenYpixels/2)]);
    Screen('FillRect', window, colorPlaceHolder, [xCenter-(stimWidth/2) yCenter-(stimHeight/2) xCenter+(stimWidth/2) yCenter+(stimHeight/2)]);
    Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
    
    %photo diode
    Screen('FillRect', window, [0 0 0], [0 0 30 30]);
    
    ITIstart = Screen('Flip', window);
    
    % initialize some things and extract from output struct (saves time to do it here during ITI and not while video is being shown)
    response=0;
    respTime = 0;
    rtMarkerTime = 0;
    currentSeq = output(trial).sequence;
    currentCatch = output(trial).catch;
    if currentCatch
        currentCatchType = catchTrialPool{3,currentCatch};
        currentCorrectResp = (2 - catchTrialPool{4,currentCatch})^3;%1 = correct = up/red, 8 = incorrect = down/blue
    end
    currentCatchTime = output(trial).time_catch;
    
    if trackeye
        % EYELINK RECORDING
        
        Eyelink('Message', 'TRIALID %d', trial);
        Eyelink('command', 'record_status_message "SUBJECT = %d ; BLOCK = %d ; TRIAL = %d"', subject_num, block, trial);
        
        % As specified before, it is better not to send to many Eyelink commands to the eye-tracker in a row.  So, after the command, we wait the pre-set time
        WaitSecs(elk.wait);
        
        % Here we start recording eyelink data (left/right gaze and pupil size), preceded by a short pause
        Eyelink('Command', 'set_idle_mode');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %transfer image to host
        [width, height]=Screen('WindowSize', screenNumber);
        imgfile= ([rootdir '\experiment\fixation.bmp']);% this is only for showing fixation cross on eyetracker PC so you get a rough idea of fixation during the experiment
        transferimginfo=imfinfo(imgfile);
        
        fprintf('img file name is %s\n',transferimginfo.Filename);
        
        % image file should be 24bit or 32bit bitmap
        % parameters of ImageTransfer:
        % imagePath, xPosition, yPosition, width, height, trackerXPosition, trackerYPosition, xferoptions
        transferStatus =  Eyelink('ImageTransfer',transferimginfo.Filename,0,0,transferimginfo.Width,transferimginfo.Height,width/2-transferimginfo.Width/2 ,height/2-transferimginfo.Height/2,1);
        if transferStatus ~= 0
            fprintf('*****Image transfer Failed*****-------\n');
        end
        
        WaitSecs(0.1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%
        WaitSecs(elk.wait);
        Eyelink('StartRecording', 1, 1, 1, 1);
        %%
        WaitSecs(elk.wait);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    % Continue with fixation cross in between trials (ITI)
    Screen('FillRect', window, bg, [xCenter-(screenXpixels/2) yCenter-(screenYpixels/2) xCenter+(screenXpixels/2) yCenter+(screenYpixels/2)]);
    Screen('FillRect', window, colorPlaceHolder, [xCenter-(stimWidth/2) yCenter-(stimHeight/2) xCenter+(stimWidth/2) yCenter+(stimHeight/2)]);
    Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
    %photo diode
    Screen('FillRect', window, [0 0 0], [0 0 30 30]);
    
    ITImiddle = Screen('Flip', window);
    WaitSecs(output(trial).ITI - (ITImiddle - ITIstart));% wait for the remaining time of the ITI
    
    % Note that use of WaitSecs is not very accurate,
    % but for the ITI it does not matter that much. It's randomly jittered anyway.
    % What is crucial is the timing of sequence onset and each frame in the sequence
    % That is why we do this with counting exact frames (see below), and by collecting photodiode data for each frame 
    testITI2 = GetSecs;
    output(trial).testITI = testITI2-ITIstart;
    %------------------------------------------------------------------
    %              START ACTUAL TRIAL WITH ACCURATE TIMING
    %------------------------------------------------------------------
    
    disco = 0;% this is the color of the small square shown at the top of the screen (i.e., outside of the video frame). The subject does not see this because it's covered by the photodiode.
    for frame = 1:ceil((VidOffset+maxTestTime) * refrate)
        
        %always draw background color and black placeholder at location of video
        Screen('FillRect', window, bg, [xCenter-(screenXpixels/2) yCenter-(screenYpixels/2) xCenter+(screenXpixels/2) yCenter+(screenYpixels/2)]);
        %         Screen('FillRect', window, colorPlaceHolder, [xCenter-(stimWidth/2) yCenter-(stimHeight/2) xCenter+(stimWidth/2) yCenter+(stimHeight/2)]);
        
        %duration of video
        if frame <= ceil(VidOffset * refrate)
            
            %photo diode
            disco = 1 - disco;%during duration of video let square for photodiode flicker between black and white
            Screen('FillRect', window, [disco disco disco], [0 0 30 30]);
            
            % Draw each frame of the video
            Screen('DrawTexture', window, stimuli{currentSeq,frame}, [], stimPos);
        end
        
        % Fixation cross
        Screen('DrawLines', window, allCoords,lineWidthPix, white, [xCenter yCenter], 2);
        
        % if catch trial, and currentCatchType = 0, fixation cross shortly turns purple
        if currentCatch && ~ischar(currentCatchType) && (frame > ceil(currentCatchTime * refrate)) && (frame <= ceil((currentCatchTime+timeCrossColor) * refrate))
            
            % Now with purple fixation cross
            Screen('DrawLines', window, allCoords,lineWidthPix, colors(3,:), [xCenter yCenter], 2);
            
            % if catch trial, and currentCatchType is a question, then video freezes and you receive the question
        elseif currentCatch && ischar(currentCatchType) && (frame > ceil(currentCatchTime * refrate)) && (frame <= ceil((VidOffset+maxTestTime) * refrate))
            % draw placeholder at location of video
            Screen('FillRect', window, bg, [xCenter-(stimWidth/2) yCenter-(stimHeight/2) xCenter+(stimWidth/2) yCenter+(stimHeight/2)]);
            
            % Draw question
            DrawFormattedText(window,currentCatchType,'center', 'center', black);
        end
        
        % TRIGGERS
        % memory display onset
        if MEG
            if frame == 1
                
                %Datapixx triggering
                triggerPulse = [1 0] .* (output(trial).trigger);
                Datapixx('StopDoutSchedule');
                Datapixx('WriteDoutBuffer', triggerPulse);
                Datapixx('SetDoutSchedule', 1.0/refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                Datapixx('StartDoutSchedule');
                
                % Eyelink triggering
                Eyelink('Message', sprintf('Video onset %d', output(trial).trigger));
                
                % White square for photodiode
                Screen('FillRect', window, [1 1 1], [0 0 30 30]);
                            
                Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
                
            elseif currentCatch && frame == ceil(currentCatchTime * refrate)+1
                                
                %Datapixx triggering
                triggerPulse = [1 0] .* (output(trial).trigger+100);
                Datapixx('StopDoutSchedule');
                Datapixx('WriteDoutBuffer', triggerPulse);
                Datapixx('SetDoutSchedule', 1.0/refrate, 1000, 2);	% Delayed trigger (1/refresh delay rate with ProPixx)
                Datapixx('StartDoutSchedule');
                
                % Eyelink triggering
                Eyelink('Message', sprintf('Video onset %d', output(trial).trigger+100));
                
                Datapixx('EnableDinDebounce');  % Set this to avoid fast oscillation in button press (if unsure use it !)
                
                % Reset and fire up the response logger
                Datapixx('SetDinLog');
                Datapixx('StartDinLog');
                
            end
            
            %Now send all instructions to the Datapixx box, as close as possible to the actual screen flip in PTB
            Datapixx('RegWrVideoSync');
            
        end

        Screen('Flip',window);
        
        if currentCatch && frame == ceil(currentCatchTime * refrate)+1
            % This has to be in its own loop after videosync and
            % screenflip, otherwise causes a delay
            Datapixx('SetMarker');
            Datapixx('RegWrRd');
            rtMarkerTime = Datapixx('GetMarker');
        end
        
        % in catch trials a response is required
        if currentCatch && (frame > ceil(currentCatchTime * refrate))
            
            if MEG
                Datapixx('RegWrRd');						% Commit changes to/from DP
                status = Datapixx('GetDinStatus');			% Get response logger status
                
                if status.newLogFrames > 0					% We've got new data in response buffer !!!
                    
                    [data, time] = Datapixx('ReadDinLog');	% Read data in
                    
                    response = bitand(data(end), responseButtonsMask);
                    respTime = time(end);	%
                    
                    % Eyelink triggering
                    Eyelink('Message', sprintf('Response %d', response));
                    
                    % Send RT/response trigger (response value + 128)
                    Datapixx('EnableDinDebounce');
                    Datapixx('StopDoutSchedule');
                    triggerPulse = [1 0] .* (response+20);
                    Datapixx('WriteDoutBuffer', triggerPulse);
                    Datapixx('SetDoutSchedule', 0, 100, 2);		% No need for programmatic delay here, no wait for projector to fire trigger
                    Datapixx('StartDoutSchedule');
                    Datapixx('RegWr');
                    
                    break
                    
                end
                
            else
                [keyIsDown,respTime,keyCode] = KbCheck;
                if keyCode(upKey)
                    response = 1;
                    if trackeye
                        Eyelink('Message', sprintf('Response %d', response));
                    end
                    break
                elseif keyCode(downKey)
                    response = 2;
                    if trackeye
                        Eyelink('Message', sprintf('Response %d', response));
                    end
                    break
                end
            end
            
        end
        
        %quick and dirty test for timing between first and last frame
        if frame == 1
            test1 = GetSecs;
        end
        if frame == ceil(VidOffset*refrate)
            test2 = GetSecs;
            output(trial).testTiming = test2-test1;
        end
        
        if currentCatch && (frame == ceil(currentCatchTime * refrate))
            testTime = GetSecs;
        end
        
        % check for esc key
        [keyIsDown,respTime,keyCode] = KbCheck;
        if keyCode(escapeKey)
            sca;
            disp('The escape key has been pressed');
            Datapixx('Close');
            return
        end
        
        if ~currentCatch && (frame > ceil(VidOffset * refrate))
            break
        end
    end
    
    if trackeye
        % Stop eyelink
        Eyelink('StopRecording');
    end

    if currentCatch
        % check RT
        output(trial).RT = respTime - rtMarkerTime;
        % respTime and testTime set to 0 at start each trial
        % if respTime is still 0 here no response was given on this trial
        % if testTime is still 0 here this trial was not a catch trial
        
        % feedback
        Screen('FillRect', window, bg, [xCenter-(screenXpixels/2) yCenter-(screenYpixels/2) xCenter+(screenXpixels/2) yCenter+(screenYpixels/2)]);
        
        if response == currentCorrectResp
            
            output(trial).correct_response = 1;
            correctCount = correctCount + 1;
            
            DrawFormattedText(window, 'Correct', 'center', 'center', colors(1,:));
            
        elseif (9-response == currentCorrectResp) && ischar(currentCatchType)
            
            output(trial).correct_response = 0;
            if practice% more elaborate feedback with errors during practice run:
                DrawFormattedText(window,...
                    'Error \n\n Pay attention to the dancer \n\n while fixating the cross',...
                    'center', 'center', colors(2,:));
            else
                DrawFormattedText(window,...
                    'Error',...
                    'center', 'center', colors(2,:));
            end
            
        elseif response == 0 && ~ischar(currentCatchType)
            
            output(trial).correct_response = 0;
            DrawFormattedText(window,...
                'Missed the color change in the fixation cross!',...
                'center', 'center', colors(2,:));
            
        elseif response == 0 && ischar(currentCatchType)
            
            output(trial).correct_response = 0;
            DrawFormattedText(window,...
                'Too slow, please respond faster next time',...
                'center', 'center', colors(2,:));
            
        end

        disp(['response :' num2str(response) ' , and RT: ' num2str(output(trial).RT)]);
            
        Screen('Flip', window);
        WaitSecs(2);
    end
    
    %Save whole output matrix after each trial so that data is stored if the experiment stops unexpectedly
    save([savedir filename], 'output');
    
end

% FEEDBACK AFTER EACH BLOCK
percentCorrect = ceil(correctCount/(sum(logical(condMatrixShuffled(2,1:trialAmountTotal))))*100);

DrawFormattedText(window,['You had ' num2str(percentCorrect) '% correct in this block'],'center', yCenter/2, black);
DrawFormattedText(window,'Please tell the experimenter that you finished this block \n\n and take a break','center', yCenter+yCenter/2, black);
Screen('Flip', window);
KbStrokeWait;

% CLOSE STUFF
% Eyelink
% Download data file
if trackeye
    try
        fprintf('Receiving data file ''%s''\n', elk.edfFile );
        status = Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2 == exist(elk.edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', elk.edfFile, savedir );
        end
    catch
        fprintf('Problem receiving data file ''%s''\n', elk.edfFile );
    end
    
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    Eyelink('ShutDown');
end

if MEG
    % Close DataPixx
    Datapixx('StopAllSchedules');		% Stop all schedules
    Datapixx('SetDoutValues', 0);       % Reset triggers to zero
    Datapixx('StopDinLog');             % Stop response buttons recording
    Datapixx('Close');                  % Close DataPixx
end

% Psychtoolbox
Priority(0);
sca

disp(['You had ' num2str(percentCorrect) '% correct in this block']);
if practice
    disp(['Refresh rate is ' num2str(refrate)]);% display in command window just as sanity check after practice run before starting real runs
end

