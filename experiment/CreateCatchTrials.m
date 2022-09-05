function [CatchTrialsKinematics, CatchTrialsSemantics, CatchTrialsKeep] = CreateCatchTrials()

%% Create pool of all potential catch trials
CatchTrialsKinematics = cell(4,14);
CatchTrialsSemantics = cell(4,14);
CatchTrialsKeep = cell(4,14);
% first dimension (4 rows) = 
%   1. sequence number (15x)
%   2. stop time (sec)
%   3. question 
%   4. correct answer (1 = correct, 0 = incorrect) 

%% kinematics
% sequence 1
CatchTrialsKinematics{1,1} = 1;
CatchTrialsKinematics{2,1} = 0.8;
CatchTrialsKinematics{3,1} = 'body moving \n\n down?';
CatchTrialsKinematics{4,1} = 1;

% sequence 3
CatchTrialsKinematics{1,2} = 3;
CatchTrialsKinematics{2,2} = 0.4;
CatchTrialsKinematics{3,2} = 'arms moving \n\n down?';
CatchTrialsKinematics{4,2} = 0;

% sequence 3 
CatchTrialsKinematics{1,3} = 3;
CatchTrialsKinematics{2,3} = 3.7;
CatchTrialsKinematics{3,3} = 'arms moving \n\n up?';
CatchTrialsKinematics{4,3} = 0;

% sequence 4
CatchTrialsKinematics{1,4} = 4;
CatchTrialsKinematics{2,4} = 4.4;
CatchTrialsKinematics{3,4} = 'left leg moving \n\n up?';
CatchTrialsKinematics{4,4} = 1;

% sequence 5
CatchTrialsKinematics{1,5} = 5;
CatchTrialsKinematics{2,5} = 3.1;
CatchTrialsKinematics{3,5} = 'body moving \n\n down?';
CatchTrialsKinematics{4,5} = 1;

% sequence 6
CatchTrialsKinematics{1,6} = 6;
CatchTrialsKinematics{2,6} = 0.5;
CatchTrialsKinematics{3,6} = 'hands moving \n\n down?';
CatchTrialsKinematics{4,6} = 0;

% sequence 7 
CatchTrialsKinematics{1,7} = 7;
CatchTrialsKinematics{2,7} = 1;
CatchTrialsKinematics{3,7} = 'hands moving \n\n up?';
CatchTrialsKinematics{4,7} = 0;

% sequence 8
CatchTrialsKinematics{1,8} = 8;
CatchTrialsKinematics{2,8} = 3.2;
CatchTrialsKinematics{3,8} = 'body moving \n\n down?';
CatchTrialsKinematics{4,8} = 1;

% sequence 9
CatchTrialsKinematics{1,9} = 9;
CatchTrialsKinematics{2,9} = 3.1;
CatchTrialsKinematics{3,9} = 'left arm moving \n\n up?';
CatchTrialsKinematics{4,9} = 1;

% sequence 9
CatchTrialsKinematics{1,10} = 9;
CatchTrialsKinematics{2,10} = 4.1;
CatchTrialsKinematics{3,10} = 'arms moving \n\n down?';
CatchTrialsKinematics{4,10} = 0;

% sequence 10
CatchTrialsKinematics{1,11} = 10;
CatchTrialsKinematics{2,11} = 1.1;
CatchTrialsKinematics{3,11} = 'left leg moving \n\n up?';
CatchTrialsKinematics{4,11} = 0;

% sequence 12
CatchTrialsKinematics{1,12} = 12;
CatchTrialsKinematics{2,12} = 0.6;
CatchTrialsKinematics{3,12} = 'right leg moving \n\n up?';
CatchTrialsKinematics{4,12} = 1;

% sequence 13
CatchTrialsKinematics{1,13} = 13;
CatchTrialsKinematics{2,13} = 4.3;
CatchTrialsKinematics{3,13} = 'body moving \n\n up?';
CatchTrialsKinematics{4,13} = 0;

% sequence 14
CatchTrialsKinematics{1,14} = 14;
CatchTrialsKinematics{2,14} = 0.8;
CatchTrialsKinematics{3,14} = 'right leg moving \n\n down?';
CatchTrialsKinematics{4,14} = 1;

%% semantics
% sequence 1
CatchTrialsSemantics{1,1} = 1;
CatchTrialsSemantics{2,1} = 3.2;
CatchTrialsSemantics{3,1} = 'pirouette?';
CatchTrialsSemantics{4,1} = 1;

% sequence 2
CatchTrialsSemantics{1,2} = 2;
CatchTrialsSemantics{2,2} = 2.7;
CatchTrialsSemantics{3,2} = 'pirouette?';
CatchTrialsSemantics{4,2} = 0;

% sequence 2
CatchTrialsSemantics{1,3} = 2;
CatchTrialsSemantics{2,3} = 4.7;
CatchTrialsSemantics{3,3} = 'pirouette?';
CatchTrialsSemantics{4,3} = 1;

% sequence 4
CatchTrialsSemantics{1,4} = 4;
CatchTrialsSemantics{2,4} = 1;
CatchTrialsSemantics{3,4} = 'jump?';
CatchTrialsSemantics{4,4} = 0;

% sequence 5
CatchTrialsSemantics{1,5} = 5;
CatchTrialsSemantics{2,5} = 3.4;
CatchTrialsSemantics{3,5} = 'bow?';
CatchTrialsSemantics{4,5} = 1;

% sequence 6
CatchTrialsSemantics{1,6} = 6;
CatchTrialsSemantics{2,6} = 4.2;
CatchTrialsSemantics{3,6} = 'bow?';
CatchTrialsSemantics{4,6} = 0;

% sequence 7
CatchTrialsSemantics{1,7} = 7;
CatchTrialsSemantics{2,7} = 4.7;
CatchTrialsSemantics{3,7} = 'swan?';
CatchTrialsSemantics{4,7} = 1;

% sequence 8
CatchTrialsSemantics{1,8} = 8;
CatchTrialsSemantics{2,8} = 1.8;
CatchTrialsSemantics{3,8} = 'pirouette?';
CatchTrialsSemantics{4,8} = 0;

% sequence 10
CatchTrialsSemantics{1,9} = 10;
CatchTrialsSemantics{2,9} = 4.5;
CatchTrialsSemantics{3,9} = 'jump?';
CatchTrialsSemantics{4,9} = 0;

% sequence 11
CatchTrialsSemantics{1,10} = 11;
CatchTrialsSemantics{2,10} = 0.5;
CatchTrialsSemantics{3,10} = 'jump?';
CatchTrialsSemantics{4,10} = 1;

% sequence 11
CatchTrialsSemantics{1,11} = 11;
CatchTrialsSemantics{2,11} = 4.6;
CatchTrialsSemantics{3,11} = 'pirouette?';
CatchTrialsSemantics{4,11} = 0;

% sequence 12
CatchTrialsSemantics{1,12} = 12;
CatchTrialsSemantics{2,12} = 4.5;
CatchTrialsSemantics{3,12} = 'passe?';
CatchTrialsSemantics{4,12} = 1;

% sequence 13
CatchTrialsSemantics{1,13} = 13;
CatchTrialsSemantics{2,13} = 1.5;
CatchTrialsSemantics{3,13} = 'passe?';
CatchTrialsSemantics{4,13} = 1;

% sequence 14
CatchTrialsSemantics{1,14} = 14;
CatchTrialsSemantics{2,14} = 2.7;
CatchTrialsSemantics{3,14} = 'pirouette?';
CatchTrialsSemantics{4,14} = 0;

%% questions to keep in the analysis, i.e. with question after video, i.e., catch trial time at t = 5 sec
% sequence 1
CatchTrialsKeep{1,1} = 1;
CatchTrialsKeep{2,1} = 5;
CatchTrialsKeep{3,1} = 'arms moving \n\n down?';
CatchTrialsKeep{4,1} = 1;

% sequence 2
CatchTrialsKeep{1,2} = 2;
CatchTrialsKeep{2,2} = 5;
CatchTrialsKeep{3,2} = 'swan?';
CatchTrialsKeep{4,2} = 0;

% sequence 3
CatchTrialsKeep{1,3} = 3;
CatchTrialsKeep{2,3} = 5;
CatchTrialsKeep{3,3} = 'pirouette?';
CatchTrialsKeep{4,3} = 0;

% sequence 4
CatchTrialsKeep{1,4} = 4;
CatchTrialsKeep{2,4} = 5;
CatchTrialsKeep{3,4} = 'pirouette?';
CatchTrialsKeep{4,4} = 0;

% sequence 5
CatchTrialsKeep{1,5} = 5;
CatchTrialsKeep{2,5} = 5;
CatchTrialsKeep{3,5} = 'pirouette?';
CatchTrialsKeep{4,5} = 0;

% sequence 6
CatchTrialsKeep{1,6} = 6;
CatchTrialsKeep{2,6} = 5;
CatchTrialsKeep{3,6} = 'jump?';
CatchTrialsKeep{4,6} = 0;

% sequence 7
CatchTrialsKeep{1,7} = 7;
CatchTrialsKeep{2,7} = 5;
CatchTrialsKeep{3,7} = 'jump?';
CatchTrialsKeep{4,7} = 0;

% sequence 8
CatchTrialsKeep{1,8} = 8;
CatchTrialsKeep{2,8} = 5;
CatchTrialsKeep{3,8} = 'left leg moving \n\n down?';
CatchTrialsKeep{4,8} = 1;

% sequence 9
CatchTrialsKeep{1,9} = 9;
CatchTrialsKeep{2,9} = 5;
CatchTrialsKeep{3,9} = 'arms moving \n\n down?';
CatchTrialsKeep{4,9} = 1;

% sequence 10
CatchTrialsKeep{1,10} = 10;
CatchTrialsKeep{2,10} = 5;
CatchTrialsKeep{3,10} = 'left arm moving \n\n down?';
CatchTrialsKeep{4,10} = 1;

% sequence 11
CatchTrialsKeep{1,11} = 11;
CatchTrialsKeep{2,11} = 5;
CatchTrialsKeep{3,11} = 'bow?';
CatchTrialsKeep{4,11} = 1;

% sequence 12 
CatchTrialsKeep{1,12} = 12;
CatchTrialsKeep{2,12} = 5;
CatchTrialsKeep{3,12} = 'arms moving \n\n down?';
CatchTrialsKeep{4,12} = 1;

% sequence 13 
CatchTrialsKeep{1,13} = 13;
CatchTrialsKeep{2,13} = 5;
CatchTrialsKeep{3,13} = 'arms moving \n\n down?';
CatchTrialsKeep{4,13} = 0;

% sequence 14 
CatchTrialsKeep{1,14} = 14;
CatchTrialsKeep{2,14} = 5;
CatchTrialsKeep{3,14} = 'pirouette?';
CatchTrialsKeep{4,14} = 1;

% save('C:\Users\ingmar.devries\Google Drive\Active projects\ProjectAction\experiment\CatchTrialPool','CatchTrialsKinematics','CatchTrialsSemantics','CatchTrialsKeep');
end



