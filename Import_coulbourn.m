function [Rat] = Import_coulbourn(varargin)
%This converts coulbourn text files to a list of times with variable
%numbers. These can be saved as .txt files or called from another script. State_Key must be an excel file in the source folder.
%The second varargin is the first letters (usually initials) of the data
%files from coulbourn
%   Detailed explanation goes here


%This should be in milliseconds
pre = -40000;
post = 40000;

%This gives the time of a Coulbourn time unit in milliseconds
time_unit = 20;

%The varargin sets the Source_Folder to be used to retrieve the data
if nargin > 1
    
    Source_Folder = varargin{1};
else    

% User_Input = inputdlg({'What folder has the data to rename (Make sure the data is in the separated stage folder)?'},'Rename the Data');
% Source_Folder = User_Input{1};
%Source_Folder = 'E:\OE\OE DIFF RECORDING\tst_files\Deval';
%Source_Folder = 'E:\unblock satiation';
%Source_Folder = 'E:\Sensory Specific Devaluation';
%Source_Folder = 'E:\OE\OE DIFF RECORDING\tst_files\Run1 TST';
%Source_Folder = 'E:\Unblocking_October_2016';
Source_Folder = 'E:\Satiety Revaluation\Run1';
%Source_Folder = 'E:\Unblocking_June_2017';

end

if nargin > 2

    String_Identifier = varargin{2};
else    
    
String_Identifier = {'CLB'};
end

%Note that StageXX MUST be on all of the files

File_End = '.txt';

%This gets all the files from the directory
AllFiles = dir(Source_Folder);

% %This sorts the files by animal and run
% [~, Ind] = sort([AllFiles(:).datenum]);
% 
% AllFiles = AllFiles(Ind);

%This puts all of the filenames into a cell array of sessions, with only
%.txt files and files which begin with group
Filenames = {AllFiles.name};

TXT_File = ~cellfun(@isempty, (strfind(Filenames, File_End)))...
    & (~cellfun(@isempty, (strfind(Filenames, String_Identifier{1})))); %| ~cellfun(@isempty, (strfind(Filenames, String_Identifier{2}))));

%TXT_Files = find(TXT_File);

%This checks to make sure StageXX is on all of the files
if any(cellfun(@isempty, (strfind(Filenames(TXT_File), 'Stage'))))
        
       
        
        if any(cellfun(@isempty, (strfind(Filenames(TXT_File), 'STAGE'))))
            
            for i = find(TXT_File)
                
                if ~isempty(strfind(Filenames{i}, 'STAGE'))
                    
                    start_i = strfind(Filenames{i}, 'STAGE');
                    New_Name = strcat(Filenames{i}(1:start_i - 1), 'Stage', Filenames{i}(start_i + 5:end));
                    
                    commandStr = sprintf('rename "%s" "%s"',[Source_Folder '\' Filenames{i}],New_Name);
                    dos(commandStr);
                end
            end    
        else
        display('Some files do not have Stage in the name')
        end
end

if any(cellfun(@(x) strfind(x, 'Stage') + 6 > numel(x) - 4, Filenames(TXT_File)))
    display('"Stage" needs to have two alphanumerics following it in the filename. Check the filenames')
    return
end

%all files need to be in the same format, this gets the protocol number positions and
%the stage number positions
protocol_pos = strfind(Filenames{find(TXT_File,1)}, String_Identifier{1}) + numel(String_Identifier{1});
protocol_pos = [protocol_pos protocol_pos + 1];

stage_pos = strfind(Filenames{find(TXT_File,1)}, 'Stage') + numel('Stage');
stage_pos = [stage_pos stage_pos + 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This gets the design from an excel sheet, the first sheet contains the
%animals run and the dates for each session and stage
[~, ~, Main_Design] = xlsread(strcat(Source_Folder, '\Design.xlsx'),1);

%First determine the number of stages and the sessions per stage. Stage
%contains a vector indicaing the stage number. Sessions is a vector with
%the number of sessions per stage. Runs is how many runs occurred within
%a given day. If runs is more than one, the protocols need to be listed for
%each animal
Stage = Main_Design(2:end, 1);
Sessions = cellfun(@(x) str2double(x), Main_Design(2:end, 2));
%Runs = cellfun(@(x) str2double(x), Main_Design(2:end, 3));
Rat_Num = Main_Design(1,3:end);

%This determines which Coulbourn state is associated with each cue for each
%stage from the design file. Each array within State_Key is the design
%information for each stage

%Incase not all of the sheets are entered, this checks how many sheets
%there are
[~, sheets] = xlsfinfo(strcat(Source_Folder, '\Design.xlsx'));
for i = 1:numel(sheets) - 1
    
    [~,~, State_Key{i}] = xlsread(strcat(Source_Folder, '\Design.xlsx'),i + 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The first loop loops through rat #
for i = 1:numel(Rat_Num)

    %This keeps the design file location in the Rat structure
    Rat(i).Design_Location = strcat(Source_Folder, '\Design.xlsx');
    
    %This sets the rat's name
    Rat(i).Name = Rat_Num{i};
    
    %The second loop goes through each stage. sheets is the number of
    %stages entered in the design file.
    for j = 1:numel(sheets) - 1
        
        %This if statement exists incase stages are entered on a
        %sheet but not in the main design sheet
        if numel(Stage) < j
            continue
        end    
        
        if j == 5
            a = 1;
        end    
        %This sets the Stage_Name
        Rat(i).Stage(j).Stage_Name = Stage{j}; %#ok<*AGROW>
        
                %Thisfinds the filename for the current rat and stage, the
                %rat is found by searching for the protocol from the design
                %by looking at the raw cell array. Since all of the
                %sessions for a given rat, stage and run should be from the
                %same protocol we can use the Main_Design to search for the file with all of the sessions.
                %The protocol number in the Main_Design is separated by
                %commas for each run.
                %This gets the protocol number from the Main_Design
                protocol_num = Main_Design{j + 1, strcmp(Rat_Num(i), Main_Design(1, :))};
                
                %This checks for the protocol number and stage number
                if ~any(cellfun(@(x) strcmp(x(protocol_pos),protocol_num), Filenames(TXT_File)) & cellfun(@(x) strcmp(x(stage_pos),Stage{j}), Filenames(TXT_File)))
                    %This skips if there is no file
                    continue
                    
                %This makes sure there is only one file with the protocol and stage number    
                elseif sum(cellfun(@(x) strcmp(x(protocol_pos),protocol_num), Filenames(TXT_File)) & cellfun(@(x) strcmp(x(stage_pos),Stage{j}), Filenames(TXT_File))) > 1
                    
                    display(sprintf('There is more than one file with the protocol number %s and the stage number %s. The stage was skipped', protocol_num, Stage{j}))
                    continue
                end
                curr_file_log = TXT_File;    
                curr_file_log(TXT_File == 1) = (cellfun(@(x) strcmp(x(protocol_pos),protocol_num), Filenames(TXT_File)) & cellfun(@(x) strcmp(x(stage_pos),Stage{j}), Filenames(TXT_File))); 
                
                curr_file = strcat(Source_Folder, '\', Filenames{curr_file_log});
                
                 %This gets the timestamps for each session
                 [TS] = get_ts(curr_file, Rat_Num{i}, protocol_num); 
                  
                 if numel(TS) > Sessions(j)
                     
                     display(sprintf('There are more sessions in the file than were set in the design file for rat %d and stage %d from file: %s', i, j, curr_file))
                 
                 elseif numel(TS) ~= Sessions(j)
                     
                     display(sprintf('There were only %d sessions imported for rat %s in stage %s', numel(TS), Rat_Num{i}, Stage{j}))
                     
                 end    
 
            %This checks the design for the trial-type and state numbers.
            %If all is written in the design, each animal does not have to
            %be entered separately
            if all(strcmpi('all',State_Key{j}(2:end,1)))% & i ~= 1
                
                trial_type = State_Key{j}(2:end, 4);
                cue_reps = cell2mat(State_Key{j}(2:end, 3));
                states = cell2mat(State_Key{j}(2:end, 2));
                
                cue_identity = State_Key{j}(2:end, 5);
                
            else
                
                %The rat number entered can either be a string or number,
                %the cellfun converts everything to string no matter
                %whether it was a sstring or a number
                curr_rat = [false; strcmpi(Rat_Num(i),cellfun(@num2str, State_Key{j}(2:end,1), 'UniformOutput', false))];
            
                trial_type = State_Key{j}(curr_rat, 4);
                cue_reps = cell2mat(State_Key{j}(curr_rat, 3));
                states = cell2mat(State_Key{j}(curr_rat, 2));
                
                cue_identity = State_Key{j}(2:end, 5);
            end
           if i == 2 && j == 7
               
               a = 1
           end    
                %This loops through all of the sessions from the file
                for k = 1:numel(TS)
                    
                    %This loops through each of the trial types
                    for l = 1:numel(trial_type)
                    
                    %This puts the times of the state entries into the
                    %data structure. If the cue occurs multiple times per
                    %trial (re-entries to the cue) than the first onset of
                    %the cue is set as time zero. All responses are
                    %referenced to this time.
                    cue_ons = TS(k).State_On(TS(k).State_On(:,2) == states(l), 1);
                    cue_onsets = cue_ons(1:cue_reps(l):numel(cue_ons));
                    
                    Rat(i).Stage(j).Session(k).(trial_type{l}).Times = (cue_onsets)*time_unit;
                    
                    cue_offs = TS(k).State_Off(TS(k).State_Off(:,2) == states(l), 1);
                    
                    for p = 1:cue_reps(l)
                    
                    if cue_reps(l) > 1
                      %This provides the offsets of each of the cue ons and cue offs from the first cue on.  
                      Rat(i).Stage(j).Session(k).(trial_type{l}).Times_On_Offsets(:,p) = (cue_ons(p:cue_reps(l):numel(cue_ons)) - cue_onsets)*time_unit;
                      
                    end
                    
                    Rat(i).Stage(j).Session(k).(trial_type{l}).Times_Off_Offsets(:,p) = (cue_offs(p:cue_reps(l):numel(cue_offs)) - cue_onsets)*time_unit;  
                    
                    end
                    
                    %This provides the cue identity
                    Rat(i).Stage(j).Session(k).(trial_type{l}).Identity = cue_identity{l};
                            
                            for m = 1:4
                                %This gets the timestamps for each of the
                                %switch inputs (on and off) for the current trial_type
                                
                                Rat(i).Stage(j).Session(k).(trial_type{l}).(sprintf('S%d_On', m)) = ...
                                    Time_Stamps((TS(k).Response_On(TS(k).Response_On(:,2) == m, 1))*time_unit, (cue_onsets)*time_unit, pre, post);
                                
                                Rat(i).Stage(j).Session(k).(trial_type{l}).(sprintf('S%d_Off', m)) = ...
                                    Time_Stamps((TS(k).Response_Off(TS(k).Response_Off(:,2) == m, 1))*time_unit, (cue_onsets)*time_unit, pre, post);
                                
                            end
                    end
                    %This puts all of the fields in alphanumeric order.
                    %This makes looping through trialtype fields for each stage be in
                    %the same order
                    if ~isempty(Rat(i).Stage(j).Session)
                        Rat(i).Stage(j).Session = orderfields(Rat(i).Stage(j).Session);
                    end
                end
    end
end

                  
end


  
function [TS] = get_ts(curr_file, rat_num, protocol_num) 

%This opens the file
fileID = fopen(curr_file);

%This gets the relevant data. Note that %d is int32
Behavior_Data = textscan(fileID, '%*s %s %s %d %*d %*d %s %d %d %d %d %d %d %d %d %d %d %d %d', 'HeaderLines', 1, 'Delimiter', '\t');        


%This gets the protocol number and checks to make sure it was correctly
%named in the filename. They're converted to numbers since the protocol
%number in the coulbourn files do not have zeros in front of 1-9
%the first cell is the protocol

if ~(str2double(Behavior_Data{2}(2)) == str2double(protocol_num))
    
    display('The protocol in the filename did not match the protocol in the file') 
end



%This determines the starts of sessions by checking for time restart (note
%that session numbers (column 3, don't always increase)
session_starts = [true; diff(Behavior_Data{5}) < 0];

%incase sessions are initiated by accident but are not run, sessions with
%less than 10 rows are removed.
session_start_ind = find(session_starts);

sessions_rem = [false; diff(session_start_ind) < 10];

session_start_ind(sessions_rem) = 0;

TS = struct;
%This loops through each session to get timestamps for the current rat
n = 0;
for i = 1:numel(session_start_ind)
    
    %This checks whether the rat at the session start is the rat of
    %interest. Note that if a rat named 1 and a rat named 10,11,12, etc.
    %are on the same file, we can't just look at the first character. For
    %this reason we need to check if the number is followed by an
    %underscore. The first if statement makes sure the name has an
    %underscore at the end of the input.
    if ~strcmp(Behavior_Data{4}{session_start_ind(i)}(end), '_')
        
        display('Theres a problem with the way the user entered the rat names for this script. The script can be updated in order to accept this input, but it needs to be done before continuing')
        return
    end    
     if ~strncmp(Behavior_Data{4}{session_start_ind(i)}, strcat(rat_num, '_'), numel(rat_num) + 1)
%         
%         
%         %This was added because the subjects were input incorrectly for an
%         %experiment. This should be removed
%         if any(str2double(rat_num) == 21:26)
%             
%             if ~(str2double(Behavior_Data{4}{session_start_ind(i)}(1)) == str2double(rat_num) - 20 && strcmp(Behavior_Data{4}{session_start_ind(i)}(2), '_'))
%         
%         %%%%%!!!!!!!!!!!!REMOVE TWO ABOVE IF STATEMENTS YET LEAVE CONTINUE         
%         continue
%               
%             end
%         else
            continue 
%         end 
%    
%      
     end
    n = n + 1;
    
    %This keeps the timestamps
    TS(n).State_On = [];
    TS(n).State_Off = [];
    TS(n).Response_On = [];
    TS(n).Response_Off = [];
    
    if i == numel(session_start_ind)
        
        sess_end = numel(session_starts);
    else 
        sess_end = session_start_ind(i + 1);
    end


    state_onoff_log = false(numel(session_starts),1);
    
    state_onoff_log(session_start_ind(i): sess_end) = Behavior_Data{7}(session_start_ind(i): sess_end) > 0;
    
    TS(n).State_On = [TS(n).State_On; [Behavior_Data{5}(state_onoff_log) Behavior_Data{7}(state_onoff_log)]];
    
    TS(n).State_Off = [TS(n).State_Off; [Behavior_Data{5}(state_onoff_log) Behavior_Data{6}(state_onoff_log)]];
    
    %This loops through the number of possible switches for Coulbourn
    for j = 1%
        
        window_log = false(numel(session_starts),1);
        window_log(session_start_ind(i): sess_end) = Behavior_Data{8 + j}(session_start_ind(i): sess_end) == 1;
        
        TS(n).Response_On = [TS(n).Response_On; [Behavior_Data{5}(window_log) int32(ones(sum(window_log), 1))]];
        
        window_log = false(numel(session_starts),1);
        window_log(session_start_ind(i): sess_end) = Behavior_Data{12 + j}(session_start_ind(i): sess_end) == 1;
        
        TS(n).Response_Off = [TS(n).Response_Off; [Behavior_Data{5}(window_log) int32(ones(sum(window_log), 1))]];
    end
    
    TS(n).Response_On = sortrows(TS(n).Response_On, 1);
    TS(n).Response_Off = sortrows(TS(n).Response_Off, 1);
    TS(n).Date = Behavior_Data{1}(session_start_ind(i));
end

fclose(fileID);
end

    
  
    
    
       


