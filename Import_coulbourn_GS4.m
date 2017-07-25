function [Rat] = Import_coulbourn_GS4(varargin)
%This converts coulbourn text files to a list of times with variable
%numbers. These can be saved as .txt files or called from another script. State_Key must be an excel file in the source folder.
%The second varargin is the first letters (usually initials) of the data
%files from coulbourn
%   Detailed explanation goes here


%This should be in milliseconds
pre = -40000;
post = 40000;

%This gives the time of a Coulbourn time unit in milliseconds
time_unit = 1;

%The varargin sets the Source_Folder to be used to retrieve the data
if nargin > 1
    
    Source_Folder = varargin{1};
else    

% User_Input = inputdlg({'What folder has the data to rename (Make sure the data is in the separated stage folder)?'},'Rename the Data');
% Source_Folder = User_Input{1};
%Source_Folder = 'E:\OE\OE DIFF RECORDING\tst_files\Deval';
%Source_Folder = 'E:\Unblock Satiation 4';
%Source_Folder = 'E:\OE\OE DIFF RECORDING\tst_files\Run1 TST';
Source_Folder = 'E:\TouchscreenTask\Data\Pavlovian Deval';

end

% if nargin > 2
% 
%     String_Identifier = varargin{2};
% else    
%     
% String_Identifier = {'CLB'};
% end

%Note that StageXX MUST be on all of the files

File_End = '.csv';

%This gets all the files from the directory
AllFiles = dir(Source_Folder);

% %This sorts the files by animal and run
% [~, Ind] = sort([AllFiles(:).datenum]);
% 
% AllFiles = AllFiles(Ind);

%This puts all of the filenames into a cell array of sessions, with only
%.txt files and files which begin with group
Filenames = {AllFiles.name};

Filenames = Filenames(~cellfun(@isempty, (strfind(Filenames, File_End))));

  %  & (~cellfun(@isempty, (strfind(Filenames, String_Identifier{1})))); %| ~cellfun(@isempty, (strfind(Filenames, String_Identifier{2}))));

%TXT_Files = find(TXT_File);

%This checks to make sure StageXX is on all of the files
% if any(cellfun(@isempty, (strfind(Filenames(TXT_File), 'Stage'))))
%         
%         display('Some files do not have Stage in the name')
% end
% 
% if any(cellfun(@(x) strfind(x, 'Stage') + 6 > numel(x) - 4, Filenames(TXT_File)))
%     display('"Stage" needs to have two alphanumerics following it in the filename. Check the filenames')
%     return
% end
% 
% %all files need to be in the same format, this gets the protocol number positions and
% %the stage number positions
% protocol_pos = strfind(Filenames{find(TXT_File,1)}, String_Identifier{1}) + numel(String_Identifier{1});
% protocol_pos = [protocol_pos protocol_pos + 1];
% 
% stage_pos = strfind(Filenames{find(TXT_File,1)}, 'Stage') + numel('Stage');
% stage_pos = [stage_pos stage_pos + 1];

%This finds the earliest date within the TXT_Files. All session and stage
%information is based on this date. It is important not to include other
%sessions within the folder which are not part of the experiment
datenums = cellfun(@(x) datenum(str2double(x(1:4)),str2double(x(6:7)),str2double(x(9:10)) ,str2double(x(13:14)), str2double(x(16:17)), str2double(x(19:20))), Filenames);

first_day = min(datenums);

%display(sprintf('The date of the first file in the source folder is %s', datestr(first_day)));

%This sorts TXT_File by date
[~, inds] = sort(datenums);

Filenames = Filenames(inds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This gets the design from an excel sheet, the first sheet contains the
%animals run and the dates for each session and stage
[~, ~, Main_Design] = xlsread(strcat(Source_Folder, '\Design.xlsx'),1);

%First determine the number of stages and the sessions per stage. Stage
%contains a vector indicating the stage number. Sessions is a vector with
%the number of sessions per stage. Sess_Num is the session number chronologically within
%a given day. Sess_Num is used to keep track of which session is which from
%the filename. 

Stage = Main_Design(2:end, 1);
Date = Main_Design(2:end, 2);
Sess_Use = cellfun(@str2double, Main_Design(2:end, 3:end));
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

%sess_analyzed keeps track of sessions which have been analyzed. This is
%needed for sessions which fall on the same day so that the same session is
%not analyzed twice
sess_analyzed = false(size(Filenames));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The first loop loops through rat #
for i = 1:numel(Rat_Num)
    
    %This keeps the design file location in the Rat structure
    Rat(i).Design_Location = strcat(Source_Folder, '\Design.xlsx');
    
    %This sets the rat's name
    Rat(i).Name = Rat_Num{i};

    %This finds the text files associated with the current rat. Note that
    %the 'numel(Rat_Num{i})' keeps track of how many characters are in the
    %rat name.
     rat_dates = cellfun(@(x) strcmp((x(end - 3 - numel(Rat_Num{i}):end - 4)),Rat_Num{i}), Filenames);
    
    %This keeps track of which session the loops are on
    sess_track = 1;
    
    %The second loop goes through each stage. sheets is the number of
    %stages entered in the design file.
    for j = 1:numel(sheets) - 1
        
        %This leaves the for loop when sess_track reaches the end of the
        %Main_Design rows
        if sess_track > numel(Stage)
            break
        end
        
        %This finds the curent stage
        curr_stage = Stage{sess_track,1};
        
        %This if statement exists incase stages are entered on a
        %sheet but not in the main design sheet
        if isempty(strfind(sheets{j + 1},curr_stage))
            display(sprintf('Stage %s is entered on the sheet, but not in the main design. This stage has been skipped.', sheets{j + 1}));
            continue
        end    
        
        %This sets the Stage_Name
        Rat(i).Stage(j).Stage_Name = curr_stage; %#ok<*AGROW>
        
                %n_stage keeps track of the number of sessions within a
                %stage
                n_stage = 1;
                
                TS = struct;
                    TS.State_On = [];
                    TS.State_Off = [];
                    TS.Response_On = [];
                    TS.Response_Off = [];
                    TS.Missing = [];
                    
                %This loops through each day of the stage, if the new
                %session is not in the current stage the while loops ends
                while (sess_track <= numel(Stage) && strcmp(Stage{sess_track}, curr_stage))  
                
                %This checks whether the current session should be analyzed
                if Sess_Use(sess_track, i) ~= 1
                    TS(n_stage).Missing = 1;
                   % display(sprintf('The File for rat %s from %s is missing, session is skipped', Rat_Num{i}, Date{sess_track}));
                    n_stage = n_stage + 1;
                    sess_track = sess_track + 1;
                    continue
                end
                
                %The filename is in format: yyyy_mm_dd, this finds the
                %correct filename. Note that sessions which occur on the
                %same day MUST be put in chrnological order in the
                %Main_Design. This allows us to track which sessions have
                %already been analyzed for a given day if two sessions come
                %from the same day

                file_date = datestr(datenum(Date{sess_track}), 'yyyy_mm_dd');
                
                filename = find(rat_dates & ~sess_analyzed & strncmp(file_date, Filenames, 10), 1);
                
                sess_analyzed(filename) = true;
                
                curr_file = strcat(Source_Folder, '\', Filenames{filename});
                
                
                if isempty(filename)
                    
                    TS(n_stage).Missing = 1;
                    display(sprintf('The File for rat %s from %s is missing, session is skipped', Rat_Num{i}, Date{sess_track}));
                    n_stage = n_stage + 1;
                    sess_track = sess_track + 1;
                    
                    continue
                end    
                
%                 %Thisfinds the filename for the current rat and stage, 
%                 protocol_num = Main_Design{j + 1, strcmp(Rat_Num(i), Main_Design(1, :))};
%                 
%                 %This checks for the protocol number and stage number
%                 if ~any(cellfun(@(x) strcmp(x(protocol_pos),protocol_num), Filenames(TXT_File)) & cellfun(@(x) strcmp(x(stage_pos),Stage{j}), Filenames(TXT_File)))
%                     %This skips if there is no file
%                     continue
%                     
%                 %This makes sure there is only one file with the protocol and stage number    
%                 elseif sum(cellfun(@(x) strcmp(x(protocol_pos),protocol_num), Filenames(TXT_File)) & cellfun(@(x) strcmp(x(stage_pos),Stage{j}), Filenames(TXT_File))) > 1
%                     
%                     display(sprintf('There is more than one file with the protocol number %s and the stage number %s. The stage was skipped', protocol_num, Stage{j}))
%                     continue
%                 end
%                 curr_file_log = TXT_File;    
%                 curr_file_log(TXT_File == 1) = (cellfun(@(x) strcmp(x(protocol_pos),protocol_num), Filenames(TXT_File)) & cellfun(@(x) strcmp(x(stage_pos),Stage{j}), Filenames(TXT_File))); 
%                 
%                 curr_file = strcat(Source_Folder, '\', Filenames{curr_file_log});
%                 

                 %This gets the timestamps for each session
                 TS(n_stage) = get_ts(curr_file, Rat_Num{i}); 
                 
                 n_stage = n_stage + 1;
                 sess_track = sess_track + 1;
                end
                
                %If the user has added in all of the dates in the design
                %file before getting all of the files, the TS structure
                %will be empty at this point, this skips the 
                
%                  if numel(TS) > Sessions(j)
%                      
%                      display(sprintf('There are more sessions in the file than were set in the design file for rat %d and stage %d from file: %s', i, j, curr_file))
%                  
%                  elseif numel(TS) ~= Sessions(j)
%                      
%                      display(sprintf('There were only %d sessions imported for rat %s in stage %s', numel(TS), Rat_Num{i}, Stage{j}))
%                      
%                  end    
 
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
        
                %This loops through all of the sessions from the file
                for k = 1:numel(TS)
 
                    if TS(k).Missing == 1
                        Rat(i).Stage(j).Session(k).Skip = true;
                        continue
                    else
                        Rat(i).Stage(j).Session(k).Skip = false;
                    end
                    
                    %This loops through each of the trial types
                    for l = 1:numel(trial_type)
                    
%                         if any(strcmp(trial_type{l}, {'V1', 'V2', 'V3'})) && j == 5
%                             
%                             continue
%                         end    
                            
                    %This puts the times of the state entries into the
                    %data structure. If the cue occurs multiple times per
                    %trial (re-entries to the cue) then the first onset of
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
                    
                    %The following lines were added for data in which there
                    %were issues with the sessions
%                     if (i == 10 || i == 18) && j == 5
%                     Rat(i).Stage(j).Session(k).(trial_type{l}).Times_Off_Offsets(:,p) = (cue_offs(p:cue_reps(l):numel(cue_offs)) - cue_onsets(1:6))*time_unit;  
%                     else
%                     Rat(i).Stage(j).Session(k).(trial_type{l}).Times_Off_Offsets(:,p) = (cue_offs(p:cue_reps(l):numel(cue_offs)) - cue_onsets)*time_unit;  
%                     end
%                     end
                    
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
                    Rat(i).Stage(j).Session = orderfields(Rat(i).Stage(j).Session);
                end
    end
end

                  
end

end
  
function [TS] = get_ts(curr_file, rat_num) 

%This opens the file
fileID = fopen(curr_file);

%This gets the relevant data. Note that %d is int32, column 12 is the first
%switch on
Behavior_Data = textscan(fileID, '%*s %s %s %d %*d %*d %s %f %d %d %d %d %d %d %d %d %d %d %d', 'HeaderLines', 1, 'Delimiter', ',');        

TS = struct;

    %This checks whether the rat at the session start is the rat of
    %interest. Note that if a rat named 1 and a rat named 10,11,12, etc.
    %are on the same file, we can't just look at the first character. For
    %this reason we need to check if the number is followed by an
    %underscore. The first if statement makes sure the name has an
    %underscore at the end of the input.
    %if ~strncmp(Behavior_Data{4}{2}, rat_num, numel(rat_num))
    if str2num(Behavior_Data{4}{2} ) ~= str2num(rat_num) 
        display('Theres a problem with the way the user entered the rat names for this script. The script can be updated in order to accept this input, but it needs to be done before continuing')
        return
    end    
    
    %This sets the time to milliseconds rather than seconds
    Behavior_Data{5} = Behavior_Data{5}*1000;
    
    %This keeps the timestamps
    TS.State_On = [];
    TS.State_Off = [];
    TS.Response_On = [];
    TS.Response_Off = [];
        
    state_onoff_log = Behavior_Data{7} > 0;
    
    TS.State_On = [TS.State_On; [Behavior_Data{5}(state_onoff_log) Behavior_Data{7}(state_onoff_log)]];
    
    TS.State_Off = [TS.State_Off; [Behavior_Data{5}(state_onoff_log) Behavior_Data{6}(state_onoff_log)]];
    
    %This loops through the number of possible switches for Coulbourn
    for j = 2%1:4
        %8
%         if any(str2num(rat_num) == [61 62 63 65 71 72 76 78])
%             j = 1;
%         end    
        window_log = Behavior_Data{8 + j} == 1;
        
        TS.Response_On = [TS.Response_On; [Behavior_Data{5}(window_log) int32(ones(sum(window_log), 1))]];
                                    %12
        window_log = Behavior_Data{12 + j} == 1;
        
        TS.Response_Off = [TS.Response_Off; [Behavior_Data{5}(window_log) int32(ones(sum(window_log), 1))]];
    end
    
    TS.Response_On = sortrows(TS.Response_On, 1);
    TS.Response_Off = sortrows(TS.Response_Off, 1);
    TS.Missing = 0;

fclose(fileID);
end

    
  
    
    
       


