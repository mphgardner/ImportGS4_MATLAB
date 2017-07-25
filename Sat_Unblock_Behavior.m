function [Results] = Sat_Unblock_Behavior(Run, cue_window, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%This should be in milliseconds
pre = -40000;
post = 40000;
%cue_window = [-30000 0];
%varargin exists incase on wants to adjust the session numbers to match
%other runs or stages. For example, in the overexpectation task, the first
%run of the experiment has 12 sessions in stage 1 whereas the second run
%has 2 sessions in stage 1. Since we'd like to include the 2 sessions of
%the second run with sessions 11 and 12 of the first run, we can declare
%here that the two sessions should be idenitified as session 11 and 12.
%This would be accomplished by inputting varargin as Session where Session
%is:
%Session{1} = {1:12, 1:4, 1} Here the first stage of run 1 has 12
%sessions, the second stage has 4, and the third stage has 1
%Session{2} = {11:12, 1:4, 1} Now the first stage of run 2 has 2 sessions but
%are considered as the 11th and 12 session

%The varargin needs to be a cell vector the length of the runs input
if nargin > 2 && (numel(varargin{1}) ~= numel(Run) || nargin > 3)
    
    display('There is an error with the varargin')
end    
    
    


%This preallocates the Results structure
% for j = 1:run
%     for i = 1:stage{j}
%         
%         if isempty(session{j}{i}) 
%             continue
%         end    
%         Results.Run(j).Stage(i).Pokes.S1 = NaN(rat, numel(Events{i}), max(Trials{i}), session{j}{i}(end));
%         Results.Run(j).Stage(i).Resp_Perc.S1 = NaN(rat, numel(Events{i}), max(Trials{i}), session{j}{i}(end));
%         Results.Run(j).Stage(i).Pokes.S2 = NaN(rat, numel(Events{i}), max(Trials{i}), session{j}{i}(end));
%         Results.Run(j).Stage(i).Resp_Perc.S2 = NaN(rat, numel(Events{i}), max(Trials{i}), session{j}{i}(end));
%         Results.Run(j).Stage(i).Pokes.Both = NaN(rat, numel(Events{i}), max(Trials{i}), session{j}{i}(end));
%         Results.Run(j).Stage(i).Resp_Perc.Both = NaN(rat, numel(Events{i}), max(Trials{i}), session{j}{i}(end));
%         
%         Results.Run(j).Stage(i).Events = Events{i};
%         
%     end
% end

%This does a short loop to find the length of the stages
% stage_length(1:numel(Run)) = 0;
% for h = 1:numel(Run)
%     for i = 1:numel(Run(h).Rat)
%         if numel(Run(h).Rat(i).Stage) > stage_length(h)
%             
%             stage_length(h) = numel(Run(h).Rat(i).Stage);
%         end
%     end
% end

fieldname_flag = [];
 
        %%%Main Loop%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First loop through runs
for h = 1:numel(Run)
    
    if isempty(Run(h))
        
        continue
    end
    
    %loop through rats
    for i = 1:numel(Run(h).Rat)
        
        if isempty(Run(h).Rat(i).Stage)
            
            continue
        end
        

        for j = 1:numel(Run(h).Rat(i).Stage)

            if isempty(Run(h).Rat(i).Stage(j).Session)
                
                continue
            end
            if j == 5
                
                a = 1
            end    
            %Here we check whether the user input special numbering for the
            %sessions. The dfault is 1 to the number of sessions
            if nargin > 2
                
                session_nums = varargin{1}{h}{j};
            
            else
                
                session_nums = 1:numel(Run(h).Rat(i).Stage(j).Session);
            end
            
                
            for k = 1:numel(Run(h).Rat(i).Stage(j).Session)
                
                if isempty(Run(h).Rat(i).Stage(j).Session(k)) %|| Run(h).Rat(i).Stage(j).Session(k).Skip == 1
                    
                    continue
                end
                
                %THese keep track of the number of multiple inpokes or
                %outpokes in a row per day. These errors should be rare.
                %This just keeps track of the rate of these errors.
                %Results.Run(h).Stage(j).S1_Error(i,session{h}{j}(k)) = 0;
                %Results.Run(h).Stage(j).S2_Error(i,session{h}{j}(k)) = 0;
                %errors1 = 0;
                %errors2 = 0;
                
                %This if statement checks to make sure all the events are
                %there for each session and reorders the fieldnames incase
                %they're in different orders for each rat. This if
                %statement is IMPORTANT don't delete. The order of events
                %aren't necessarily always the same
                Temp_Events = fieldnames(Run(h).Rat(i).Stage(j).Session(k));
                Temp_Events = Temp_Events(cellfun(@isempty, (strfind(Temp_Events, 'Times'))) & cellfun(@isempty, (strfind(Temp_Events, 'Identity')))  & cellfun(@isempty, (strfind(Temp_Events, 'Skip'))));
                
                if numel(fieldname_flag) < j || isempty(Temp_Events)
                    fieldname_flag(j) = 1;
                    Events{j} = Temp_Events;
                    
                elseif ~all(strcmp(Temp_Events,Events{j}))
                    display(sprintf('The session from rat %d, stage %d, session %d, does not have matching events with other sessions from the same stage.', i,j,k))
                end
                %This just removes any fields which contain time
                %information or the Identity of the cue
                if i == 1 && j == 7
                    a = 1;
                end
                
                for l = 1:numel(Events{j})
                    
                    if isempty(Run(h).Rat(i).Stage(j).Session(k).(Events{j}{l}))
                        
                        continue
                    end
                    
                    %This loops through the four possible switch inputs
                    for p = 1:4
                    
                    S_On = cell(1);
                    S_Off = cell(1);

                    %This for loop removes any double ins or outs, and sets
                    %an in at the beginning if the first event is an out.
                    %Likewise if the last event is an in.

                    poketrain = [];
                    for m = 1:numel(Run(h).Rat(i).Stage(j).Session(k).(Events{j}{l}).(sprintf('S%d_On', p)))
                        
                        %REMOVE
                        if k == 3 && m == 16 && l == 2
                           a = 1; % break
                       end
                        
                        [poketrain(m,:)] = Check_On_Off(Run(h).Rat(i).Stage(j).Session(k).(Events{j}{l}).(sprintf('S%d_On', p)){m}, Run(h).Rat(i).Stage(j).Session(k).(Events{j}{l}).(sprintf('S%d_Off', p)){m}, pre, post);
                       % errors1 = errors1 + err1;
                        
                        if flag == 1
                            display(sprintf('Error occured for rat %d, stage %d, session %d, event %s, on switch %d on trial %d', i,j,k,  Events{j}{l}, p, m))
                        end
                        
                        
                    end
                      
                    %this sets 1 millisecond bins
                    %[response] = signal_on_off(S_On, S_Off, cue_window(1), cue_window(2), 1);
                    
                    response = poketrain(:, -pre + cue_window(1) + 1: -pre + cue_window(2));
                    
                    
                   % responseboth = response1 | response2;
                 
                   % Results.Run(h).Stage(j).S1_Error(i,session{h}{j}(k)) = Results.Run(h).Stage(j).S1_Error(i,session{h}{j}(k)) + errors1;
                   % Results.Run(h).Stage(j).S2_Error(i,session{h}{j}(k)) = Results.Run(h).Stage(j).S2_Error(i,session{h}{j}(k)) + errors2;
                    
                    trial_num = size(response,1);
                    
                    %If the animal is poking at the beginning of the bin,
                    %it can be counted as a poke or not. The first line
                    %below is not counted, the second is counted
                    
                    %Results.Run(h).Stage(j).Pokes.(sprintf('S%d', p))(i,l,1:trial_num,session_nums(k))= sum(diff(response,1,2) > 0, 2)*(60000/(cue_window(2) - cue_window(1)));
                    
                    Results.Run(h).Stage(j).Pokes.(sprintf('S%d', p))(i,l,1:trial_num,session_nums(k))= sum(diff([false(size(response,1),1) response],1,2) > 0, 2)*(60000/(cue_window(2) - cue_window(1)));
                    
                    Results.Run(h).Stage(j).Resp_Perc.(sprintf('S%d', p))(i,l,1:trial_num,session_nums(k)) = (sum(response,2)/size(response,2))*100;
                    Results.Run(h).Stage(j).Event_Names{l} = Events{j}{l};
                   % Results.Run(h).Stage(j).Pokes.Both(i,l,1:trial_num,session{h}{j}(k)) = sum(diff(responseboth,1,2) < 0, 2);
                   % Results.Run(h).Stage(j).Resp_Perc.Both(i,l,1:trial_num,session{h}{j}(k)) = (sum(responseboth,2)/size(responseboth,2))*100;
                    
                    %This determines the period in which the nosepoke
                    %detector is detecting at a high frequency.
                    bin = 1000;
                    freq_check = 2;
                    high_freq = [];
                    
                    Ons = [false(size(response,1),1) diff(response,1,2) < 0];
                    
                    for t = 1:size(response,2)/bin
                    
                        high_freq(:,t) = sum(Ons(:, (t-1)*bin + 1: t*bin),2) > freq_check;
                    end
                    
                    Results.Run(h).Stage(j).High_Freq.(sprintf('S%d', p))(i,l,1:trial_num,session_nums(k)) = sum(high_freq, 2);
                   
                    
                    clear S_On S_Off
                    end    
                end              
            end
        end
    end
end

end


% function [inpoke, outpoke, flag, errors]= Check_On_Off(inpoke, outpoke, pre, post)
% errors = 0;
% flag = 0;
% if (isempty(inpoke) &&  isempty(outpoke)) || (numel(inpoke) == numel(outpoke) && all(inpoke == outpoke))
%     return
% end
%     
% if ~isempty(inpoke) &&  isempty(outpoke) 
%     outpoke = post;  
%     
% elseif isempty(inpoke) &&  ~isempty(outpoke) 
%     inpoke = pre;
% end
% 
% if size(inpoke,1) > 1
%     inpoke = inpoke';
% end
% 
% if size(outpoke,1) > 1
%     outpoke = outpoke';
% end        
% if inpoke(1) > outpoke(1)
%     inpoke = [pre inpoke];
% end
% 
% if inpoke(end) > outpoke(end)
%     outpoke = [outpoke post];
% end
% 
% not_resolving = 0;
% 
% if numel(inpoke) ~= numel(outpoke)
%     %errors = 1;
% end
% 
% %multiple inpokes or outpokes in a row need to be removed. If this occurs
% %for both inpokes and outpokes, the number of elements of the two may be
% %the same. A systematic way of checking for multiple inpokes and outpokes
% %is performed below
% 
% %while numel(inpoke) ~= numel(outpoke) || not_resolving == 10
% flag_out = 0;
% while flag_out == 0 && not_resolving < 50
%     
%     for j = 1: numel(outpoke)
%         if j <= numel(inpoke) && outpoke(j) - inpoke(j) <= 0
%             
%             if j == 1
%                 outpoke(j) = [];
%             else
%                 outpoke(j - 1) = [];
%                 
%             end
%             errors = errors + 1;
%             break
%             
%         elseif j < numel(inpoke) && inpoke(j + 1) - outpoke(j) <= 0
%             
%             inpoke(j + 1) = [];
%             errors = errors + 1;
%             break
%             
%         elseif j > numel(inpoke)
%             outpoke(j-1) = [];
%                   
%             break
%             
%         end     
% 
%     end
%     if numel(inpoke) > numel(outpoke)
%              
%         outpoke(end) = post;
%         errors = errors + 1;
%         
%     elseif numel(outpoke) == numel(inpoke) && any(outpoke < inpoke) 
%         not_resolving = not_resolving + 1;
%         
%     elseif numel(outpoke) == numel(inpoke) 
%         flag_out = 1;
%     else    
%         not_resolving = not_resolving + 1;
%     end
% end
% if not_resolving == 50
%          flag = 1;
% end
% end



function [response]= Check_On_Off(inpoke, outpoke, pre, post)

%Set the nosepokes to all false (unpoked)
response = false(1, post - pre);

if (isempty(inpoke) &&  isempty(outpoke)) || (numel(inpoke) == numel(outpoke) && all(inpoke == outpoke))
    return
end
    
if ~isempty(inpoke) &&  isempty(outpoke) 
    outpoke = post;  
    
elseif isempty(inpoke) &&  ~isempty(outpoke) 
    inpoke = pre;
end

if size(inpoke,2) > 1
    inpoke = inpoke';
end

if size(outpoke,2) > 1
    outpoke = outpoke';
end        
if inpoke(1) > outpoke(1)
    inpoke = [pre; inpoke];
end

if inpoke(end) == post
    
    inpoke(end) = [];
end

if outpoke(end) == post
    
    outpoke(end) = [];
end

ins = [inpoke ones(size(inpoke))];
outs = [outpoke 2*ones(size(outpoke))];

allp = sortrows([ins; outs],1);

%last keeps track of whether the prior time was on off or an on (its
%defaulted to off)
last = false;
last_ind = 1;

inds = allp(:,1) - pre + 1;
for i = 1:size(allp,1)
    
    response(last_ind: inds(i)) = last;
    if i < numel(inds)

        last = allp(i,2) == 1;
        last_ind = inds(i);
        
    else
  
        response(inds(i):end) = allp(i,2) == 1;
    end
end
end
        
        
    
    
