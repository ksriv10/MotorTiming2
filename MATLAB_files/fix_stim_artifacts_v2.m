function [] = fix_stim_artifacts_v2()

if isunix
        !ls **.cbin > batchfile
        %!ls *1701.19*.not.mat > batchfile
    else
        !dir /B **.cbin > batchfile
    end
fid=fopen('batchfile','r');

while 1
    fn=fgetl(fid);if (~ischar(fn));break;end
    disp(fn);
    [dat,fs]=evsoundin('',fn,'obs0');
    [trig_dat,fs]=evsoundin('',fn,'obs1');
    % time=1/fs:1/fs:length(dat)/fs;
    % filt_dat=bandpass_filtfilt(dat,fs,1,50,'hanningfir');
    % figure(1)
    % ax(1)=subplot(4,1,1);
    % plot(time,dat)
    % hold on
    ind_big_increases=find(diff(trig_dat)>5000);
    ind_big_decreases=find(diff(trig_dat)<-5000);
    if ~isempty(ind_big_increases)
        ind_big_increases=[ind_big_increases(diff(ind_big_increases)>5);ind_big_increases(end)];
    end
    if ~isempty(ind_big_decreases)
        ind_big_decreases=[ind_big_decreases(diff(ind_big_decreases)>5);ind_big_decreases(end)];
    end
%     if ~isempty(ind_big_increases)
%         ind_big_increases=[ind_big_increases(1); ind_big_increases(ind_temp+1)];
%         if ind_big_increases(1)==1
%             ind_big_increases=ind_big_increases(2:end);
%         end
%     end
% 
%     ind_temp=find(diff(ind_big_decreases)>5);
%     if ~isempty(ind_big_decreases)
%         ind_big_decreases=[ind_big_decreases(1); ind_big_decreases(ind_temp+1)];
%         if ind_big_decreases(1)==1
%             ind_big_decreases=ind_big_decreases(2:end);
%         end
%     end
%     
%     for k=1:length(ind_big_increases)
%         temp=find(ind_big_decreases>ind_big_increases(k) & ind_big_decreases<ind_big_increases(k)+20, 1);
% %         temp_ind_big_decreases=ind_big_decreases(ind_big_decreases>ind_big_increases(k));
% %         temp=find((ind_big_increases(k)-temp_ind_big_decreases)>-20,1);
%         if isempty(temp)
%             ind_big_decreases=[ind_big_decreases; ind_big_increases(k)+16];
%         end
%     end
%     ind_big_decreases=sort(ind_big_decreases);
%     
%     for k=1:length(ind_big_decreases)
%         temp=find(ind_big_increases<ind_big_decreases(k) & ind_big_increases>ind_big_decreases(k)-20, 1);
% %         temp=find(abs(ind_big_increases-ind_big_decreases(k))<20,1);
%         if isempty(temp)
%             ind_big_increases=[ind_big_increases; ind_big_decreases(k)-16];
%         end
%     end
%     ind_big_decreases=sort(ind_big_decreases);
%     if length(ind_big_decreases)<length(ind_big_increases)

%        temp=find((ind_big_increases(1:length(ind_big_decreases))-ind_big_decreases)<-20,1);
% %        t_problem=time(ind_big_increases(temp));
%        while 1
%            ind_big_decreases=[ind_big_decreases; ind_big_increases(temp)+16];
%            ind_big_decreases=sort(ind_big_decreases);
%            if length(ind_big_increases)==length(ind_big_decreases)
%                break;
%            else
%                temp=find((ind_big_increases(1:length(ind_big_decreases))-ind_big_decreases)<-20,1);
%                if isempty(temp)
%                    break;
%                end
%            end
%        end
%     elseif length(ind_big_decreases)>length(ind_big_increases)
%        temp=find((ind_big_increases-ind_big_decreases(1:length(ind_big_increases)))<-20,1);
% %        t_problem=time(ind_big_increases(temp));    
%        while 1
%           ind_big_increases=[ind_big_increases; ind_big_decreases(temp)-16];
%           ind_big_increases=sort(ind_big_increases);
%           if length(ind_big_increases)==length(ind_big_decreases)
%                break;
%           else
%               temp=find((ind_big_increases-ind_big_decreases(1:length(ind_big_increases)))<-20,1);
%               if isempty(temp)
%                   break;
%               end
%           end
%        end
%     end
    % plot(time(ind_big_increases),dat(ind_big_increases),'ro',time(ind_big_decreases),dat(ind_big_decreases),'go')

    % ax(2)=subplot(4,1,2);
    % plot(time,filt_dat)
    % hold on
    % plot(time(ind_big_increases),filt_dat(ind_big_increases),'ro',time(ind_big_decreases),filt_dat(ind_big_decreases),'go')
%        time=1/fs:1/fs:length(dat)/fs; 
%        figure;
%        plot(time,dat)
%        hold on
%        plot(time(ind_big_increases),dat(ind_big_increases),'ro',time(ind_big_decreases),dat(ind_big_decreases),'go')

    for j=1:length(ind_big_increases)
%         dat(ind_big_increases(j)+1:ind_big_decreases(j))=dat(ind_big_increases(j)+1:ind_big_decreases(j))-(dat(ind_big_decreases(j))-dat(ind_big_increases(j)));
%         dat(ind_big_increases(j)+1:ind_big_increases(j)+3)=mean([dat(ind_big_increases(j)) dat(ind_big_increases(j)+4)]);
%         dat(ind_big_decreases(j)+1:ind_big_decreases(j)+3)=mean([dat(ind_big_decreases(j)) dat(ind_big_decreases(j)+4)]);
        temp=ind_big_decreases(j)+10-ind_big_increases(j)+1;
        dat(ind_big_increases(j):ind_big_decreases(j)+10)=interp1([1 temp],[dat(ind_big_increases(j)) dat(ind_big_decreases(j)+10)],1:1:temp,'linear');
%         dat(ind_big_increases(j):ind_big_decreases(j)+1)-(dat(ind_big_decreases(j)-1)-dat(ind_big_increases(j)-1));
%         dat(ind_big_increases(j):ind_big_increases(j)+4)=mean([dat(ind_big_increases(j)-1) dat(ind_big_increases(j)+5)]);
%         dat(ind_big_decreases(j):ind_big_decreases(j)+4)=mean([dat(ind_big_decreases(j)-1) dat(ind_big_decreases(j)+5)]);
    end
    % ax(3)=subplot(4,1,3);
    % plot(time,dat)
    % hold on
    % plot(time(ind_big_increases),dat(ind_big_increases),'ro',time(ind_big_decreases),dat(ind_big_decreases),'go')
%     time=1/fs:1/fs:length(dat)/fs;
%     plot(time,dat)
    % filt_dat=bandpass_filtfilt(dat,fs,1,50,'hanningfir');
    % ax(4)=subplot(4,1,4);
    % plot(time,filt_dat)
    % hold on
    % plot(time(ind_big_increases),filt_dat(ind_big_increases),'ro',time(ind_big_decreases),filt_dat(ind_big_decreases),'go')
    % linkaxes(ax,'x')
    savename=[fn(1:end-4) 'mat'];
    if exist(savename,'file')==2
        save(savename,'dat','trig_dat','fs','-v7.3','-append');
    else
        save(savename,'dat','trig_dat','fs','-v7.3');
    end
end
end

