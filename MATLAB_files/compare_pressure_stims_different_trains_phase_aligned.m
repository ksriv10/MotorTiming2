function [] = compare_pressure_stims_different_trains_phase_aligned(load_dat)
%compares the pressure waveform when no stim is delivered to when stim is
%delivered. When you rerun code with saved file, you can specify delay of
%stims to focus on for waveforms

if nargin==0
    load_dat=0;
end
pre_time=0.3; %in secs  %for viewing waveform
post_onset_time=.55; %in secs
fs=32000;
size_guess=10000;
clip_dur=fix(fs*(pre_time+post_onset_time))+1;
current_dir=pwd;
list=dir([current_dir '\*.csv']);
stimtrain_csv=list.name;
all_stim_trains_file=csvread(stimtrain_csv);
all_stim_trains=all_stim_trains_file(:,1);
all_stim_times=all_stim_trains_file(:,2);
all_stim_times_within_file=all_stim_trains_file(:,3);
all_stim_file_numbers=all_stim_trains_file(:,4);
stim_train_counter=1;
% stim_trains_saved=[];
stim_trains_saved=zeros(size_guess,1);
% pressure_period_catch=[];
pressure_period_catch=zeros(size_guess,1);
% pressure_clips_stim=[];
pressure_clips_stim=zeros(size_guess,clip_dur);
% stim_latencies=[]; %gonna be in secs
stim_latencies=zeros(size_guess,1); %gonna be in secs
% stim_phases=[]; %gonna be in radians
stim_phases=zeros(size_guess,1); %gonna be in radians
% pressure_clips_catch=[];
pressure_clips_catch=zeros(size_guess,clip_dur);
% pressure_clips_catch_w_prestim=[];
% catch_comparison_index=[];
catch_comparison_index=zeros(size_guess,1);
catch_counter=1;
stim_counter=1;

% if ~isempty(strfind(stimtrain_csv,'100ms'))
%     latency=0.1; %in seconds,
% else
%     latency=0.05;
% end
latency=0.025;

pre_time_stim=50; %in ms  %cutoffs relative to expiration onset to be classified as stim
post_onset_time_stim=300; %in ms
total_data_points=0;
load('calibration_later_file.mat')
cal_matrix=[0 pressure_zero;pressure_cal pressure_one];
slope=(cal_matrix(1,1)-cal_matrix(2,1))/(cal_matrix(1,2)-cal_matrix(2,2));

if ~load_dat
    if isunix
        !ls **.not.mat > batchfile
        %!ls *1701.19*.not.mat > batchfile
    else
        !dir /B **.not.mat > batchfile
    end
    fid=fopen('batchfile','r');
    
    while 1
        fn=fgetl(fid);dat_fn=fn(1:end-8);if (~ischar(fn));break;end
        load(fn);disp(fn);
        ind_per=strfind(dat_fn,'.');
        file_number=str2num(dat_fn(ind_per(1)+1:ind_per(2)-1));
        load(dat_fn,'dat','fs')
        dat=dat';
        
%         stim_times_old=compute_stim_times(dat_fn);
%         stim_times=(all_stim_times(all_stim_file_numbers==(file_number+1))-total_data_points)*1000/fs;
        stim_times=(all_stim_times_within_file(all_stim_file_numbers==(file_number+1)))*1000/fs;
%         file_stim_trains=all_stim_trains(stim_train_counter:stim_train_counter+length(stim_times)-1);
        file_stim_trains=all_stim_trains(all_stim_file_numbers==(file_number+1));
        stim_train_counter=stim_train_counter+length(stim_times);
        total_data_points=total_data_points+length(dat);
        pressure=0.249174*slope*(dat-cal_matrix(1,2));
%         [pressure,sp,t,f,filtsong]=SmoothData_pressure(pressure,fs,1,0,512,0.8,2.0,1,50);
%         pressure=bandpass_filtfilt(pressure,fs,1,50,'hanningfir');
        for j=1:length(onsets)
            samp_on=fix((onsets(j)/1000-pre_time)*fs);
            samp_off=fix(samp_on+(pre_time+post_onset_time)*fs);
            if samp_on>0 && samp_off<length(pressure)
                if isempty(stim_times(stim_times>onsets(j)-pre_time_stim & stim_times<onsets(j)+post_onset_time_stim));
                    pressure_clips_catch(catch_counter,:)=pressure(samp_on:samp_off);
                    if j<length(onsets)
                        pressure_period_catch(catch_counter)=(onsets(j+1)-onsets(j))/1000;
                    else
                        pressure_period_catch(catch_counter)=NaN;
                    end
%                     pressure_clips_catch_w_prestim=[pressure_clips_catch_w_prestim;pressure(samp_on:samp_off)];
                    catch_counter=catch_counter+1;
                else
                    ind_temp=find(stim_times>onsets(j)-pre_time_stim & stim_times<onsets(j)+post_onset_time_stim);
                    stim_time_temp=stim_times(ind_temp);
                    if length(stim_time_temp)==1
                        pressure_clips_stim(stim_counter,:)=pressure(samp_on:samp_off);
                        stim_trains_saved(stim_counter)=file_stim_trains(ind_temp);
                        stim_latencies(stim_counter)=(stim_time_temp-onsets(j))/1000;
                        if j<length(onsets)
                            stim_phases(stim_counter)=stim_latencies(stim_counter)*2*pi*1000/(onsets(j+1)-onsets(j));
                        else
                            stim_phases(stim_counter)=stim_latencies(stim_counter)*2*pi*1000/(onsets(j)-onsets(j-1));
                        end
%                         samp_off1=fix(stim_time_temp/1000*fs)-1;
%                         catch_comparison_index(end+1)=max(1,size(pressure_clips_catch_w_prestim,1));
                        catch_comparison_index(stim_counter)=max(1,catch_counter-1);
%                         pressure_clips_catch_w_prestim=[pressure_clips_catch_w_prestim;[pressure(samp_on:samp_off1) NaN(1,(samp_off-samp_off1))]];
                        stim_counter=stim_counter+1;
                    end
                end
            end
        end
    end
    pressure_clips_catch=pressure_clips_catch(1:catch_counter-1,:);
    pressure_period_catch=pressure_period_catch(1:catch_counter-1);
    pressure_clips_stim=pressure_clips_stim(1:stim_counter-1,:);
    stim_trains_saved=stim_trains_saved(1:stim_counter-1);
    stim_latencies=stim_latencies(1:stim_counter-1);
    stim_phases=stim_phases(1:stim_counter-1);
    catch_comparison_index=catch_comparison_index(1:stim_counter-1);
%     pressure_clips_catch=pressure_clips_catch_w_prestim;
    pre_aligned_time=0.01;
    time2=-pre_aligned_time:1/fs:0.1;
    stim_aligned=NaN(length(stim_latencies),length(time2));
    stim_aligned_single_catch_sub=NaN(length(stim_latencies),length(time2));
%     z=size(pressure_clips_catch_w_prestim,1);
    stim_effect_latencies=zeros(size(stim_latencies));
    for k=1:length(stim_latencies)
%         temp=pressure_clips_stim(k,:)-nanmean(pressure_clips_catch_w_prestim);
        temp=pressure_clips_stim(k,:)-nanmean(pressure_clips_catch);
        id_start=fix((pre_time+stim_latencies(k)-pre_aligned_time)*fs);
        temp2=temp(id_start:min(id_start+length(time2)-1,length(temp)));
        stim_aligned(k,1:length(temp2))=temp2;
        i=catch_comparison_index(k);
%         p=0;
%         while 1
%             if sum(isnan(pressure_clips_catch(i,id_start:min(id_start+3300,length(temp)))))==0
                if i<6
                    ind_catch_sub=1:10;
                elseif i>(size(pressure_clips_catch,1)-5)
                    ind_catch_sub=(size(pressure_clips_catch,1)-10):(size(pressure_clips_catch,1)-1);
                else
                    ind_catch_sub=i-4:i+5;
                end
                pressure_neighbors=NaN(10,length(time2));
                for m=1:10
                    samp_on=fix((pre_time+stim_phases(k)*pressure_period_catch(ind_catch_sub(m))/(2*pi))*fs);
                    samp_off=samp_on+length(time2)-1;
                    if samp_on>0 && samp_off<length(pressure_clips_catch(ind_catch_sub(m),:))
                        pressure_neighbors(m,:)=pressure_clips_catch(ind_catch_sub(m),samp_on:samp_off);
                    end
                end
                control_mean=nanmean(pressure_neighbors);
                temp3=pressure_clips_stim(k,id_start:min(id_start+length(time2)-1,length(temp)));
                temp3=temp3-control_mean(1:length(temp3));
                stim_aligned_single_catch_sub(k,1:length(temp3))=temp3;

                
%                 break;
%             elseif i==1
%                 i=catch_comparison_index(k)+1;
%                 p=1;
%             elseif p==0
%                 i=i-1;
%             elseif i==z
%                 break;
%             else
%                 i=i+1;
%             end
%                 
%         end
        disp(k);
        
%         stim_line=pressure_clips_stim(k,fix((pre_time+stim_latencies(k))*fs):end)-pressure_clips_stim(k,fix((pre_time+stim_latencies(k))*fs));
%         pressure_clips_catch_temp=pressure_clips_catch(:,fix((pre_time+stim_latencies(k))*fs):end)...
%             -repmat(pressure_clips_catch(:,fix((pre_time+stim_latencies(k))*fs)),1,length(pressure_clips_catch(1,fix((pre_time+stim_latencies(k))*fs):end)));
%         thresh=prctile(pressure_clips_catch_temp,95);
%         m=1;
%         while m<length(stim_line)
%             if stim_line(m)>thresh(m)
%                 stim_effect_latencies(k)=m/fs;
%                 break;
%             end
%             m=m+1;
%         end
    end
    num_unique_trains=length(unique(all_stim_trains));
    catch_comparison_index_sorted=cell(1,num_unique_trains);
    stim_aligned_single_catch_sub_sorted=cell(1,num_unique_trains);
    stim_phases_sorted=cell(1,num_unique_trains);
    stim_aligned_sorted=cell(1,num_unique_trains);
    stim_effect_latencies_sorted=cell(1,num_unique_trains);
    pressure_clips_stim_sorted=cell(1,num_unique_trains);
    stim_latencies_sorted=cell(1,num_unique_trains);
    for p=1:num_unique_trains
        ind_train=find(stim_trains_saved==(p-1));
        catch_comparison_index_sorted{p}=catch_comparison_index(ind_train);
        stim_aligned_single_catch_sub_sorted{p}=stim_aligned_single_catch_sub(ind_train,:);
        stim_phases_sorted{p}=stim_phases(ind_train);
        stim_aligned_sorted{p}=stim_aligned(ind_train,:);
        stim_effect_latencies_sorted{p}=stim_effect_latencies(ind_train);
        pressure_clips_stim_sorted{p}=pressure_clips_stim(ind_train,:);
        stim_latencies_sorted{p}=stim_latencies(ind_train);
    end
    
        
    
    save('Pressure_stim_different_trains_phase_aligned.mat','pressure_clips_catch','num_unique_trains','all_stim_trains',...
        'stim_trains_saved','catch_comparison_index','stim_aligned_single_catch_sub','stim_phases','time2','stim_aligned','stim_effect_latencies',...
        'pressure_clips_stim','stim_latencies','fs','pre_time','post_onset_time','pre_time_stim','post_onset_time_stim','catch_comparison_index_sorted',...
        'stim_aligned_single_catch_sub_sorted','stim_phases_sorted','stim_aligned_sorted','stim_effect_latencies_sorted','pressure_clips_stim_sorted',...
        'stim_latencies_sorted','-v7.3')
else
    load('Pressure_stim_different_trains_phase_aligned.mat','stim_aligned_single_catch_sub_sorted','num_unique_trains','time2','stim_latencies_sorted')
end

plot_colors=zeros(num_unique_trains,3);
a=jet;
plot_colors(1,:)=a(1,:);
plot_colors(end,:)=a(end,:);
for j=2:num_unique_trains-1
    plot_colors(j,:)=a(fix(1+(j-1)*(length(a)-1)/(num_unique_trains-1)),:);
end



% stim_latencies=stim_latencies(stim_latencies>-.206 & stim_latencies<.206);
% stim_effect_latencies=stim_effect_latencies(stim_latencies>-.206 & stim_latencies<.206);
% figure(2)
% plot(stim_latencies,stim_effect_latencies,'bo')
% [rho,pval]=corr(stim_latencies,stim_effect_latencies)
% 
% period=0.412;
% nbins=12;
% stim_latencies_radians=zeros(size(stim_latencies));
% for l=1:length(stim_latencies_radians)
%     if stim_latencies(l)<0
%         stim_latencies_radians(l)=stim_latencies(l)*2*pi/period+2*pi;
%         
%     elseif stim_latencies(l)>period
%         stim_latencies_radians(l)=stim_latencies(l)*2*pi/period-2*pi;
%         
%     else
%         stim_latencies_radians(l)=stim_latencies(l)*2*pi/period;
%         
%     end
% end
% 
% bin_means=zeros(nbins,1);
% bin_centers=zeros(nbins,1);
% for j=1:nbins
%     ind=(stim_latencies_radians>=(j-1)*(2*pi)/nbins & stim_latencies_radians<j*(2*pi)/nbins);
%     bin_means(j)=mean(stim_effect_latencies(ind));
%     bin_centers(j)=(j-0.5)*2*pi/nbins;
% end
% 
% figure(3)
% polar(stim_latencies_radians,stim_effect_latencies,'ro')
% hold on
% polar(bin_centers,bin_means,'bo')
% 
% [rho pval]=circ_corrcl(stim_latencies_radians',stim_effect_latencies)


% pressure_clips_catch=pressure_clips_catch_w_prestim;
if 0
    for j=1:num_unique_trains
        id_delay=find(stim_latencies_sorted{j}>-0.05 & stim_latencies_sorted{j}<0);
        stim_latencies_sorted{j}=stim_latencies_sorted{j}(id_delay);
        stim_aligned_single_catch_sub_sorted{j}=stim_aligned_single_catch_sub_sorted{j}(id_delay,:);
    end
end
% time=-pre_time:1/fs:post_onset_time;
% figure; hold on
% ax(1)=subplot(3,1,1);
% plot(time,nanmean(pressure_clips_catch),'b-',time,nanmean(pressure_clips_catch)+2*nanstd(pressure_clips_catch),'b--',time,nanmean(pressure_clips_catch)-2*nanstd(pressure_clips_catch),'b--')
% hold on
% plot(time,mean(pressure_clips_stim),'r-')
% ylabel('Mean Pressure (kPa)')
% ax(2)=subplot(3,1,3);
% plot(time,pressure_clips_catch,'b-')
% hold on
% for j=1:length(stim_latencies)
%     plot(time(fix((pre_time+stim_latencies(j))*fs):min(fix((pre_time+stim_latencies(j)+0.1)*fs),size(pressure_clips_stim,2))),...
%         pressure_clips_stim(j,fix((pre_time+stim_latencies(j))*fs):min(fix((pre_time+stim_latencies(j)+0.1)*fs),size(pressure_clips_stim,2))),'r-')
% %     plot(time(fix((pre_time+stim_latencies(j))*fs):end),...
% %         pressure_clips_stim(j,fix((pre_time+stim_latencies(j))*fs):end),'r-')
%     hold on
% end
% ylabel('Individual File Pressures (kPa)')
% xlabel('Time Relative to Expiration Onset (s)')
% 
% ax(3)=subplot(3,1,2);
% plot(time,zeros(size(time)),'b-',time,2*nanstd(pressure_clips_catch),'b--',time,2*nanstd(pressure_clips_catch),'b--')
% hold on
% plot(time,mean(pressure_clips_stim)-nanmean(pressure_clips_catch),'r-')
% ylabel('Mean-Catch-Subtracted Mean Pressure (kPa)')

% ax(4)=subplot(2,2,4);
% plot(time,pressure_clips_catch-repmat(nanmean(pressure_clips_catch),size(pressure_clips_catch,1),1),'b-')
% hold on
% mean_pressure_clips_catch=nanmean(pressure_clips_catch);
% for j=1:length(stim_latencies)
%     plot(time(fix((pre_time+stim_latencies(j))*fs):min(fix((pre_time+stim_latencies(j)+0.1)*fs),size(pressure_clips_stim,2))),...
%         pressure_clips_stim(j,fix((pre_time+stim_latencies(j))*fs):min(fix((pre_time+stim_latencies(j)+0.1)*fs),size(pressure_clips_stim,2)))...
%         -mean_pressure_clips_catch(fix((pre_time+stim_latencies(j))*fs):min(fix((pre_time+stim_latencies(j)+0.1)*fs),size(pressure_clips_stim,2))),'r-')
% %     plot(time(fix((pre_time+stim_latencies(j))*fs):end),...
% %         pressure_clips_stim(j,fix((pre_time+stim_latencies(j))*fs):end),'r-')
%     hold on
% end
% ylabel('Mean-Catch-Subtracted Individual File Pressures (kPa)')
% xlabel('Time Relative to Expiration Onset (s)')
% % linkaxes(ax,'x')
% fs=20000;
% time2=time2/30000*fs;
% tot_time=.11;
% no_stim_mean=nanmean(stim_aligned_single_catch_sub_sorted{end-1});
h=figure(1);
% pulses=2;
if ~isempty(strfind(stimtrain_csv,'3pulses_20ms'))
    pulse_intervals=[10 10;12 8;14 6;16 4;18 2;2 18;4 16;6 14;8 12];
elseif ~isempty(strfind(stimtrain_csv,'3pulses_40ms'))
    pulse_intervals=[12 28;16 24;20 20;24 16;28 12;32 8;36 4;4 36;8 32];
elseif ~isempty(strfind(stimtrain_csv,'curare'))||~isempty(strfind(stimtrain_csv,'saline'))
    pulse_intervals=[10 10;2 0;20 0];
elseif ~isempty(strfind(stimtrain_csv,'3pulses'))
%     pulse_intervals=[10 10;11 9;12 8;13 7;14 6;15 5;16 4;4 16;5 15;6 14;7 13;8 12;9 11];
    pulse_intervals=[10 10;12 8;14 6;16 4;18 2;2 18;4 16;6 14;8 12];
%     pulse_intervals=[1 19;10 10;11 9;12 8;13 7;14 6;15 5;16 4;17 3;18 2;19 1;2 18;20 0;3 17;4 16;5 15;6 14;7 13;8 12;9 11];
%     pulse_intervals=[10 10;11 9;12 8;13 7;14 6;15 5;16 4;17 3;18 2;19 1;1 19;20 0;2 18;3 17;4 16;5 15;6 14;7 13;8 12;9 11];
%     pulse_intervals=[12 28;16 24;20 20;24 16;28 12;32 8;36 4;4 36;8 32];
elseif ~isempty(strfind(stimtrain_csv,'2pulses_subset'))
        pulse_intervals=[2;40];
elseif ~isempty(strfind(stimtrain_csv,'2pulses'))
%     pulse_intervals=[10;12;14;16;18;20;25;30;4;40;5;6;60;7;8;9];
    pulse_intervals=[1;10;12;14;16;18;2;20;25;3;30;4;40;5;6;60;7;8;9];
%     pulse_intervals=[1;10;14;18;2;20;4;40;6;60;8];
%     pulse_intervals=[10;2;20;4;40;6;60;8];
end
[y,i]=sort(pulse_intervals(:,1));
for k=1:length(i)
    stim_aligned_single_catch_sub_sorted_temp{k}=stim_aligned_single_catch_sub_sorted{i(k)};
end
stim_aligned_single_catch_sub_sorted_temp{length(i)+1}=stim_aligned_single_catch_sub_sorted{length(i)+1};
stim_aligned_single_catch_sub_sorted_temp{length(i)+2}=stim_aligned_single_catch_sub_sorted{length(i)+2};
stim_aligned_single_catch_sub_sorted=stim_aligned_single_catch_sub_sorted_temp;
nostim_mean=nanmean(stim_aligned_single_catch_sub_sorted{end-1});
% nostim_mean=zeros(1,length(time2));
for j=1:num_unique_trains
    
    % time2=0:1/fs:0.050;
    % plot(time2,zeros(1,0.050*fs+1),'b-')
    % plot(time2,zeros(1,tot_time*fs+1),'b-')
    if j==1
        plot(time2,zeros(1,length(time2)),'b-')
%         hold on
%         plot(time2,nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch)))/sqrt(sum(sum(~isnan(pressure_clips_catch)))),'b--')
%         plot(time2,-nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch)))/sqrt(sum(sum(~isnan(pressure_clips_catch)))),'b--')
    end
    hold on
     

    % plot(time2,nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch))),'b--')
    % plot(time2,-nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch))),'b--')
   
    % stim_aligned=zeros(length(stim_latencies),length(time2));
    % for m=1:length(stim_latencies)
    %     temp=pressure_clips_stim(m,:)-nanmean(pressure_clips_catch);
    %     temp2=temp(fix((pre_time+stim_latencies(m))*fs):min(fix((pre_time+stim_latencies(m))*fs)+0.05*fs,length(temp)));
    %     stim_aligned(m,1:length(temp2))=temp2;
    %     disp(m);
    % %     plot(time2,stim_aligned(m,:),'r-')
    % end
%     plot(time2,nanmean(stim_aligned_sorted{j}),'color',plot_colors(j,:),'LineWidth',2,'LineStyle','-')
    temp_mean=nanmean(stim_aligned_single_catch_sub_sorted{j});
%     plot(time2,temp_mean-temp_mean(1),'color',plot_colors(j,:),'LineWidth',2,'LineStyle','-')
%     plot(time2,temp_mean-nostim_mean,'color',plot_colors(j,:),'LineWidth',2,'LineStyle','-')
%     plot(time2,nanmean(stim_aligned_single_catch_sub_sorted{j})-no_stim_mean,'color',plot_colors(j,:),'LineWidth',2,'LineStyle','-')
    hold on
    % plot(time2,nanmean(stim_aligned)+nanstd(stim_aligned),'k--',time2,nanmean(stim_aligned)-nanstd(stim_aligned),'k--')
%     plot(time2,nanmean(stim_aligned_sorted{j})+nanstd(stim_aligned_sorted{j})/sqrt(size(stim_aligned_sorted{j},2)),'color',plot_colors(j,:),'LineStyle','--',time2,nanmean(stim_aligned_sorted{j})-nanstd(stim_aligned_sorted{j})/sqrt(size(stim_aligned_sorted{j},2)),'color',plot_colors(j,:),'LineStyle','--')
%     plot(time2,nanmean(stim_aligned_single_catch_sub_sorted{j})+nanstd(stim_aligned_single_catch_sub_sorted{j})/sqrt(size(stim_aligned_single_catch_sub_sorted{j},2)),'color',plot_colors(j,:),'LineStyle',':',...
%         time2,nanmean(stim_aligned_single_catch_sub_sorted{j})-nanstd(stim_aligned_single_catch_sub_sorted{j})/sqrt(size(stim_aligned_single_catch_sub_sorted{j},2)),'color',plot_colors(j,:),'LineStyle',':')
end
xlabel('Aligned time from Stim onset')
ylabel('Pressure (kPa)')
title('mean and stand err for aligned stims')
if ~isempty(strfind(stimtrain_csv,'2pulses_subset'))
    legend('Baseline','2','40','nostim','singlepulse')
elseif strcmp(stimtrain_csv,'bl48wh78_250uA_3pulses_200815_1146.stimtrains.csv')
    legend('Baseline','1 19','2 18','3 17','4 16','5 15','6 14','7 13','8 12','9 11','10 10','11 9','12 8','13 7','14 6','15 5','16 4','17 3','18 2','19 1','20','nostim','singlepulse');
elseif strcmp(stimtrain_csv,'gr73or13_250uA_3pulses_041215_1253.stimtrains.csv')||strcmp(stimtrain_csv,'gy13wh6_100uA_3pulses_241115_1156.stimtrains.csv')||strcmp(stimtrain_csv,'gy13wh6_50uA_3pulses_241115_0916.stimtrains.csv')||strcmp(stimtrain_csv,'bl48wh78_250uA_3pulses_100ms_latency_200815_1552.stimtrains.csv')||strcmp(stimtrain_csv,'bl69wh93_3pulses_20ms_250uA_030915_1351.stimtrains.csv')||strcmp(stimtrain_csv,'rd48wh18_250uA_3pulses_20ms_090915_1038.stimtrains.csv')||strcmp(stimtrain_csv,'gr60yw60_250uA_3pulses_20ms_150915_1026.stimtrains.csv')
    legend('Baseline','2 18','4 16','6 14','8 12','10 10','12 8','14 6','16 4','18 2','nostim','singlepulse');
elseif strcmp(stimtrain_csv,'lb61rd26_3pulses_40ms_250uA_020915_1103.stimtrains.csv')||strcmp(stimtrain_csv,'bl69wh93_3pulses_40ms_250uA_030915_1056.stimtrains.csv')||strcmp(stimtrain_csv,'rd48wh18_250uA_3pulses_40ms_090915_1208.stimtrains.csv')||strcmp(stimtrain_csv,'gr60yw60_250uA_3pulses_40ms_150915_1118.stimtrains.csv')
    legend('Baseline','4 36','8 32','12 28','16 24','20 20','24 16','28 12','32 8','36 4','nostim','singlepulse');
elseif strcmp(stimtrain_csv,'bl48wh78_250uA_2pulses_100ms_latency_200815_1658.stimtrains.csv')
    legend('Baseline','2','4','6','8','10','20','40','60','nostim','singlepulse');
elseif ~isempty(strfind(stimtrain_csv,'curare'))|| ~isempty(strfind(stimtrain_csv,'saline'))
    legend('Baseline','2','10 10','20','nostim','singlepulse');
elseif strcmp(stimtrain_csv,'bl48wh78_250uA_2pulses_200815_1332.stimtrains.csv')||strcmp(stimtrain_csv,'bl69wh93_2pulses_250uA_030915_1152.stimtrains.csv')
    legend('Baseline','1','2','3','4','5','6','7','8','9','10','12','14','16','18','20','25','30','40','60','nostim','singlepulse');
elseif strcmp(stimtrain_csv,'gy13wh6_100uA_2pulses_241115_1307.stimtrains.csv')||strcmp(stimtrain_csv,'gy13wh6_50uA_2pulses_241115_1027.stimtrains.csv')||strcmp(stimtrain_csv,'rd48wh18_250uA_2pulses_090915_1623.stimtrains.csv')||strcmp(stimtrain_csv,'rd48wh18_250uA_2pulses_090915_1338.stimtrains.csv')||strcmp(stimtrain_csv,'gr60yw60_250uA_2pulses_150915_1226.stimtrains.csv')
    legend('Baseline','1','2','4','6','8','10','14','18','20','40','60','nostim','singlepulse');
end
% if ~isempty(strfind(stimtrain_csv,'3pulses'))
% %     legend('Baseline','10 10','11 9','13 7','15 5','17 3','3 17','5 15','7 13','9 11');
% %     legend('Baseline','10 10','11 9','12 8','13 7','14 6','15 5','16 4','4 16','5 15','6 14','7 13','8 12','9 11','single pulse');
%     legend('Baseline','10 10','12 8','14 6','16 4','18 2','2 18','4 16','6 14','8 12','nostim','singlepulse');
% elseif ~isempty(strfind(stimtrain_csv,'2pulses'))
%     legend('Baseline','10','20','3','4','40','5','6','60','8');
% elseif ~isempty(strfind(stimtrain_csv,'1pulse'))
%     legend('Baseline','Single Pulse')
% end

for j=1:num_unique_trains
% for j=5:6
%     if j==1
%         hold on
% %         plot(time2,nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch)))/sqrt(sum(sum(~isnan(pressure_clips_catch)))),'b--')
% %         plot(time2,-nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch)))/sqrt(sum(sum(~isnan(pressure_clips_catch)))),'b--')
%         temp=nanmean(pressure_clips_catch);
%         mean_sub_catch_clips=zeros(size(pressure_clips_catch));
%         for k=1:size(pressure_clips_catch)
%             mean_sub_catch_clips(k,:)=pressure_clips_catch(k,:)-temp;
%         end
% %         plot(time2,nanstd(reshape(mean_sub_catch_clips,1,numel(mean_sub_catch_clips))),'b--')
% %         plot(time2,-nanstd(reshape(mean_sub_catch_clips,1,numel(mean_sub_catch_clips))),'b--')
%         plot(time2,nanstd(reshape(mean_sub_catch_clips,1,numel(mean_sub_catch_clips)))/sqrt(sum(sum(~isnan(mean_sub_catch_clips)))),'b--')
%         plot(time2,-nanstd(reshape(mean_sub_catch_clips,1,numel(mean_sub_catch_clips)))/sqrt(sum(sum(~isnan(mean_sub_catch_clips)))),'b--')
%     end
    temp_mean=nanmean(stim_aligned_single_catch_sub_sorted{j});
    hold on
%     plot(time2,nanmean(stim_aligned_sorted{j})+nanstd(stim_aligned_sorted{j})/sqrt(size(stim_aligned_sorted{j},2)),'color',plot_colors(j,:),'LineStyle',':')
%     plot(time2,nanmean(stim_aligned_sorted{j})-nanstd(stim_aligned_sorted{j})/sqrt(size(stim_aligned_sorted{j},2)),'color',plot_colors(j,:),'LineStyle',':')
%     plot(time2,temp_mean-temp_mean(1)+nanstd(stim_aligned_single_catch_sub_sorted{j})/sqrt(size(stim_aligned_single_catch_sub_sorted{j},1)),'color',plot_colors(j,:),'LineStyle','--')
%     plot(time2,temp_mean-temp_mean(1)-nanstd(stim_aligned_single_catch_sub_sorted{j})/sqrt(size(stim_aligned_single_catch_sub_sorted{j},1)),'color',plot_colors(j,:),'LineStyle','--')

%     plot(time2,temp_mean-nostim_mean+nanstd(stim_aligned_single_catch_sub_sorted{j})/sqrt(size(stim_aligned_single_catch_sub_sorted{j},1)),'color',plot_colors(j,:),'LineStyle','--')
%     plot(time2,temp_mean-nostim_mean-nanstd(stim_aligned_single_catch_sub_sorted{j})/sqrt(size(stim_aligned_single_catch_sub_sorted{j},1)),'color',plot_colors(j,:),'LineStyle','--')
    X=[time2 fliplr(time2)];
    y1=temp_mean-nostim_mean-nanstd(stim_aligned_single_catch_sub_sorted{j})/sqrt(size(stim_aligned_single_catch_sub_sorted{j},1));
    y2=temp_mean-nostim_mean+nanstd(stim_aligned_single_catch_sub_sorted{j})/sqrt(size(stim_aligned_single_catch_sub_sorted{j},1));
    Y=[y1 fliplr(y2)];
    p=fill(X,Y,plot_colors(j,:));
    set(p,'EdgeColor','none');
%     plot(time2,temp_mean-nostim_mean,'k')
%     plot(time2,nanmean(stim_aligned_single_catch_sub_sorted{j})-no_stim_mean+nanstd(stim_aligned_single_catch_sub_sorted{j})/sqrt(size(stim_aligned_single_catch_sub_sorted{j},1)),'color',plot_colors(j,:),'LineStyle','--')
%     plot(time2,nanmean(stim_aligned_single_catch_sub_sorted{j})-no_stim_mean-nanstd(stim_aligned_single_catch_sub_sorted{j})/sqrt(size(stim_aligned_single_catch_sub_sorted{j},1)),'color',plot_colors(j,:),'LineStyle','--')
end

    
% % figure(2)
% % for j=1:length(stim_latencies)
% %     
% % end
% stim_phases=1000*stim_phases;
% figure
% ax1(1)=subplot(2,1,1);
% plot(stim_phases,nanmean(stim_aligned(:,0.01*fs+1:end),2),'bo')
% xlabel('Stim phase (radians)')
% ylabel('Mean value over 50ms stimulation period')
% 
% ax1(2)=subplot(2,1,2);
% plot(stim_phases,nanmax(stim_aligned(:,0.01*fs+1:end)'),'ro')
% xlabel('Stim phase (radians)')
% ylabel('Max value over 50ms stimulation period')
% linkaxes(ax1,'x');
% 
% 
% nbins=7;
% bin_size=pi/7;
% stim_aligned_divided=cell(1,nbins);
% for h=1:length(stim_phases)
%     if ceil(stim_phases(h)/bin_size)<nbins+1 && ceil(stim_phases(h)/bin_size)>0
%         stim_aligned_divided{ceil(stim_phases(h)/bin_size)}=[stim_aligned_divided{ceil(stim_phases(h)/bin_size)}; stim_aligned(h,:)];
%     end
% end
% 
% figure
% % plot(time2,zeros(1,0.050*fs+1),'b-')
% 
% hold on
% 
% 
% col_str='rgbkymc';
% for m=1:nbins
%     plot(time2,nanmean(stim_aligned_divided{m}),col_str(m))
%     hold on
% end
% legend('1st pi/7','2nd pi/7','3rd pi/7','4th pi/7','5th pi/7','6th pi/7','7th pi/7')
% xlabel('Time relative to stimulation (s)')
% ylabel('Pressure (kPa)')
% % plot(time2,nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch))),'b--')
% % plot(time2,-nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch))),'b--')
% plot(time2,nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch)))/sqrt(sum(sum(~isnan(pressure_clips_catch)))),'b--')
% plot(time2,-nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch)))/sqrt(sum(sum(~isnan(pressure_clips_catch)))),'b--')
% hold on
% for m=1:nbins
%     plot(time2,nanmean(stim_aligned_divided{m})+nanstd(stim_aligned_divided{m})/sqrt(size(stim_aligned_divided{m},1)),[col_str(m) '--'],...
%         time2,nanmean(stim_aligned_divided{m})-nanstd(stim_aligned_divided{m})/sqrt(size(stim_aligned_divided{m},1)),[col_str(m) '--'])
% end
% % plot(time2,zeros(1,tot_time*fs+1),'b-')
% plot(time2,zeros(1,length(time2)),'b-')
% 

% figure(100)
% hold on
% % plot(time2,zeros(1,tot_time*fs+1),'b-')
% plot(time2,zeros(1,length(time2)),'b-')
% hold on
% % % plot(time2,nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch))),'b--')
% % % plot(time2,-nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch))),'b--')
% % plot(time2,nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch)))/sqrt(sum(sum(~isnan(pressure_clips_catch)))),'b--')
% % plot(time2,-nanstd(reshape(pressure_clips_catch,1,numel(pressure_clips_catch)))/sqrt(sum(sum(~isnan(pressure_clips_catch)))),'b--')
% col_str2='b';
% plot(time2,nanmean(stim_aligned_divided{1}),col_str2)
% plot(time2,nanmean(stim_aligned_divided{1})+nanstd(stim_aligned_divided{1})/sqrt(size(stim_aligned_divided{1},1)),[col_str2 '--'])
% plot(time2,nanmean(stim_aligned_divided{1})-nanstd(stim_aligned_divided{1})/sqrt(size(stim_aligned_divided{1},1)),[col_str2 '--'])
% xlabel('Time Relative to stimulation (s)')
% ylabel('Pressure (kPa)')
% singlepulsemean=nanmean(stim_aligned_single_catch_sub_sorted{1});
% singlepulsesterr=nanstd(stim_aligned_single_catch_sub_sorted{1})/sqrt(size(stim_aligned_single_catch_sub_sorted{1},1));
% save('singlepulse.mat','singlepulsemean','singlepulsesterr')

ind=strfind(stimtrain_csv,'_');
figsavename=[stimtrain_csv(1:ind(3)-1) '_phase_aligned'];
saveas(h,figsavename,'fig');
saveas(h,figsavename,'jpg');

end

function stim_times=compute_stim_times(dat_fn)
load(dat_fn)
if ~exist('stim_times','var')
    trig=trig_dat;
    stim_times=[];
    t=3/100*fs+1;
    while t<length(trig)
        if trig(t)>mean(trig(t-(3/100*fs):t-1))+50*std(trig(t-(3/100*fs):t-1))
            stim_times=[stim_times t*1000/fs];   %convert to ms to match onsets time
            t=t+80/1000*fs;  %stim is 20 ms long, so this jumps over the stim for future searching
    %         disp(['Computed stim_time = ' num2str(stim_times(end))])
        else
            t=t+1;
        end
        percent=t/length(trig)*100;
        if mod(t,5*fs)==0
            disp(['Computed Stim times for ' num2str(round(percent)) '% of file'])
        end
    end
    save(dat_fn,'stim_times','-append');
end
end