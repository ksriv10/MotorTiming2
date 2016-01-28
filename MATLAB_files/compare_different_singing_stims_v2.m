function [] = compare_different_singing_stims_v2(dirs,syl)
%use this script to compare stim results with different

% dir1='T:\Kyles_Right_Computer\bl22wh20\STIM3\9_3_IPIS';
% dir2='T:\Kyles_Right_Computer\bl22wh20\STIM3\3_9_IPIS';
% dir1='C:\DATA\useful_birds\singing_stim_different_trains\lb61rd26\STIM\500uA\9_3_IPIS';
% dir2='C:\DATA\useful_birds\singing_stim_different_trains\lb61rd26\STIM\500uA\3_9_IPIS';
% dir1='C:\DATA\useful_birds\singing_stim_different_trains\lb61rd26\STIM\2mA\9_3_IPIS\temp';
% dir2='C:\DATA\useful_birds\singing_stim_different_trains\lb61rd26\STIM\2mA\3_9_IPIS\temp';
% dir1='C:\DATA\useful_birds\singing_stim_different_trains\gr16lb166\round4\6_3_IPIS';
% dir2='C:\DATA\useful_birds\singing_stim_different_trains\gr16lb166\round4\3_6_IPIS';
% dir3='C:\DATA\useful_birds\singing_stim_different_trains\gr16lb166\round4\1pulse2';
% dir1='C:\DATA\useful_birds\singing_stim_different_trains\gr97or27\STIM\500uA\9_3_IPIS';
% dir2='C:\DATA\useful_birds\singing_stim_different_trains\gr97or27\STIM\500uA\3_9_IPIS';
% dir1='C:\DATA\useful_birds\singing_stim_different_trains\gr97or27\STIM\2mA_pt1\9_3_IPIS';
% dir2='C:\DATA\useful_birds\singing_stim_different_trains\gr97or27\STIM\2mA_pt1\3_9_IPIS';
% dir3='C:\DATA\useful_birds\singing_stim_different_trains\gr97or27\STIM\2mA_pt1\1pulse2';
% dir1='C:\DATA\useful_birds\singing_stim_different_trains\gr97or27\STIM\2mA_pt2\9_3_IPIS\for_analysis';
% dir2='C:\DATA\useful_birds\singing_stim_different_trains\gr97or27\STIM\2mA_pt2\3_9_IPIS\for_analysis';
% dir3='C:\DATA\useful_birds\singing_stim_different_trains\gr97or27\STIM\2mA_pt2\1pulse2_39\for_analysis';


% dir1='C:\DATA\useful_birds\singing_stim_different_trains\or41wh181\STIM\6_3_IPIS';
% dir2='C:\DATA\useful_birds\singing_stim_different_trains\or41wh181\STIM\3_6_IPIS';
% dir3='C:\DATA\useful_birds\singing_stim_different_trains\or41wh181\STIM\1pulse2_36';
% dir1='C:\DATA\useful_birds\singing_stim_different_trains\gr49or99\STIM\6_3_IPIS';
% dir2='C:\DATA\useful_birds\singing_stim_different_trains\gr49or99\STIM\3_6_IPIS';
% dir3='C:\DATA\useful_birds\singing_stim_different_trains\gr49or99\STIM\1pulse2_36';
% dir1='C:\DATA\useful_birds\singing_stim_different_trains\or58wh3\STIM2\6_3_IPIS';
% dir2='C:\DATA\useful_birds\singing_stim_different_trains\or58wh3\STIM2\3_6_IPIS';
% dir3='C:\DATA\useful_birds\singing_stim_different_trains\or58wh3\STIM2\1pulse2_36';
% dir1='C:\DATA\useful_birds\singing_stim_different_trains\or48wh188\STIM2\6_3_IPIS';
% dir2='C:\DATA\useful_birds\singing_stim_different_trains\or48wh188\STIM2\3_6_IPIS';
% dir3='C:\DATA\useful_birds\singing_stim_different_trains\or48wh188\STIM2\1pulse2_36';
% dir1='C:\DATA\useful_birds\singing_stim_different_trains\gr55pu56_LDTB\STIM\9_3_IPIS';
% dir2='C:\DATA\useful_birds\singing_stim_different_trains\gr55pu56_LDTB\STIM\3_9_IPIS';
% dir3='C:\DATA\useful_birds\singing_stim_different_trains\gr55pu56_LDTB\STIM\1pulse2_39';


cd(dirs{1}); 
% syl='c';
filename=['TEMP_OUTPUT_SE_256_PC_' syl '.mat'];
id=strfind(dirs{1},'\');
birdid=dirs{1}((id(4)+1):(id(5)-1));
if ~isempty(strfind(birdid,'LDTB'))
    id=strfind(birdid,'LDTB');
    birdid=birdid(1:id-2);
end
if isempty(strfind(dirs{1},'9_3'))
    savename=[birdid '_3663_' filename(end-4:end)];
    legend_str=sprintf('legend(''Catch'',''6 3'',''3 6'',''1 pulse'',''Stim time'')');

else
    savename=[birdid '_3993_' filename(end-4:end)];
    legend_str=sprintf('legend(''Catch'',''9 3'',''3 9'',''1 pulse'',''Stim time'')');
end
load(filename)
fs=32000;


if strcmp(birdid,'gr16lb166')
    id_delay2=1:1:50;
else
    id_delay2=1:1:1000;
end

if 1
    %     id_delay2=1:1:(length(amp_stimmed_syl)/2);
    %     id_delay2=1:1:75;

    if strcmp(birdid,'gr16lb166')
%         id_delay=find(delay_from_stim_to_onset_of_stim_syl>(0-pq_offset) & delay_from_stim_to_onset_of_stim_syl<(50-pq_offset)); % specify delay in msec
        id_delay=find(delay_from_stim_to_onset_of_stim_syl<(-60) & delay_from_stim_to_onset_of_stim_syl>(-95)); % specify delay in msec
    elseif strcmp(birdid,'gr97or27')
        if syl=='x'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-30)); % specify delay in msec
        elseif syl=='y'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-30)); % specify delay in msec
        end
    elseif strcmp(birdid,'or41wh181')
        if syl=='a'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-20)); % specify delay in msec
        elseif syl=='b'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        end 
    elseif strcmp(birdid,'gr49or99')
        id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
    elseif strcmp(birdid,'or58wh3')
        if syl=='a'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-25)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        end
    elseif strcmp(birdid,'or48wh188')
        id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-30)); % specify delay in msec
    elseif strcmp(birdid,'gr55pu56')
        if syl=='a'
%             id_delay=find(delay_from_stim_to_onset_of_stim_syl>(-50-pq_offset) & delay_from_stim_to_onset_of_stim_syl<(50-pq_offset)); % specify delay in msec
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-100)); % specify delay in msec
        elseif syl=='f'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-40)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-40)); % specify delay in msec
        end
    end
        id_delay=intersect(id_delay,id_delay2);
    %     id_delay3=find(delay_from_stim_to_onset_of_stim_syl>median(delay_from_stim_to_onset_of_stim_syl));
    %     id_delay=intersect(id_delay3,id_delay);
        amp_stimmed_syl=amp_stimmed_syl(id_delay);
        amp_rms_stimmed_syl=amp_rms_stimmed_syl(id_delay);
        stim_clips_smoothed=stim_clips_smoothed(id_delay,:);
        pitch_stimmed_syl=pitch_stimmed_syl(id_delay);
        se_stimmed_syl=se_stimmed_syl(id_delay);
        P1_allfreqs_stim=P1_allfreqs_stim(id_delay,:);
        P1_stim=P1_stim(id_delay,:);
        F1_allfreqs_stim=F1_allfreqs_stim(id_delay,:);
        F1_stim=F1_stim(id_delay,:);
        tmp=cell(0);
        tmp2=[];
        tmp3=[];
        for x=1:length(id_delay)
            tmp{x}=stim_clips{id_delay(x)};
            tmp2(:,x)=pitch_contours_stim(:,id_delay(x));
            tmp3(:,x)=se_contours_stim(:,id_delay(x));
        end
        stim_clips=tmp;
        pitch_contours_stim=tmp2;
        se_contours_stim=tmp3;
        delay_from_stim_to_onset_of_stim_syl=delay_from_stim_to_onset_of_stim_syl(id_delay);
end

if 1
    id_good_amp_catch=find(amp_catch_syl~=0);
    id_good_amp_stim=find(amp_stimmed_syl~=0);

%     id_good_amp_catch=intersect(id_good_amp_catch,1:100);
%     id_good_amp_stim=intersect(id_good_amp_stim,1:100);
    amp_stimmed_syl=amp_stimmed_syl(id_good_amp_stim);
    amp_rms_stimmed_syl=amp_rms_stimmed_syl(id_good_amp_stim);
    catch_clips_smoothed=catch_clips_smoothed(id_good_amp_catch,:);
    stim_clips_smoothed=stim_clips_smoothed(id_good_amp_stim,:);
    pitch_stimmed_syl=pitch_stimmed_syl(id_good_amp_stim);
    se_stimmed_syl=se_stimmed_syl(id_good_amp_stim);
    P1_allfreqs_stim=P1_allfreqs_stim(id_good_amp_stim,:);
    P1_stim=P1_stim(id_good_amp_stim,:);
    F1_allfreqs_stim=F1_allfreqs_stim(id_good_amp_stim,:);
    F1_stim=F1_stim(id_good_amp_stim,:);
    tmp=cell(0);
    tmp2=[];
    tmp3=[];
    for x=1:length(id_good_amp_stim)
        tmp{x}=stim_clips{id_good_amp_stim(x)};
        tmp2(:,x)=pitch_contours_stim(:,id_good_amp_stim(x));
        tmp3(:,x)=se_contours_stim(:,id_good_amp_stim(x));
    end
    stim_clips=tmp;
    pitch_contours_stim=tmp2;
    se_contours_stim=tmp3;
    delay_from_stim_to_onset_of_stim_syl=delay_from_stim_to_onset_of_stim_syl(id_good_amp_stim);
    
    amp_catch_syl=amp_catch_syl(id_good_amp_catch);
    amp_rms_catch_syl=amp_rms_catch_syl(id_good_amp_catch);
    pitch_catch_syl=pitch_catch_syl(id_good_amp_catch);
    se_catch_syl=se_catch_syl(id_good_amp_catch);
    P1_allfreqs_catch=P1_allfreqs_catch(id_good_amp_catch,:);
    P1_catch=P1_catch(id_good_amp_catch,:);
    F1_allfreqs_catch=F1_allfreqs_catch(id_good_amp_catch,:);
    F1_catch=F1_catch(id_good_amp_catch,:);
    tmp=cell(0);
    tmp2=[];
    tmp3=[];
    for x=1:length(id_good_amp_catch)
        tmp{x}=catch_clips{id_good_amp_catch(x)};
        tmp2(:,x)=pitch_contours_catch(:,id_good_amp_catch(x));
        tmp3(:,x)=se_contours_catch(:,id_good_amp_catch(x));
    end
    catch_clips=tmp;
    pitch_contours_catch=tmp2;
    se_contours_catch=tmp3;
end
if 1
    amp_catch_syl=amp_rms_catch_syl;
    amp_stimmed_syl=amp_rms_stimmed_syl;
end
disp('Converting amplitudes to log amplitudes')
amp_catch_syl=20*log10(amp_catch_syl/pref)+94;
amp_stimmed_syl=20*log10(amp_stimmed_syl/pref)+94;

amp_catch_syl1=amp_catch_syl;
amp_stimmed_syl1=amp_stimmed_syl;
se_catch_syl1=se_catch_syl;
se_stimmed_syl1=se_stimmed_syl;
pitch_catch_syl1=pitch_catch_syl;
pitch_stimmed_syl1=pitch_stimmed_syl;
pitch_contours_catch1=pitch_contours_catch;
pitch_contours_stim1=pitch_contours_stim;
stim_clips_smoothed1=stim_clips_smoothed;
catch_clips_smoothed1=catch_clips_smoothed;

se_contours_catch1=se_contours_catch;
se_contours_stim1=se_contours_stim;
delay_from_stim_to_onset_of_stim_syl1=delay_from_stim_to_onset_of_stim_syl;

cd(dirs{2}); 
load(filename)

if 1
%     id_delay2=1:1:(length(amp_stimmed_syl)/2);
%     id_delay2=1:1:75;   
    if strcmp(birdid,'gr16lb166')
%         id_delay=find(delay_from_stim_to_onset_of_stim_syl>(0-pq_offset) & delay_from_stim_to_onset_of_stim_syl<(50-pq_offset)); % specify delay in msec
        id_delay=find(delay_from_stim_to_onset_of_stim_syl<(-60) & delay_from_stim_to_onset_of_stim_syl>(-95)); % specify delay in msec
    elseif strcmp(birdid,'gr97or27')
        if syl=='x'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-30)); % specify delay in msec
        elseif syl=='y'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-30)); % specify delay in msec
        end
    elseif strcmp(birdid,'or41wh181')
        if syl=='a'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-20)); % specify delay in msec
        elseif syl=='b'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        end 
    elseif strcmp(birdid,'gr49or99')
        id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
    elseif strcmp(birdid,'or58wh3')
        if syl=='a'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-25)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        end
    elseif strcmp(birdid,'or48wh188')
        id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-30)); % specify delay in msec
    elseif strcmp(birdid,'gr55pu56')
        if syl=='a'
%             id_delay=find(delay_from_stim_to_onset_of_stim_syl>(-50-pq_offset) & delay_from_stim_to_onset_of_stim_syl<(50-pq_offset)); % specify delay in msec
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-100)); % specify delay in msec
        elseif syl=='f'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-40)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-40)); % specify delay in msec
        end
    end
    id_delay=intersect(id_delay,id_delay2);
%     id_delay3=find(delay_from_stim_to_onset_of_stim_syl>median(delay_from_stim_to_onset_of_stim_syl));
%     id_delay=intersect(id_delay3,id_delay);
    amp_stimmed_syl=amp_stimmed_syl(id_delay);
    amp_rms_stimmed_syl=amp_rms_stimmed_syl(id_delay);
    stim_clips_smoothed=stim_clips_smoothed(id_delay,:);
    pitch_stimmed_syl=pitch_stimmed_syl(id_delay);
    se_stimmed_syl=se_stimmed_syl(id_delay);
    P1_allfreqs_stim=P1_allfreqs_stim(id_delay,:);
    P1_stim=P1_stim(id_delay,:);
    F1_allfreqs_stim=F1_allfreqs_stim(id_delay,:);
    F1_stim=F1_stim(id_delay,:);
    tmp=cell(0);
    tmp2=[];
    tmp3=[];
    for x=1:length(id_delay)
        tmp{x}=stim_clips{id_delay(x)};
        tmp2(:,x)=pitch_contours_stim(:,id_delay(x));
        tmp3(:,x)=se_contours_stim(:,id_delay(x));
    end
    stim_clips=tmp;
    pitch_contours_stim=tmp2;
    se_contours_stim=tmp3;
    delay_from_stim_to_onset_of_stim_syl=delay_from_stim_to_onset_of_stim_syl(id_delay);
end

if 1
    id_good_amp_catch=find(amp_catch_syl~=0);
    id_good_amp_stim=find(amp_stimmed_syl~=0);

%     id_good_amp_catch=intersect(id_good_amp_catch,1:100);
%     id_good_amp_stim=intersect(id_good_amp_stim,1:100);
    amp_stimmed_syl=amp_stimmed_syl(id_good_amp_stim);
    amp_rms_stimmed_syl=amp_rms_stimmed_syl(id_good_amp_stim);
    catch_clips_smoothed=catch_clips_smoothed(id_good_amp_catch,:);
    stim_clips_smoothed=stim_clips_smoothed(id_good_amp_stim,:);
    pitch_stimmed_syl=pitch_stimmed_syl(id_good_amp_stim);
    se_stimmed_syl=se_stimmed_syl(id_good_amp_stim);
    P1_allfreqs_stim=P1_allfreqs_stim(id_good_amp_stim,:);
    P1_stim=P1_stim(id_good_amp_stim,:);
    F1_allfreqs_stim=F1_allfreqs_stim(id_good_amp_stim,:);
    F1_stim=F1_stim(id_good_amp_stim,:);
    tmp=cell(0);
    tmp2=[];
    tmp3=[];
    for x=1:length(id_good_amp_stim)
        tmp{x}=stim_clips{id_good_amp_stim(x)};
        tmp2(:,x)=pitch_contours_stim(:,id_good_amp_stim(x));
        tmp3(:,x)=se_contours_stim(:,id_good_amp_stim(x));
    end
    stim_clips=tmp;
    pitch_contours_stim=tmp2;
    se_contours_stim=tmp3;
    delay_from_stim_to_onset_of_stim_syl=delay_from_stim_to_onset_of_stim_syl(id_good_amp_stim);
    
    amp_catch_syl=amp_catch_syl(id_good_amp_catch);
    amp_rms_catch_syl=amp_rms_catch_syl(id_good_amp_catch);
    pitch_catch_syl=pitch_catch_syl(id_good_amp_catch);
    se_catch_syl=se_catch_syl(id_good_amp_catch);
    P1_allfreqs_catch=P1_allfreqs_catch(id_good_amp_catch,:);
    P1_catch=P1_catch(id_good_amp_catch,:);
    F1_allfreqs_catch=F1_allfreqs_catch(id_good_amp_catch,:);
    F1_catch=F1_catch(id_good_amp_catch,:);
    tmp=cell(0);
    tmp2=[];
    tmp3=[];
    for x=1:length(id_good_amp_catch)
        tmp{x}=catch_clips{id_good_amp_catch(x)};
        tmp2(:,x)=pitch_contours_catch(:,id_good_amp_catch(x));
        tmp3(:,x)=se_contours_catch(:,id_good_amp_catch(x));
    end
    catch_clips=tmp;
    pitch_contours_catch=tmp2;
    se_contours_catch=tmp3;
end
if 1
    amp_catch_syl=amp_rms_catch_syl;
    amp_stimmed_syl=amp_rms_stimmed_syl;
end
disp('Converting amplitudes to log amplitudes')
amp_catch_syl=20*log10(amp_catch_syl/pref)+94;
amp_stimmed_syl=20*log10(amp_stimmed_syl/pref)+94;

amp_catch_syl2=amp_catch_syl;
amp_stimmed_syl2=amp_stimmed_syl;
se_catch_syl2=se_catch_syl;
se_stimmed_syl2=se_stimmed_syl;
pitch_catch_syl2=pitch_catch_syl;
pitch_stimmed_syl2=pitch_stimmed_syl;
pitch_contours_catch2=pitch_contours_catch;
pitch_contours_stim2=pitch_contours_stim;

stim_clips_smoothed2=stim_clips_smoothed;
catch_clips_smoothed2=catch_clips_smoothed;

se_contours_catch2=se_contours_catch;
se_contours_stim2=se_contours_stim;
delay_from_stim_to_onset_of_stim_syl2=delay_from_stim_to_onset_of_stim_syl;

cd(dirs{3}); 
load(filename)

if 1
%     id_delay2=1:1:(length(amp_stimmed_syl)/2);
%     id_delay2=1:1:75;   
    if strcmp(birdid,'gr16lb166')
%         id_delay=find(delay_from_stim_to_onset_of_stim_syl>(0-pq_offset) & delay_from_stim_to_onset_of_stim_syl<(50-pq_offset)); % specify delay in msec
        id_delay=find(delay_from_stim_to_onset_of_stim_syl<(-60) & delay_from_stim_to_onset_of_stim_syl>(-95)); % specify delay in msec
    elseif strcmp(birdid,'gr97or27')
        if syl=='x'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-30)); % specify delay in msec
        elseif syl=='y'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-30)); % specify delay in msec
        end
    elseif strcmp(birdid,'or41wh181')
        if syl=='a'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-20)); % specify delay in msec
        elseif syl=='b'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        end 
    elseif strcmp(birdid,'gr49or99')
        id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
    elseif strcmp(birdid,'or58wh3')
        if syl=='a'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-25)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-35)); % specify delay in msec
        end
    elseif strcmp(birdid,'or48wh188')
        id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-30)); % specify delay in msec
    elseif strcmp(birdid,'gr55pu56')
        if syl=='a'
%             id_delay=find(delay_from_stim_to_onset_of_stim_syl>(-50-pq_offset) & delay_from_stim_to_onset_of_stim_syl<(50-pq_offset)); % specify delay in msec
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-100)); % specify delay in msec
        elseif syl=='f'
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-40)); % specify delay in msec
        else
            id_delay=find(delay_from_stim_to_onset_of_stim_syl<(0) & delay_from_stim_to_onset_of_stim_syl>(-40)); % specify delay in msec
        end
    end
    id_delay=intersect(id_delay,id_delay2);%     id_delay3=find(delay_from_stim_to_onset_of_stim_syl>median(delay_from_stim_to_onset_of_stim_syl));
%     id_delay=intersect(id_delay3,id_delay);
    amp_stimmed_syl=amp_stimmed_syl(id_delay);
    amp_rms_stimmed_syl=amp_rms_stimmed_syl(id_delay);
    stim_clips_smoothed=stim_clips_smoothed(id_delay,:);
    pitch_stimmed_syl=pitch_stimmed_syl(id_delay);
    se_stimmed_syl=se_stimmed_syl(id_delay);
    P1_allfreqs_stim=P1_allfreqs_stim(id_delay,:);
    P1_stim=P1_stim(id_delay,:);
    F1_allfreqs_stim=F1_allfreqs_stim(id_delay,:);
    F1_stim=F1_stim(id_delay,:);
    tmp=cell(0);
    tmp2=[];
    tmp3=[];
    for x=1:length(id_delay)
        tmp{x}=stim_clips{id_delay(x)};
        tmp2(:,x)=pitch_contours_stim(:,id_delay(x));
        tmp3(:,x)=se_contours_stim(:,id_delay(x));
    end
    stim_clips=tmp;
    pitch_contours_stim=tmp2;
    se_contours_stim=tmp3;
    delay_from_stim_to_onset_of_stim_syl=delay_from_stim_to_onset_of_stim_syl(id_delay);
end

if 1
    id_good_amp_catch=find(amp_catch_syl~=0);
    id_good_amp_stim=find(amp_stimmed_syl~=0);

%     id_good_amp_catch=intersect(id_good_amp_catch,1:100);
%     id_good_amp_stim=intersect(id_good_amp_stim,1:100);
    amp_stimmed_syl=amp_stimmed_syl(id_good_amp_stim);
    amp_rms_stimmed_syl=amp_rms_stimmed_syl(id_good_amp_stim);
    catch_clips_smoothed=catch_clips_smoothed(id_good_amp_catch,:);
    stim_clips_smoothed=stim_clips_smoothed(id_good_amp_stim,:);
    pitch_stimmed_syl=pitch_stimmed_syl(id_good_amp_stim);
    se_stimmed_syl=se_stimmed_syl(id_good_amp_stim);
    P1_allfreqs_stim=P1_allfreqs_stim(id_good_amp_stim,:);
    P1_stim=P1_stim(id_good_amp_stim,:);
    F1_allfreqs_stim=F1_allfreqs_stim(id_good_amp_stim,:);
    F1_stim=F1_stim(id_good_amp_stim,:);
    tmp=cell(0);
    tmp2=[];
    tmp3=[];
    for x=1:length(id_good_amp_stim)
        tmp{x}=stim_clips{id_good_amp_stim(x)};
        tmp2(:,x)=pitch_contours_stim(:,id_good_amp_stim(x));
        tmp3(:,x)=se_contours_stim(:,id_good_amp_stim(x));
    end
    stim_clips=tmp;
    pitch_contours_stim=tmp2;
    se_contours_stim=tmp3;
    delay_from_stim_to_onset_of_stim_syl=delay_from_stim_to_onset_of_stim_syl(id_good_amp_stim);
    
    amp_catch_syl=amp_catch_syl(id_good_amp_catch);
    amp_rms_catch_syl=amp_rms_catch_syl(id_good_amp_catch);
    pitch_catch_syl=pitch_catch_syl(id_good_amp_catch);
    se_catch_syl=se_catch_syl(id_good_amp_catch);
    P1_allfreqs_catch=P1_allfreqs_catch(id_good_amp_catch,:);
    P1_catch=P1_catch(id_good_amp_catch,:);
    F1_allfreqs_catch=F1_allfreqs_catch(id_good_amp_catch,:);
    F1_catch=F1_catch(id_good_amp_catch,:);
    tmp=cell(0);
    tmp2=[];
    tmp3=[];
    for x=1:length(id_good_amp_catch)
        tmp{x}=catch_clips{id_good_amp_catch(x)};
        tmp2(:,x)=pitch_contours_catch(:,id_good_amp_catch(x));
        tmp3(:,x)=se_contours_catch(:,id_good_amp_catch(x));
    end
    catch_clips=tmp;
    pitch_contours_catch=tmp2;
    se_contours_catch=tmp3;
end
if 1
    amp_catch_syl=amp_rms_catch_syl;
    amp_stimmed_syl=amp_rms_stimmed_syl;
end
disp('Converting amplitudes to log amplitudes')
amp_catch_syl=20*log10(amp_catch_syl/pref)+94;
amp_stimmed_syl=20*log10(amp_stimmed_syl/pref)+94;

amp_catch_syl3=amp_catch_syl;
amp_stimmed_syl3=amp_stimmed_syl;
se_catch_syl3=se_catch_syl;
se_stimmed_syl3=se_stimmed_syl;
pitch_catch_syl3=pitch_catch_syl;
pitch_stimmed_syl3=pitch_stimmed_syl;
pitch_contours_catch3=pitch_contours_catch;
pitch_contours_stim3=pitch_contours_stim;

stim_clips_smoothed3=stim_clips_smoothed;
catch_clips_smoothed3=catch_clips_smoothed;

se_contours_catch3=se_contours_catch;
se_contours_stim3=se_contours_stim;
delay_from_stim_to_onset_of_stim_syl3=delay_from_stim_to_onset_of_stim_syl;

[h,p]=ttest2(amp_catch_syl1,amp_catch_syl2);
if h==1 && mean(amp_catch_syl1)>mean(amp_catch_syl2)
    fprintf('Group 1 catch amplitude is greater than Group 2 catch amplitude (p = %0.3g)\n',p);
elseif h==1 && mean(amp_catch_syl1)<mean(amp_catch_syl2)
    fprintf('Group 2 catch amplitude is greater than Group 1 catch amplitude (p = %0.3g)\n',p);
else
    fprintf('Group 1 catch amplitude is not significantly different from Group 2 catch amplitude (p = %0.3g)\n',p);
end

[h,p]=ttest2(amp_stimmed_syl1,amp_stimmed_syl2);
if h==1 && mean(amp_stimmed_syl1)>mean(amp_stimmed_syl2)
    fprintf('Group 1 stimmed amplitude is greater than Group 2 stimmed amplitude (p = %0.3g)\n',p);
elseif h==1 && mean(amp_stimmed_syl1)<mean(amp_stimmed_syl2)
    fprintf('Group 2 stimmed amplitude is greater than Group 1 stimmed amplitude (p = %0.3g)\n',p);
else
    fprintf('Group 1 stimmed amplitude is not significantly different from Group 2 stimmed amplitude (p = %0.3g)\n',p);
end

[h,p]=ttest2(se_catch_syl1,se_catch_syl2);
if h==1 && mean(se_catch_syl1)>mean(se_catch_syl2)
    fprintf('Group 1 catch SE is greater than Group 2 catch SE (p = %0.3g)\n',p);
elseif h==1 && mean(se_catch_syl1)<mean(se_catch_syl2)
    fprintf('Group 2 catch SE is greater than Group 1 catch SE (p = %0.3g)\n',p);
else
    fprintf('Group 1 catch SE is not significantly different from Group 2 catch SE (p = %0.3g)\n',p);
end

[h,p]=ttest2(se_stimmed_syl1,se_stimmed_syl2);
if h==1 && mean(se_stimmed_syl1)>mean(se_stimmed_syl2)
    fprintf('Group 1 stimmed SE is greater than Group 2 stimmed SE (p = %0.3g)\n',p);
elseif h==1 && mean(se_stimmed_syl1)<mean(se_stimmed_syl2)
    fprintf('Group 2 stimmed SE is greater than Group 1 stimmed SE (p = %0.3g)\n',p);
else
    fprintf('Group 1 stimmed SE is not significantly different from Group 2 stimmed SE (p = %0.3g)\n',p);
end

[h,p]=ttest2(pitch_catch_syl1,pitch_catch_syl2);
if h==1 && mean(pitch_catch_syl1)>mean(pitch_catch_syl2)
    fprintf('Group 1 catch pitch is greater than Group 2 catch pitch (p = %0.3g)\n',p);
elseif h==1 && mean(pitch_catch_syl1)<mean(pitch_catch_syl2)
    fprintf('Group 2 catch pitch is greater than Group 1 catch pitch (p = %0.3g)\n',p);
else
    fprintf('Group 1 catch pitch is not significantly different from Group 2 catch pitch (p = %0.3g)\n',p);
end

[h,p]=ttest2(pitch_stimmed_syl1,pitch_stimmed_syl2);
if h==1 && mean(pitch_stimmed_syl1)>mean(pitch_stimmed_syl2)
    fprintf('Group 1 stimmed pitch is greater than Group 2 stimmed pitch (p = %0.3g)\n',p);
elseif h==1 && mean(pitch_stimmed_syl1)<mean(pitch_stimmed_syl2)
    fprintf('Group 2 stimmed pitch is greater than Group 1 stimmed pitch (p = %0.3g)\n',p);
else
    fprintf('Group 1 stimmed pitch is not significantly different from Group 2 stimmed pitch (p = %0.3g)\n',p);
end

delay_from_stim_to_onset_of_stim_syl=[delay_from_stim_to_onset_of_stim_syl1 delay_from_stim_to_onset_of_stim_syl2 delay_from_stim_to_onset_of_stim_syl3];
%plot mean pitch contours
figure
mean_pc_stim1=mean(pitch_contours_stim1,2);
% mean_pc_catch1=mean(pitch_contours_catch1,2);
stderr_pc_stim1=std(pitch_contours_stim1,0,2)/sqrt(length(pitch_stimmed_syl1));
mean_pc_stim2=mean(pitch_contours_stim2,2);
% mean_pc_catch2=mean(pitch_contours_catch2,2);
stderr_pc_stim2=std(pitch_contours_stim2,0,2)/sqrt(length(pitch_stimmed_syl2));
mean_pc_stim3=mean(pitch_contours_stim3,2);
% mean_pc_catch3=mean(pitch_contours_catch3,2);
stderr_pc_stim3=std(pitch_contours_stim3,0,2)/sqrt(length(pitch_stimmed_syl3));
all_pc_catch=[pitch_contours_catch1 pitch_contours_catch2 pitch_contours_catch3];
mean_pc_catch=mean(all_pc_catch,2);
mean_pc_catch1=mean(pitch_contours_catch1,2);
mean_pc_catch2=mean(pitch_contours_catch2,2);
mean_pc_catch3=mean(pitch_contours_catch3,2);
stderr_pc_catch=std(all_pc_catch,0,2)/sqrt(size(all_pc_catch,2));
% time_pc=-4*(length(mean_pc_catch)-1)/(3*fs):4/fs:8*(length(mean_pc_catch)-1)/(3*fs);
% time_pc=linspace(-(length(pstim_clips{1})-1020)^2/(2*fs*4*length(mean_pc_catch1)),(length(pstim_clips{1})-1020)^2/(2*fs*4*length(mean_pc_catch1)),length(mean_pc_catch1));
% time_pc=linspace(-0.08,0.4,length(mean_pc_stim1));
time_pc=time_pc-0.08;

% plot(time_pc,mean_pc_stim1-mean_pc_catch1,'b-',time_pc,mean_pc_stim1+stderr_pc_stim1-mean_pc_catch1,'b--',time_pc,mean_pc_stim1-stderr_pc_stim1-mean_pc_catch1,'b--')
plot(time_pc,mean_pc_catch,'b-','LineWidth',2);
hold on
plot(time_pc,mean_pc_stim1,'k-','LineWidth',2);
plot(time_pc,mean_pc_stim2,'r-','LineWidth',2);
plot(time_pc,mean_pc_stim3,'g-','LineWidth',2);
yl=ylim;
stim_time=-1*mean([delay_from_stim_to_onset_of_stim_syl1 delay_from_stim_to_onset_of_stim_syl2])/1000;
% hqp_time=pq_offset/1000;
plot([stim_time stim_time],yl,'g-','LineWidth',2)
% plot(time_pc,mean_pc_catch1,'y-','LineWidth',2);
% plot(time_pc,mean_pc_catch2,'m-','LineWidth',2);

eval(legend_str)
hold on
plot(time_pc,mean_pc_catch+stderr_pc_catch,'b--','LineWidth',2)
plot(time_pc,mean_pc_catch-stderr_pc_catch,'b--','LineWidth',2)
plot(time_pc,mean_pc_stim1+stderr_pc_stim1,'k--','LineWidth',2)
plot(time_pc,mean_pc_stim1-stderr_pc_stim1,'k--','LineWidth',2)
plot(time_pc,mean_pc_stim2+stderr_pc_stim2,'r--','LineWidth',2)
plot(time_pc,mean_pc_stim2-stderr_pc_stim2,'r--','LineWidth',2)
plot(time_pc,mean_pc_stim3+stderr_pc_stim3,'g--','LineWidth',2)
plot(time_pc,mean_pc_stim3-stderr_pc_stim3,'g--','LineWidth',2)
% plot(time_pc,pitch_contours_catch,'b-')
% plot(time_pc,mean_pc_catch,'b-')
hold on
% plot(time_pc,mean_pc_stim,'r-')
% plot(time_pc,mean_pc_stim2-mean_pc_catch2,'r-',time_pc,mean_pc_stim2+stderr_pc_stim2-mean_pc_catch2,'r--',time_pc,mean_pc_stim2-stderr_pc_stim2-mean_pc_catch2,'r--')

title('Pitch Contours During Stimulation pattern 1 vs 2 vs 1pulse')
xlabel('Time Aligned to Syllable Onset (sec)')
ylabel('Pitch (Hz)')
% plot(-1/1000*mean(delay_from_stim_to_onset_of_stim_syl),yl(1)+0.5*(yl(2)-yl(1)),'go')
% plot([-1/1000*(mean(delay_from_stim_to_onset_of_stim_syl)-std(delay_from_stim_to_onset_of_stim_syl)) -1/1000*(mean(delay_from_stim_to_onset_of_stim_syl)+std(delay_from_stim_to_onset_of_stim_syl))], [yl(1)+0.5*(yl(2)-yl(1)) yl(1)+0.5*(yl(2)-yl(1))],'g-')
yl=get(gca,'ylim');
plot(-1/1000*mean(delay_from_stim_to_onset_of_stim_syl),yl(1)+0.5*(yl(2)-yl(1)),'go')
plot([-1/1000*(mean(delay_from_stim_to_onset_of_stim_syl)-std(delay_from_stim_to_onset_of_stim_syl)) -1/1000*(mean(delay_from_stim_to_onset_of_stim_syl)+std(delay_from_stim_to_onset_of_stim_syl))], [yl(1)+0.5*(yl(2)-yl(1)) yl(1)+0.5*(yl(2)-yl(1))],'g-')

pre_time=fix(0.01/(time_pc(2)-time_pc(1))); %holds 10 ms before (in samples)
post_time=fix(0.1/(time_pc(2)-time_pc(1))); %holds 100ms after (in samples)
time_stim_effects=-0.01:(time_pc(2)-time_pc(1)):0.1;
figure
stim_effects1=nan(size(pitch_contours_stim1));
for k=1:size(pitch_contours_stim1,2)
    [y,i]=min(abs(time_pc+delay_from_stim_to_onset_of_stim_syl1(k)/1000));
    temp=pitch_contours_stim1(i-pre_time:i+post_time,k)-mean_pc_catch1(i-pre_time:i+post_time);
    stim_effects1(1:length(temp),k)=temp;
end
mean_stim_effects1=nanmean(stim_effects1,2);
mean_stim_effects1=mean_stim_effects1(1:length(time_stim_effects));
stderr_stim_effects1=nanstd(stim_effects1,0,2)/sqrt(length(pitch_stimmed_syl1));
stderr_stim_effects1=stderr_stim_effects1(1:length(time_stim_effects));

stim_effects2=nan(size(pitch_contours_stim2));
for k=1:size(pitch_contours_stim2,2)
    [y,i]=min(abs(time_pc+delay_from_stim_to_onset_of_stim_syl2(k)/1000));
    temp=pitch_contours_stim2(i-pre_time:i+post_time,k)-mean_pc_catch2(i-pre_time:i+post_time);
    stim_effects2(1:length(temp),k)=temp;
end
mean_stim_effects2=nanmean(stim_effects2,2);
mean_stim_effects2=mean_stim_effects2(1:length(time_stim_effects));
stderr_stim_effects2=nanstd(stim_effects2,0,2)/sqrt(length(pitch_stimmed_syl2));
stderr_stim_effects2=stderr_stim_effects2(1:length(time_stim_effects));

stim_effects3=nan(size(pitch_contours_stim3));
for k=1:size(pitch_contours_stim3,2)
    [y,i]=min(abs(time_pc+delay_from_stim_to_onset_of_stim_syl3(k)/1000));
    temp=pitch_contours_stim3(i-pre_time:i+post_time,k)-mean_pc_catch3(i-pre_time:i+post_time);
    stim_effects3(1:length(temp),k)=temp;
end
mean_stim_effects3=nanmean(stim_effects3,2);
mean_stim_effects3=mean_stim_effects3(1:length(time_stim_effects));
stderr_stim_effects3=nanstd(stim_effects3,0,2)/sqrt(length(pitch_stimmed_syl3));
stderr_stim_effects3=stderr_stim_effects3(1:length(time_stim_effects));

stim_effects2_save=stim_effects2(1:length(time_stim_effects),:);
stim_effects1_save=stim_effects1(1:length(time_stim_effects),:);
stim_effects3_save=stim_effects3(1:length(time_stim_effects),:);
save(savename,'stim_effects1_save','stim_effects2_save','stim_effects3_save','time_stim_effects');


plot(time_stim_effects,zeros(size(time_stim_effects)),'b-','LineWidth',2)
hold on
plot(time_stim_effects,mean_stim_effects1,'k-','LineWidth',2)
plot(time_stim_effects,mean_stim_effects2,'r-','LineWidth',2)
plot(time_stim_effects,mean_stim_effects3,'g-','LineWidth',2)
hold on
eval(legend_str)
plot(time_stim_effects,mean_stim_effects1+stderr_stim_effects1,'k--','LineWidth',2)
plot(time_stim_effects,mean_stim_effects1-stderr_stim_effects1,'k--','LineWidth',2)
plot(time_stim_effects,mean_stim_effects2+stderr_stim_effects2,'r--','LineWidth',2)
plot(time_stim_effects,mean_stim_effects2-stderr_stim_effects2,'r--','LineWidth',2)
plot(time_stim_effects,mean_stim_effects3+stderr_stim_effects3,'g--','LineWidth',2)
plot(time_stim_effects,mean_stim_effects3-stderr_stim_effects3,'g--','LineWidth',2)
title('Catch Subtracted Pitch effects During Stimulation pattern 1 vs 2 vs 1pulse')
xlabel('Time Aligned to Stim Onset (sec)')
ylabel('Relative Pitch (Hz)')

%plot mean SE contours
figure
mean_sec_stim1=mean(se_contours_stim1,2);
stderr_sec_stim1=std(se_contours_stim1,0,2)/sqrt(length(pitch_stimmed_syl1));
mean_sec_stim2=mean(se_contours_stim2,2);
stderr_sec_stim2=std(se_contours_stim2,0,2)/sqrt(length(pitch_stimmed_syl2));
mean_sec_stim3=mean(se_contours_stim3,2);
stderr_sec_stim3=std(se_contours_stim3,0,2)/sqrt(length(pitch_stimmed_syl3));
all_sec_catch=[se_contours_catch1 se_contours_catch2 se_contours_catch3];
mean_sec_catch1=mean(se_contours_catch1,2);
mean_sec_catch2=mean(se_contours_catch2,2);
mean_sec_catch3=mean(se_contours_catch3,2);
mean_sec_catch=mean(all_sec_catch,2);
stderr_sec_catch=std(all_sec_catch,0,2)/sqrt(size(all_sec_catch,2));
% time_pc=-4*(length(mean_pc_catch)-1)/(3*fs):4/fs:8*(length(mean_pc_catch)-1)/(3*fs);
% time_pc=linspace(-(length(pstim_clips{1})-1020)^2/(2*fs*4*length(mean_pc_catch)),(length(pstim_clips{1})-1020)^2/(2*fs*4*length(mean_pc_catch)),length(mean_pc_catch));
% time_pc=time_pc+0.02;
% time_pc=linspace(-0.08,0.4,length(mean_pc_catch));

% plot(time_pc,mean_pc_stim1-mean_pc_catch1,'b-',time_pc,mean_pc_stim1+stderr_pc_stim1-mean_pc_catch1,'b--',time_pc,mean_pc_stim1-stderr_pc_stim1-mean_pc_catch1,'b--')
plot(time_pc,mean_sec_catch,'b-',time_pc,mean_sec_stim1,'k-',time_pc,mean_sec_stim2,'r-',time_pc,mean_sec_stim3,'g-')
hold on
yl=ylim;
stim_time=-1*mean([delay_from_stim_to_onset_of_stim_syl1 delay_from_stim_to_onset_of_stim_syl2 delay_from_stim_to_onset_of_stim_syl3])/1000;
% hqp_time=pq_offset/1000;
plot([stim_time stim_time],yl,'g-')

eval(legend_str)
hold on
plot(time_pc,mean_sec_catch+stderr_sec_catch,'b--',time_pc,mean_sec_catch-stderr_sec_catch,'b--',time_pc,mean_sec_stim1+stderr_sec_stim1,'k--',...
    time_pc,mean_sec_stim1-stderr_sec_stim1,'k--',time_pc,mean_sec_stim2+stderr_sec_stim2,'r--',time_pc,mean_sec_stim2-stderr_sec_stim2,'r--',...
    time_pc,mean_sec_stim3+stderr_sec_stim3,'g--',time_pc,mean_sec_stim3-stderr_sec_stim3,'g--')
% plot(time_pc,pitch_contours_catch,'b-')
% plot(time_pc,mean_pc_catch,'b-')
hold on
% plot(time_pc,mean_pc_stim,'r-')
% plot(time_pc,mean_pc_stim2-mean_pc_catch2,'r-',time_pc,mean_pc_stim2+stderr_pc_stim2-mean_pc_catch2,'r--',time_pc,mean_pc_stim2-stderr_pc_stim2-mean_pc_catch2,'r--')

title('SE Contours During Stimulation pattern 1 vs 2 vs 1pulse')
xlabel('Time Aligned to Syllable Onset (sec)')
ylabel('SE (arb units)')
yl=get(gca,'ylim');
plot(-1/1000*mean(delay_from_stim_to_onset_of_stim_syl),yl(1)+0.5*(yl(2)-yl(1)),'go')
plot([-1/1000*(mean(delay_from_stim_to_onset_of_stim_syl)-std(delay_from_stim_to_onset_of_stim_syl)) -1/1000*(mean(delay_from_stim_to_onset_of_stim_syl)+std(delay_from_stim_to_onset_of_stim_syl))], [yl(1)+0.5*(yl(2)-yl(1)) yl(1)+0.5*(yl(2)-yl(1))],'g-')

figure

stim_effects_SE1=nan(size(se_contours_stim1));
for k=1:size(se_contours_stim1,2)
    [y,i]=min(abs(time_pc+delay_from_stim_to_onset_of_stim_syl1(k)/1000));
    temp=se_contours_stim1(i-pre_time:i+post_time,k)-mean_sec_catch1(i-pre_time:i+post_time);
    stim_effects_SE1(1:length(temp),k)=temp;
end
mean_stim_effects_SE1=nanmean(stim_effects_SE1,2);
mean_stim_effects_SE1=mean_stim_effects_SE1(1:length(time_stim_effects));
stderr_stim_effects_SE1=nanstd(stim_effects_SE1,0,2)/sqrt(length(pitch_stimmed_syl1));
stderr_stim_effects_SE1=stderr_stim_effects_SE1(1:length(time_stim_effects));

stim_effects_SE2=nan(size(se_contours_stim2));
for k=1:size(se_contours_stim2,2)
    [y,i]=min(abs(time_pc+delay_from_stim_to_onset_of_stim_syl2(k)/1000));
    temp=se_contours_stim2(i-pre_time:i+post_time,k)-mean_sec_catch2(i-pre_time:i+post_time);
    stim_effects_SE2(1:length(temp),k)=temp;
end
mean_stim_effects_SE2=nanmean(stim_effects_SE2,2);
mean_stim_effects_SE2=mean_stim_effects_SE2(1:length(time_stim_effects));
stderr_stim_effects_SE2=nanstd(stim_effects_SE2,0,2)/sqrt(length(pitch_stimmed_syl2));
stderr_stim_effects_SE2=stderr_stim_effects_SE2(1:length(time_stim_effects));

stim_effects_SE3=nan(size(se_contours_stim3));
for k=1:size(se_contours_stim3,2)
    [y,i]=min(abs(time_pc+delay_from_stim_to_onset_of_stim_syl3(k)/1000));
    temp=se_contours_stim3(i-pre_time:i+post_time,k)-mean_sec_catch3(i-pre_time:i+post_time);
    stim_effects_SE3(1:length(temp),k)=temp;
end
mean_stim_effects_SE3=nanmean(stim_effects_SE3,2);
mean_stim_effects_SE3=mean_stim_effects_SE3(1:length(time_stim_effects));
stderr_stim_effects_SE3=nanstd(stim_effects_SE3,0,2)/sqrt(length(pitch_stimmed_syl3));
stderr_stim_effects_SE3=stderr_stim_effects_SE3(1:length(time_stim_effects));


plot(time_stim_effects,zeros(size(time_stim_effects)),'b-','LineWidth',2)
hold on
plot(time_stim_effects,mean_stim_effects_SE1,'k-','LineWidth',2)
plot(time_stim_effects,mean_stim_effects_SE2,'r-','LineWidth',2)
plot(time_stim_effects,mean_stim_effects_SE3,'g-','LineWidth',2)
hold on
eval(legend_str)
plot(time_stim_effects,mean_stim_effects_SE1+stderr_stim_effects_SE1,'k--','LineWidth',2)
plot(time_stim_effects,mean_stim_effects_SE1-stderr_stim_effects_SE1,'k--','LineWidth',2)
plot(time_stim_effects,mean_stim_effects_SE2+stderr_stim_effects_SE2,'r--','LineWidth',2)
plot(time_stim_effects,mean_stim_effects_SE2-stderr_stim_effects_SE2,'r--','LineWidth',2)
plot(time_stim_effects,mean_stim_effects_SE3+stderr_stim_effects_SE3,'g--','LineWidth',2)
plot(time_stim_effects,mean_stim_effects_SE3-stderr_stim_effects_SE3,'g--','LineWidth',2)
title('Catch Subtracted SE effects During Stimulation pattern 1 vs 2 vs 1pulse')
xlabel('Time Aligned to Stim Onset (sec)')
ylabel('Relative SE')



figure
% fs=32;
stim_clips_smoothed1=20*log10(stim_clips_smoothed1./pref)+94;
mean_stim_clips_smoothed1=mean(stim_clips_smoothed1,1);
stim_clips_smoothed2=20*log10(stim_clips_smoothed2./pref)+94;
mean_stim_clips_smoothed2=mean(stim_clips_smoothed2,1);
stim_clips_smoothed3=20*log10(stim_clips_smoothed3./pref)+94;
mean_stim_clips_smoothed3=mean(stim_clips_smoothed3,1);
catch_clips_smoothed1=20*log10(catch_clips_smoothed1./pref)+94;
mean_catch_clips_smoothed1=mean(catch_clips_smoothed1,1);
catch_clips_smoothed2=20*log10(catch_clips_smoothed2./pref)+94;
mean_catch_clips_smoothed2=mean(catch_clips_smoothed2,1);
catch_clips_smoothed3=20*log10(catch_clips_smoothed3./pref)+94;
mean_catch_clips_smoothed3=mean(catch_clips_smoothed3,1);
catch_clips_smoothed=[catch_clips_smoothed1;catch_clips_smoothed2;catch_clips_smoothed3];
[r1,c1]=size(stim_clips_smoothed1);
[r2,c2]=size(stim_clips_smoothed2);
[r3,c3]=size(stim_clips_smoothed3);
[r,c]=size(catch_clips_smoothed);
time=-0.168:1/fs:0.168;
plot(time,mean(catch_clips_smoothed,1),'b-',time,mean(stim_clips_smoothed1,1),'k-',time,mean(stim_clips_smoothed2,1),'r-',time,mean(stim_clips_smoothed3,1),'g-')
hold on
yl=ylim;
plot([stim_time stim_time],yl,'g-')
eval(legend_str)
plot(time,mean(stim_clips_smoothed1,1)+std(stim_clips_smoothed1,0,1)/sqrt(r1),'k--',time,mean(stim_clips_smoothed1,1)-std(stim_clips_smoothed1,0,1)/sqrt(r1),'k--')
plot(time,mean(stim_clips_smoothed2,1)+std(stim_clips_smoothed2,0,1)/sqrt(r2),'r--',time,mean(stim_clips_smoothed2,1)-std(stim_clips_smoothed2,0,1)/sqrt(r2),'r--')
plot(time,mean(stim_clips_smoothed3,1)+std(stim_clips_smoothed3,0,1)/sqrt(r3),'g--',time,mean(stim_clips_smoothed3,1)-std(stim_clips_smoothed3,0,1)/sqrt(r3),'g--')
plot(time,mean(catch_clips_smoothed,1)+std(catch_clips_smoothed,0,1)/sqrt(r),'b--',time,mean(catch_clips_smoothed,1)-std(catch_clips_smoothed,0,1)/sqrt(r),'b--')
% plot(time,max(catch_clips_smoothed,1),'b--',time,min(catch_clips_smoothed,1),'b--')
% plot(time,max(stim_clips_smoothed,1),'b--',time,min(stim_clips_smoothed,1),'b--')
xlabel('Time relative to syllable onset (ms)')
ylabel('Smooth Rectified Amplitude (dB)')
% axis([-160 75 min(mean(stim_clips_smoothed(:,fix(end/2)),1))-.2 max(mean(stim_clips_smoothed,1))+.2])
yl=get(gca,'ylim');
plot(-1/1000*mean(delay_from_stim_to_onset_of_stim_syl),yl(1)+0.5*(yl(2)-yl(1)),'go')
plot([-1/1000*(mean(delay_from_stim_to_onset_of_stim_syl)-std(delay_from_stim_to_onset_of_stim_syl)) -1/1000*(mean(delay_from_stim_to_onset_of_stim_syl)+std(delay_from_stim_to_onset_of_stim_syl))], [yl(1)+0.5*(yl(2)-yl(1)) yl(1)+0.5*(yl(2)-yl(1))],'g-')


figure
pre_time_amp=0.01*fs;
post_time_amp=0.1*fs;
time_stim_effects_amp=-0.01:(time(2)-time(1)):0.1;
stim_effects_amp1=nan(size(stim_clips_smoothed1));
for k=1:size(stim_clips_smoothed1,1)
    [y,i]=min(abs(time+delay_from_stim_to_onset_of_stim_syl1(k)/1000));
    temp=stim_clips_smoothed1(k,i-pre_time_amp:min(i+post_time_amp,size(stim_clips_smoothed1,2)))-mean_catch_clips_smoothed1(i-pre_time_amp:min(i+post_time_amp,size(stim_clips_smoothed1,2)));
    stim_effects_amp1(k,1:length(temp))=temp;
end
mean_stim_effects_amp1=nanmean(stim_effects_amp1,1);
mean_stim_effects_amp1=mean_stim_effects_amp1(1:length(time_stim_effects_amp));
stderr_stim_effects_amp1=nanstd(stim_effects_amp1,1)/sqrt(length(pitch_stimmed_syl1));
stderr_stim_effects_amp1=stderr_stim_effects_amp1(1:length(time_stim_effects_amp));

stim_effects_amp2=nan(size(stim_clips_smoothed2));
for k=1:size(stim_clips_smoothed2,1)
    [y,i]=min(abs(time+delay_from_stim_to_onset_of_stim_syl2(k)/1000));
    temp=stim_clips_smoothed2(k,i-pre_time_amp:min(i+post_time_amp,size(stim_clips_smoothed2,2)))-mean_catch_clips_smoothed2(i-pre_time_amp:min(i+post_time_amp,size(stim_clips_smoothed2,2)));
    stim_effects_amp2(k,1:length(temp))=temp;
end
mean_stim_effects_amp2=nanmean(stim_effects_amp2,1);
mean_stim_effects_amp2=mean_stim_effects_amp2(1:length(time_stim_effects_amp));
stderr_stim_effects_amp2=nanstd(stim_effects_amp2,1)/sqrt(length(pitch_stimmed_syl2));
stderr_stim_effects_amp2=stderr_stim_effects_amp2(1:length(time_stim_effects_amp));

stim_effects_amp3=nan(size(stim_clips_smoothed3));
for k=1:size(stim_clips_smoothed3,1)
    [y,i]=min(abs(time+delay_from_stim_to_onset_of_stim_syl3(k)/1000));
    temp=stim_clips_smoothed3(k,i-pre_time_amp:min(i+post_time_amp,size(stim_clips_smoothed3,2)))-mean_catch_clips_smoothed3(i-pre_time_amp:min(i+post_time_amp,size(stim_clips_smoothed3,2)));
    stim_effects_amp3(k,1:length(temp))=temp;
end
mean_stim_effects_amp3=nanmean(stim_effects_amp3,1);
mean_stim_effects_amp3=mean_stim_effects_amp3(1:length(time_stim_effects_amp));
stderr_stim_effects_amp3=nanstd(stim_effects_amp3,1)/sqrt(length(pitch_stimmed_syl3));
stderr_stim_effects_amp3=stderr_stim_effects_amp3(1:length(time_stim_effects_amp));


plot(time_stim_effects_amp,zeros(size(time_stim_effects_amp)),'b-','LineWidth',2)
hold on
plot(time_stim_effects_amp,mean_stim_effects_amp1,'k-','LineWidth',2)
plot(time_stim_effects_amp,mean_stim_effects_amp2,'r-','LineWidth',2)
plot(time_stim_effects_amp,mean_stim_effects_amp3,'g-','LineWidth',2)
hold on
eval(legend_str)
plot(time_stim_effects_amp,mean_stim_effects_amp1+stderr_stim_effects_amp1,'k--','LineWidth',2)
plot(time_stim_effects_amp,mean_stim_effects_amp1-stderr_stim_effects_amp1,'k--','LineWidth',2)
plot(time_stim_effects_amp,mean_stim_effects_amp2+stderr_stim_effects_amp2,'r--','LineWidth',2)
plot(time_stim_effects_amp,mean_stim_effects_amp2-stderr_stim_effects_amp2,'r--','LineWidth',2)
plot(time_stim_effects_amp,mean_stim_effects_amp3+stderr_stim_effects_amp3,'g--','LineWidth',2)
plot(time_stim_effects_amp,mean_stim_effects_amp3-stderr_stim_effects_amp3,'g--','LineWidth',2)
title('Catch Subtracted amp effects During Stimulation pattern 1 vs 2 vs 1pulse')
xlabel('Time Aligned to Stim Onset (sec)')
ylabel('Relative amp')


disp(mean(delay_from_stim_to_onset_of_stim_syl));