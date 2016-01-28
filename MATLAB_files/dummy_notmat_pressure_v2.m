function dummy_notmat_pressure_v2(rhd)
%make_notmat
% generates make .not.mat files
%
% syntax:
% make_notmat(min_dur,min_int,threshold,sm_win)
%
% input arguments: segmenting parameters
%   min_dur = minimum syllable duration in ms, default is 20
%   min_int = minimum inter-syllable interval in ms, default is 2
%   threshold = that amplitude crosses, default is 1k
%   sm_win = size of smoothing window, default is 2
%
% last edited: David Nicholson 06/20/12. commented and added ability to
% specifiy segmenting parameters

    min_dur=100;
    min_int=20;
    threshold=0;
    sm_win=2;

pressure_channel=1;
    
if rhd    
    if isunix
        !ls **.rhd > batchfile
    else
        !dir /B **.rhd > batchfile
    end
else
    if isunix
        !ls **.cbin > batchfile
    else
        !dir /B **.cbin > batchfile
    end
end
load('calibration_later_file.mat')
cal_matrix=[0 pressure_zero;pressure_cal pressure_one];
slope=(cal_matrix(1,1)-cal_matrix(2,1))/(cal_matrix(1,2)-cal_matrix(2,2));
ct=0;n_files=0;fid=fopen('batchfile');while 1;fname=fgetl(fid);if (~ischar(fname));break;end;n_files=n_files+1;end;fclose(fid);
% mean_pressure=[];
% norm_mean_pressure=[];
fid=fopen('batchfile');
while 1
    fname=fgetl(fid);
    if (~ischar(fname));break;end
    if rhd
        fname_mat=[fname(1:end-3) 'mat'];
        load(fname_mat,'board_adc_data','frequency_parameters')
        dat=board_adc_data(pressure_channel,:);
        fs=frequency_parameters.board_adc_sample_rate;
    else
        fname_mat=[fname(1:end-4) 'mat'];
        load(fname_mat,'dat','fs')
    end
    dat=dat/mean(dat);
%     mean_pressure(end+1)=mean(dat);
    dat=0.249174*slope*(dat*first_file_mean-cal_matrix(1,2));  %coefficient converts from inH2O to kPa
%     dat=dat*first_file_mean;
%     norm_mean_pressure(end+1)=mean(dat);
%     Fs=frequency_parameters.board_adc_sample_rate;
%     [dat,Fs]=evsoundin('',fname,'obs0');
    
    DOFILT=1;
%     sm=SmoothData_pressure(dat,Fs,DOFILT,0,512,0.8,2.0,1,50);
    sm=bandpass_filtfilt(dat,fs,1,50,'hanningfir');
%     sm=SmoothData_EMG(dat,Fs,DOFILT);
    sm(1)=0.0;sm(end)=0.0;
    
    [onsets, offsets]=SegmentNotes_local(sm, fs, min_int, min_dur, threshold);
    labels = char(ones([1,length(onsets)])*fix('Q'));
    
    onsets=onsets*1000;  % put in msec to match .not.mat output of evsonganaly.m
    offsets= offsets*1000;
    
    notmat_fname=[fname_mat '.not.mat'];
    save(notmat_fname,'fs','fname','labels','min_dur','min_int','threshold','sm_win','onsets','offsets','-v7.3')
    ct=ct+1;disp(['Processed ' num2str(ct) ' of ' num2str(n_files) ' files'])
end
% save('temp_mat','mean_pressure','norm_mean_pressure')
fclose(fid);


%  This is copied from SmoothData.m, part of the evsonganaly.m package
% modified by removing spectrogram component
%
function smooth=SmoothData_local(rawsong,Fs,DOFILT,nfft,olap,sm_win,F_low,F_High);
% [smooth,spec,t,f]=evsmooth(rawsong,Fs,DOFILT,nfft,olap,sm_win,F_low,F_High,DOFILT);
% returns the smoothed waveform/envelope + the spectrum
%


if (~exist('F_low','var'));F_low  = 500.0;end
if (~exist('F_high','var'));F_high = 10000.0;end
if (~exist('nfft','var'));nfft = 512;end
if (~exist('olap','var'));olap = 0.5;end
if (~exist('sm_win','var'));sm_win = 2.0;end % msec
if (~exist('DOFILT','var'));DOFILT=1;end

filter_type = 'hanningfir';

if (DOFILT==1)
    %filtsong=bandpass(rawsong,Fs,F_low,F_high,filter_type);
    %disp('SmoothData.m called - now using filtfilt')
    filtsong=bandpass_filtfilt(rawsong,Fs,F_low,F_high,filter_type);
else
    filtsong=rawsong;
end

squared_song = filtsong.^2;

%smooth the rectified song
len = round(Fs*sm_win/1000);
h   = ones(1,len)/len;
smooth = conv(h, squared_song);
offset = round((length(smooth)-length(filtsong))/2);
smooth=smooth(1+offset:length(filtsong)+offset);


%  This is copied from SegmentNotes.m, part of the evsonganaly.m package
%
function [onsets, offsets]=SegmentNotes_local(smooth, Fs, min_int, min_dur, threshold)
% [ons,offs]=evsegment(smooth,Fs,min_int,min_dur,threshold);
% segment takes smoothed filtered song and returns vectors of note
% onsets and offsets values are in seconds

h=[1 -1];

%threshold input
%notetimes=abs(diff(smooth))>threshold;
notetimes=double(smooth>threshold);

%extract index values for note onsets and offsets
trans=conv(h,notetimes);
t_onsets  = find(trans>0);
t_offsets = find(trans<0);

onsets = t_onsets;offsets=t_offsets;
if ((length(onsets)<1)|(length(offsets)<1))
    onsets=[];offsets=[];
    return;
end

if (length(t_onsets) ~= length(t_offsets))
    disp('number of note onsets and offsets do not match')
else
    %eliminate short intervals
    temp_int=(onsets(2:length(onsets))-offsets(1:length(offsets)-1))*1000/Fs;
    real_ints=temp_int>min_int;
    onsets=[onsets(1); nonzeros(onsets(2:length(onsets)).*real_ints)];
    offsets=[nonzeros(offsets(1:length(offsets)-1).*real_ints); offsets(length(offsets))];
    
    %eliminate short notes
    temp_dur=(offsets-onsets)*1000/Fs;
    real_durs=temp_dur>min_dur;
    onsets=[nonzeros((onsets).*real_durs)];
    offsets=[nonzeros((offsets).*real_durs)];
    
    %convert to ms: peculiarities here are to prevent rounding problem
    % if t_ons is simply replaced with onsets, everything gets rounded
    onsets = onsets/Fs; % all in seconds
    offsets = offsets/Fs; %all in seconds
end

