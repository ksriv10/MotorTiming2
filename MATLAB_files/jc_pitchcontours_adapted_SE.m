function [pitch_data, SE_data, T1]=jc_pitchcontours_adapted_SE(fvals,N,OVERLAP,sigma,F_low,F_high,harms,filetype)
%F_low and F_high are the boundary frequencies in which to look for a peak
    %in the autocorrelation function (i.e. in which to calculate the pitch at
    %each point in time).

% sigma of 1.0 is generally good for BF song (with harmonic spacing >>1500Hz)
% sigma of 1.5-2.0 is generally good for ZF song (with harmonic spacing <<1500Hz)
    % this is based on the consensus function from Tim Gardner

% fvals is the output of one of the findwnote functions.  If you don't use
% these, note that each column of shifted (see line ) is about 100ms of
% song aligned to the beginning of the syllable of interest.

% harms allows you to weight the  pitch estimate at each time point based 
% on the power at each harmonic.  Usually, I just choose the 
% single harmonic that tends to have the most power. Either way gives a similar solution.
% Looking only at the first harmonic
    % pitchcontours_pre=jc_pitchcontours(fvals_pre,1024,1020,1,2000,3000,[1],'obs0');
% Looking at all three harmonics below 10000
    % pitchcontours_pre=jc_pitchcontours(fvals_pre,1024,1020,1,2000,3000,[1 2 3],'obs0');
    
    %kyle's example: pitch_contours_stim=jc_pitchcontours_adapted2(pstim_clips,1024,1020,1,2500,6500,[1],'obs0');
    %where pstim_clips is an array with each cell being the sound data for
    %a different "trial". Pitch_contours_stim is a matrix of the different
    %pitch contours. remember that sampling rate of output is different
    %than that of input. I used this for other data: time_pc=linspace(-(length(pstim_clips{1})-1020)^2/(2*fs*4*length(mean_pc_catch)),(length(pstim_clips{1})-1020)^2/(2*fs*4*length(mean_pc_catch)),length(mean_pc_catch));

   
shifted=zeros(length(fvals),length(fvals{1}));    
for i=1:length(fvals)
    shifted(i,:)=fvals{i};
end
m=size(shifted);
NumberOfNotes=m(1);

%Parameters that one may want to change.
if strcmp(filetype,'obs0')
    SAMPLING=32000; 
else SAMPLING=44100;
end

%sonodvA program
t=-N/2+1:N/2;
sigma=(sigma/1000)*SAMPLING;
%Gaussian and first derivative as windows.
w=exp(-(t/sigma).^2);
timebins=floor((size(shifted,2)-(N))/(N-OVERLAP))+1;
freqbins=(N/2)+1;

Nyquist=SAMPLING/2;
step=Nyquist/freqbins; %This is the width (in ms) between bins.
mini=round(F_high/step);
maxi=round(F_low/step);
highest_harmonic=1; 

sonogram=zeros(freqbins,timebins);
pitch_data=zeros(timebins,NumberOfNotes);
SE_data=zeros(timebins,NumberOfNotes);
for found_note=1:NumberOfNotes
    if mod(fix(100*found_note/NumberOfNotes),5)==0
        disp([num2str(ceil(100*found_note/NumberOfNotes)) '% complete'])
    end
    if found_note==1
        [S1,F1,T1,P1]=spectrogram(shifted(found_note,:),w,OVERLAP,N,SAMPLING);
        sonogram=abs(flipdim(S1,1)); %gaussian windowed spectrogram    
        
    else
        sonogram=abs(flipdim(spectrogram(shifted(found_note,:),w,OVERLAP,N,SAMPLING),1)); %gaussian windowed spectrogram
    end

    %This loop determines the pitch (fundamental frequency) at each time bin by
    %taking the peak of the autocorrelation function within a range specified
    %by the user
    for currenttime_bin=1:size(sonogram,2)
        slice=sonogram(:,currenttime_bin); %power at each frequency for the current time bin
        Powerful(1)=0;
        Freqbinest(1)=0;
        for i=1:length(harms)
            minimum=freqbins-mini*harms(i);  
            maximum=freqbins-maxi*harms(i);
            freq_window=slice(minimum:maximum);
            [Pow,Ind]=max(freq_window);  
            % Interpolation--subtraction of 0.5 because we care about the
            % central frequency of the bin, not the upper boundary.
            if Ind==1 || Ind==length(freq_window)
%                 if Ind==1
%                     disp([extra_str 'Peak at LOWER extreme of range - looking for convex peak between limits'])
%                 elseif Ind==length(freq_window)
%                     disp([extra_str 'Peak at UPPER extreme of range - looking for convex peak between limits'])
%                 end
                id_neg_2nd_der=find(sign(diff(sign(diff(freq_window))))<0)+1;   % this is length(spect_slice)-2
                max_p_id=id_neg_2nd_der(find(freq_window(id_neg_2nd_der)==max(freq_window(id_neg_2nd_der))));
                if isempty(max_p_id)
                    Indest=Ind;
                else
                    Indest=pinterp([max_p_id-1;max_p_id;max_p_id+1], [freq_window(max_p_id-1);freq_window(max_p_id);freq_window(max_p_id+1)]);
                end
            else
                Indest=pinterp([Ind-1;Ind;Ind+1], [freq_window(Ind-1);freq_window(Ind);freq_window(Ind+1)]);
            end
            Real_Index=minimum+Indest-1; % -1 to account for window;
            Freqbinest(i)=(freqbins-Real_Index)/harms(i);
            Powerful(i)=Pow;
        end
        abs_S1 = abs(freq_window);
        abs_S1_rat = abs_S1/sum(abs_S1);
        normalizer=sum(Powerful);
        normpower=Powerful/normalizer;
        freqbin_estimate(currenttime_bin)=dot(normpower,Freqbinest);
        se_estimate(currenttime_bin)=-sum(abs_S1_rat.*log10(abs_S1_rat));
    end
    pitch_data(:,found_note)=freqbin_estimate*(Nyquist/(freqbins-1));
    SE_data(:,found_note)=se_estimate;
end