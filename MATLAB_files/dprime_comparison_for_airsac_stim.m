function [] = dprime_comparison_for_airsac_stim( pressure_file,old,load_file)
%performs d prime analysis between two stim trains effects on air sac
%pressure; old is a binary indicating whether the file 
if nargin==2
    load_file=0;
end
title_str='10_10_vs_12_8_phase_aligned';
% title_str='2_18_vs_4_16_phase_aligned';
if ~load_file
    numBootstrapTrials=1000;
    comparison_patterns=[1 2]; %not a great way to automate this, so select which cell in the array you want to use for comparison

    load(pressure_file,'stim_aligned_single_catch_sub_sorted','time2','stim_latencies_sorted');
    if 0
        for j=1:2
            id_delay=find(stim_latencies_sorted{comparison_patterns(j)}>-0.05 & stim_latencies_sorted{comparison_patterns(j)}<0);
            stim_latencies_sorted{comparison_patterns(j)}=stim_latencies_sorted{comparison_patterns(j)}(id_delay);
            stim_aligned_single_catch_sub_sorted{comparison_patterns(j)}=stim_aligned_single_catch_sub_sorted{comparison_patterns(j)}(id_delay,:);
        end
    end

    if old
        trials1=stim_aligned_single_catch_sub_sorted{comparison_patterns(1)};
        trials2=stim_aligned_single_catch_sub_sorted{comparison_patterns(2)};
    else
        trials1=stim_aligned_single_catch_sub_sorted{comparison_patterns(1)}-repmat(nanmean(stim_aligned_single_catch_sub_sorted{end-1}),size(stim_aligned_single_catch_sub_sorted{comparison_patterns(1)},1),1);
        trials2=stim_aligned_single_catch_sub_sorted{comparison_patterns(2)}-repmat(nanmean(stim_aligned_single_catch_sub_sorted{end-1}),size(stim_aligned_single_catch_sub_sorted{comparison_patterns(2)},1),1);
    end
    
    mean_trials1=nanmean(trials1);
    mean_trials2=nanmean(trials2);    
    sterr_trials1=nanstd(trials1)/sqrt(size(trials1,1));
    sterr_trials2=nanstd(trials2)/sqrt(size(trials2,1));

    timebins=length(time2);

    numOne = size(trials1,1);
    numTwo = size(trials2,1);

    replacedPressure1 = zeros(numOne,timebins);
    replacedPressure2 = zeros(numTwo,timebins);


    comparisonByWord = zeros(numBootstrapTrials,timebins);


    for i = 1:numBootstrapTrials;
        if mod(100*i/numBootstrapTrials,10)==0
            disp([num2str(100*i/numBootstrapTrials) '% done'])
        end
        samplesChosen1 = randsample(1:numOne, numOne, true);
        samplesChosen2 = randsample(1:numTwo, numTwo, true);
        for j = 1: numOne;
            replacedPressure1(j, :) = trials1(samplesChosen1(j),:);
        end
        for j = 1: numTwo;
            replacedPressure2(j, :) = trials2(samplesChosen2(j),:);
        end

        meanOne = mean(replacedPressure1);
        stdOne = std(replacedPressure1);
        meanTwo = mean(replacedPressure2);
        stdTwo = std(replacedPressure2);

        for j = 1:timebins;
        comparisonByWord(i,j) = ((meanOne(j)-meanTwo(j))/(sqrt((stdOne(j)^2) + (stdTwo(j)^2))));
        end
    end

    stdCompByWord = std(comparisonByWord);
    meanCompByWord = mean(comparisonByWord);
    id=strfind(pressure_file,'_');
    savename=[pressure_file(1:id) 'dprime_' title_str];
    save(savename,'comparisonByWord','stdCompByWord','meanCompByWord','time2','mean_trials1','mean_trials2','sterr_trials1','sterr_trials2');
else
    load(pressure_file)
end

% load('Pressure_dprime_2_18_vs_4_16.mat')
figure(1)
% ax(1)=subplot(2,1,1);
% time2=time2(time2<0.06);
x=time2;
% y1=meanCompByWord(time2<0.06)+stdCompByWord(time2<0.06);
% y2=meanCompByWord(time2<0.06)-stdCompByWord(time2<0.06);
y1=meanCompByWord+stdCompByWord;
y2=meanCompByWord-stdCompByWord;
X=[x fliplr(x)];
Y=[y1 fliplr(y2)];

color='b';
hold on
fill(X,Y,color)
% plot(time2,meanCompByWord,'b-',time2,meanCompByWord+stdCompByWord,'b--',time2,meanCompByWord-stdCompByWord,'b--')
% hold on
plot([min(time2) max(time2)],[0 0],'k-')
xlabel('Time relative to stim onset (s)')
ylabel('d-prime')
title(strrep(title_str,'_',' '))
% title('2 18 vs 4 16')
% 
% load('Pressure_dprime_10_10_vs_12_8.mat')
% figure(1)
% ax(2)=subplot(2,1,2);
% 
% y1=meanCompByWord(time2<0.06)+stdCompByWord(time2<0.06);
% y2=meanCompByWord(time2<0.06)-stdCompByWord(time2<0.06);
% X=[x fliplr(x)];
% Y=[y1 fliplr(y2)];
% 
% hold on
% fill(X,Y,color)
% % plot(time2,meanCompByWord,'b-',time2,meanCompByWord+stdCompByWord,'b--',time2,meanCompByWord-stdCompByWord,'b--')
% % hold on
% plot([min(time2) max(time2)],[0 0],'k-')
% xlabel('Time relative to stim onset (s)')
% ylabel('d-prime')
% % title(strrep(title_str,'_',' '))
% title('10 10 vs 12 8')
% 
% linkaxes(ax,'x')




end

