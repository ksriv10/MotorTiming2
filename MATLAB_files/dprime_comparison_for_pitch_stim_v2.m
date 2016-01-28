function [] = dprime_comparison_for_pitch_stim_v2( pc_file,load_file)
%performs d prime analysis between two stim trains effects on air sac
%pressure; old is a binary indicating whether the file 
if nargin==1
    load_file=0;
end
id=strfind(pc_file,'.mat');
syl_str=pc_file(id-2:id-1);
if ~isempty(strfind(pc_file,'3663'))
    title_str=['6_3_vs_3_6' syl_str];
    legend_str=sprintf('legend(''6 3'',''3 6'',''1 pulse'')');
else
    title_str=['9_3_vs_3_9' syl_str];
%     legend_str=sprintf('legend(''9 3'',''3 9'',''1 pulse'')');
    legend_str=sprintf('legend(''9 3'',''3 9'')');
end
if ~load_file
    numBootstrapTrials=1000;
    

    load(pc_file);

    trials_1pulse=stim_effects3_save';
    trials1=stim_effects1_save';
    trials2=stim_effects2_save';
    figure(1);subplot(2,1,1);plot(trials1')
    subplot(2,1,2);plot(trials2')
    lower_bound=input('What is lower bound index on x-axis for clean pitch measurement?');
    upper_bound=input('What is upper bound index on x-axis for clean pitch measurement?');
    temp_stds1=std(trials1(:,lower_bound:upper_bound),0,2);
    temp_stds2=std(trials2(:,lower_bound:upper_bound),0,2);
    temp_stds3=std(trials_1pulse(:,lower_bound:upper_bound),0,2);
    close all
    hist([temp_stds1;temp_stds2;temp_stds3],100);
    cutoff1=input('What is the separation point for stds in the histogram?');
    close all
%     hist(temp_stds2,100);
%     cutoff2=input('What is the separation point for stds in the histogram?');
    
    trials1=trials1(temp_stds1<cutoff1,:);
    trials2=trials2(temp_stds2<cutoff1,:);
    trials_1pulse=trials_1pulse(temp_stds3<cutoff1,:);
    size(trials1,1)
    size(trials2,1)
    size(trials_1pulse,1)
    
    
    mean_trials1=nanmean(trials1);
    mean_trials2=nanmean(trials2);    
    mean_trials_1pulse=nanmean(trials_1pulse); 
    sterr_trials1=nanstd(trials1)/sqrt(size(trials1,1));
    sterr_trials2=nanstd(trials2)/sqrt(size(trials2,1));
    sterr_trials_1pulse=nanstd(trials_1pulse)/sqrt(size(trials_1pulse,1));
    

    

    timebins=length(time_stim_effects);

    numOne = size(trials1,1);
    numTwo = size(trials2,1);

    replacedPC1 = zeros(numOne,timebins);
    replacedPC2 = zeros(numTwo,timebins);


    comparisonByWord = zeros(numBootstrapTrials,timebins);


    for i = 1:numBootstrapTrials;
        if mod(100*i/numBootstrapTrials,10)==0
            disp([num2str(100*i/numBootstrapTrials) '% done'])
        end
        samplesChosen1 = randsample(1:numOne, numOne, true);
        samplesChosen2 = randsample(1:numTwo, numTwo, true);
        for j = 1: numOne;
            replacedPC1(j, :) = trials1(samplesChosen1(j),:);
        end
        for j = 1: numTwo;
            replacedPC2(j, :) = trials2(samplesChosen2(j),:);
        end

        meanOne = mean(replacedPC1);
        stdOne = std(replacedPC1);
        meanTwo = mean(replacedPC2);
        stdTwo = std(replacedPC2);

        for j = 1:timebins;
        comparisonByWord(i,j) = ((meanOne(j)-meanTwo(j))/(sqrt((stdOne(j)^2) + (stdTwo(j)^2))));
        end
    end

    stdCompByWord = nanstd(comparisonByWord);
    meanCompByWord = nanmean(comparisonByWord);
    id=strfind(pc_file,'_');
    savename=[pc_file(1:id) 'dprime_v2_' title_str];
    save(savename,'stdCompByWord','meanCompByWord','mean_trials1','mean_trials2','mean_trials_1pulse','sterr_trials_1pulse','sterr_trials1','sterr_trials2','time_stim_effects');
else
    load(pc_file)
end
figure
hold on
ax(1)=subplot(2,1,1);
x=time_stim_effects(1:3*end/4);
y1=mean_trials1+sterr_trials1;
y2=mean_trials1-sterr_trials1;
X=[x fliplr(x)];
Y=[y1(1:3*end/4) fliplr(y2(1:3*end/4))];
p=fill(X,Y,'b');
set(p,'EdgeColor','none');
hold on
y1=mean_trials2+sterr_trials2;
y2=mean_trials2-sterr_trials2;
X=[x fliplr(x)];
Y=[y1(1:3*end/4) fliplr(y2(1:3*end/4))];
p=fill(X,Y,'r');
set(p,'EdgeColor','none');
plot(x,mean_trials1(1:3*end/4),'k-',x,mean_trials2(1:3*end/4),'k-')
% hold on
% y1=mean_trials_1pulse+sterr_trials_1pulse;
% y2=mean_trials_1pulse-sterr_trials_1pulse;
% X=[x fliplr(x)];
% Y=[y1(1:3*end/4) fliplr(y2(1:3*end/4))];
% fill(X,Y,'g')
eval(legend_str)
% axis([-0.005 0.035 -150 300])

% plot(time_stim_effects,meanCompByWord,'b-',time_stim_effects,meanCompByWord+stdCompByWord,'b--',time_stim_effects,meanCompByWord-stdCompByWord,'b--')
% hold on
% plot([min(time2) max(time2)],[0 0],'k-')
% xlabel('Time relative to stim onset (s)')
% ylabel('d-prime')
% title(strrep(title_str,'_',' '))

plot([min(time_stim_effects) max(time_stim_effects)],[0 0],'k-')
ylabel('Pitch Change (Hz)')
title(strrep(title_str,'_',' '))
ax(2)=subplot(2,1,2);

y1=-1*meanCompByWord+stdCompByWord;
y2=-1*meanCompByWord-stdCompByWord;
X=[x fliplr(x)];
Y=[y1(1:3*end/4) fliplr(y2(1:3*end/4))];
p=fill(X,Y,'g');
set(p,'EdgeColor','none');
% plot(time2,meanCompByWord_N,'b-',time2,meanCompByWord_N+stdCompByWord_N,'b--',time2,meanCompByWord_N-stdCompByWord_N,'b--')
hold on
plot(x,-1*meanCompByWord(1:3*end/4),'k-')
plot([min(time_stim_effects) max(time_stim_effects)],[0 0],'k-')
ylabel('d-prime using Pitch Change')

xlabel('Time relative to stim onset (s)')


linkaxes(ax,'x')
xlim([-0.005 0.04])

% subplot(2,1,2);
% axis([-0.005 0.035 -1 1])
[y,ind1]=min(abs(time_stim_effects));
[y,ind2]=min(abs(time_stim_effects-0.035));
[y,ind]=max(abs(meanCompByWord(ind1:ind2)));
max_dprime=meanCompByWord(ind1+ind)
std_dprime=stdCompByWord(ind1+ind)

end

