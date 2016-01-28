function [] = dprime_comparison_for_IV_force_v2(force_file,load_file)
%performs d prime analysis between two stim trains effects on air sac
%pressure; old is a binary indicating whether the file 
if nargin==1
    load_file=0;
end
% title_str='2_18_vs_4_16';
% title_str='10_10_vs_12_8';
% title_str='3_6_vs_6_3';
title_str='3_9_vs_9_3';
% legend_str=sprintf('legend(''single pulse'',''2 18'',''4 16'')');
% legend_str=sprintf('legend(''single pulse'',''10 10'',''12 8'')');
% legend_str=sprintf('legend(''single pulse'',''3 6'',''6 3'')');
legend_str=sprintf('legend(''single pulse'',''3 9'',''9 3'')');
if ~load_file
    numBootstrapTrials=1000;
    comparison_patterns=[4 5]; %not a great way to automate this, so select which cell in the array you want to use for comparison

    load(force_file,'WORD','WORD_N');
    trials_1pulse=WORD{1};
    trials_1pulse_N=WORD_N{1};
    trials1=WORD{comparison_patterns(1)};
    trials2=WORD{comparison_patterns(2)};
    trials1_N=WORD_N{comparison_patterns(1)};
    trials2_N=WORD_N{comparison_patterns(2)};
%     trials_1pulse=WORD{1}(1:120,:);
%     trials_1pulse_N=WORD_N{1}(1:120,:);
%     trials1=WORD{comparison_patterns(1)}(1:120,:);
%     trials2=WORD{comparison_patterns(2)}(1:120,:);
%     trials1_N=WORD_N{comparison_patterns(1)}(1:120,:);
%     trials2_N=WORD_N{comparison_patterns(2)}(1:120,:);
    
    mean_trials_1pulse=mean(trials_1pulse);
    mean_trials_1pulse_N=mean(trials_1pulse_N);
    sterr_trials_1pulse=std(trials_1pulse)/sqrt(size(trials_1pulse,1));
    sterr_trials_1pulse_N=std(trials_1pulse_N)/sqrt(size(trials_1pulse_N,1));
    mean_trials1=mean(trials1);
    mean_trials2=mean(trials2);    
    mean_trials1_N=mean(trials1_N);
    mean_trials2_N=mean(trials2_N);
    sterr_trials1=std(trials1)/sqrt(size(trials1,1));
    sterr_trials2=std(trials2)/sqrt(size(trials2,1));
    sterr_trials1_N=std(trials1_N)/sqrt(size(trials1_N,1));
    sterr_trials2_N=std(trials2_N)/sqrt(size(trials2_N,1));


    timebins=size(trials1,2);

    numOne = size(trials1,1);
    numTwo = size(trials2,1);

    replacedForce1 = zeros(numOne,timebins);
    replacedForce2 = zeros(numTwo,timebins);
    replacedForce1_N = zeros(numOne,timebins);
    replacedForce2_N = zeros(numTwo,timebins);

    comparisonByWord = zeros(numBootstrapTrials,timebins);
    comparisonByWord_N = zeros(numBootstrapTrials,timebins);


    for i = 1:numBootstrapTrials;
        if mod(100*i/numBootstrapTrials,10)==0
            disp([num2str(100*i/numBootstrapTrials) '% done'])
        end
        samplesChosen1 = randsample(1:numOne, numOne, true);
        samplesChosen2 = randsample(1:numTwo, numTwo, true);
        for j = 1: numOne;
            replacedForce1(j, :) = trials1(samplesChosen1(j),:);
            replacedForce1_N(j, :) = trials1_N(samplesChosen1(j),:);
        end
        for j = 1: numTwo;
            replacedForce2(j, :) = trials2(samplesChosen2(j),:);
            replacedForce2_N(j, :) = trials2_N(samplesChosen2(j),:);
        end

        meanOne = mean(replacedForce1);
        stdOne = std(replacedForce1);
        meanTwo = mean(replacedForce2);
        stdTwo = std(replacedForce2);
        
        meanOne_N = mean(replacedForce1_N);
        stdOne_N = std(replacedForce1_N);
        meanTwo_N = mean(replacedForce2_N);
        stdTwo_N = std(replacedForce2_N);

        for j = 1:timebins;
        comparisonByWord(i,j) = ((meanOne(j)-meanTwo(j))/(sqrt((stdOne(j)^2) + (stdTwo(j)^2))));
        comparisonByWord_N(i,j) = ((meanOne_N(j)-meanTwo_N(j))/(sqrt((stdOne_N(j)^2) + (stdTwo_N(j)^2))));
        end
    end

    stdCompByWord = std(comparisonByWord);
    meanCompByWord = mean(comparisonByWord);
    stdCompByWord_N = std(comparisonByWord_N);
    meanCompByWord_N = mean(comparisonByWord_N);
    id=strfind(force_file,'_');
    savename=[force_file(1:id(2)) 'dprime_force_' title_str];
    save(savename,'stdCompByWord','meanCompByWord','stdCompByWord_N','meanCompByWord_N','mean_trials1','mean_trials2',...
        'mean_trials1_N','mean_trials2_N','sterr_trials1','sterr_trials2','sterr_trials1_N','sterr_trials2_N',...
        'mean_trials_1pulse','mean_trials_1pulse_N','sterr_trials_1pulse','sterr_trials_1pulse_N');
else
    load(force_file)
end
fs=20000;
time2=1/fs:1/fs:length(meanCompByWord)/fs;
time2=time2-0.025;
figure
ax(1)=subplot(2,1,1);
% plot(time2,meanCompByWord,'b-',time2,meanCompByWord+stdCompByWord,'b--',time2,meanCompByWord-stdCompByWord,'b--')
% hold on
% plot([min(time2) max(time2)],[0 0],'k-')
% ylabel('d-prime')
x=time2;
y1=mean_trials_1pulse_N+sterr_trials_1pulse_N;
y2=mean_trials_1pulse_N-sterr_trials_1pulse_N;
X=[x fliplr(x)];
Y=[y1 fliplr(y2)];
fill(X,Y,'g')
hold on
y1=mean_trials1_N+sterr_trials1_N;
y2=mean_trials1_N-sterr_trials1_N;
X=[x fliplr(x)];
Y=[y1 fliplr(y2)];
fill(X,Y,'b')
% plot(time2,mean_trials1_N,'b-',time2,mean_trials2_N,'r-')
hold on
y1=mean_trials2_N+sterr_trials2_N;
y2=mean_trials2_N-sterr_trials2_N;
X=[x fliplr(x)];
Y=[y1 fliplr(y2)];
fill(X,Y,'r')
eval(legend_str)

% ylim([0 0.3])

% plot(time2,mean_trials1_N+sterr_trials1_N,'b--',time2,mean_trials1_N-sterr_trials1_N,'b--',...
%     time2,mean_trials2_N+sterr_trials2_N,'r--',time2,mean_trials2_N-sterr_trials2_N,'r--')
plot([min(time2) max(time2)],[0 0],'k-')
ylabel('Normalized Force')
title(strrep(title_str,'_',' '))
ax(2)=subplot(2,1,2);
y1=meanCompByWord_N+stdCompByWord_N;
y2=meanCompByWord_N-stdCompByWord_N;
X=[x fliplr(x)];
Y=[y1 fliplr(y2)];
fill(X,Y,'g')
% plot(time2,meanCompByWord_N,'b-',time2,meanCompByWord_N+stdCompByWord_N,'b--',time2,meanCompByWord_N-stdCompByWord_N,'b--')
hold on
plot([min(time2) max(time2)],[0 0],'k-')
ylabel('d-prime using normalized force')

xlabel('Time relative to stim onset (s)')


linkaxes(ax,'x')


xlim([-0.005 0.06])
% ylim([-2 12])

[y,ind1]=min(abs(time2));
[y,ind2]=min(abs(time2-0.035));
[y,ind]=max(abs(meanCompByWord_N(ind1:ind2)));
max_dprime=meanCompByWord_N(ind1+ind)
std_dprime=stdCompByWord_N(ind1+ind)


end


