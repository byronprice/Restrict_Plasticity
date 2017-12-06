function [] = RestrictSRP_Analysis()
% RestrictSRP_Analysis.m

%  Analyze data from the restricted SRP experiment that Cambria has been
%  running. Briefly, animals view a phase-reversing sinusoidal grating that
%  is restricted to a small region of visual space (20 degree diameter).
%  Training runs for 3 days, then on test day, the animal sees 4 different
%  stimuli: trained position, trained angle; trained position, novel
%  angle; novel position, trained angle; novel position, novel angle. There
%  are two electrode penetrations per animal: target electrode site and
%  non-target electrode site ... the trained position is at the target
%  electrode site.

% Created: 2017/09/06 at 24 Cummington, Boston
%  Byron Price
% Updated: 2017/09/28
%  By: Byron Price
load('RestrictSRPVars.mat','numStimuli','spatFreq','orientation','degreeRadius');

numChans = 2;
filesDay4 = dir('RestrictSRPDataDay4*.mat');

AnimalNames = [];
RFoverlap = [];
VEPmagAtOtherCenter = [];
for ii=1:length(filesDay4)
   fileName4 = filesDay4(ii).name;

   index = regexp(fileName4,'_');
   AnimalName = str2double(fileName4(index+1:end-4));
   cd('~/CloudStation/ByronExp/Retino/');
   fileName = strcat('RetinoMapBayes*',num2str(AnimalName),'.mat');
   files = dir(fileName);
   
   load(files(end).name);
   
   [numChans,~,numSamples] = size(posteriorSample);
   
   x = 1:2:w_pixels;y=1:2:h_pixels;
   [X,Y] = meshgrid(x,y);
   
   meanIms = zeros(length(y),length(x),numChans);
   footprint = zeros(length(y),length(x),numChans);
   result = zeros(numChans,1);
   centerPos = zeros(numChans,2);
   minPos = zeros(numChans,2);
   
   for jj=1:numChans
       N = 1000;
       finalIm = zeros(length(y),length(x),N);
       samples = squeeze(posteriorSample(jj,:,:));
       for ll=1:N
           index = random('Discrete Uniform',numSamples);
           parameterVec = samples(:,index);
           b = [parameterVec(1),parameterVec(4),parameterVec(5),parameterVec(6)];
           distX = X-parameterVec(2);distY = Y-parameterVec(3);
           finalIm(:,:,ll) = b(1)*exp(-(distX.^2)./(2*b(2)*b(2))-...
               (distY.^2)./(2*b(3)*b(3)))+b(4);
       end
       meanIm = mean(finalIm,3);
       meanIms(:,:,jj) = meanIm;
       
       [~,maxInd] = max(meanIm(:));
       [maxrow,maxcol] = ind2sub(size(meanIm),maxInd);

       centerPos(jj,1) = maxrow;centerPos(jj,2) = maxcol;
       
       [~,minInd] = min(meanIm(:));
       [minrow,mincol] = ind2sub(size(meanIm),minInd);
       
       minPos(jj,1) = minrow;minPos(jj,2) = mincol;
       
       datamax = squeeze(finalIm(maxrow,maxcol,:));
       datamin = squeeze(finalIm(minrow,mincol,:));
       
       alpha = 0.05;
       qmax = quantile(datamax,alpha/2);
       qmin = quantile(datamin,1-alpha/2);
%        figure();histogram(datamax);hold on;histogram(datamin);
       result(jj) = (qmax-qmin)>0;
       
       
       temp = permute(finalIm,[3,1,2]);
       temp2 = quantile(temp,alpha/2);
       footprint(:,:,jj) = squeeze(temp2)>qmin;
   end
   
   if result(1) == 1 && result(2) == 1
       AnimalNames = [AnimalNames,AnimalName];
       rf1 = squeeze(footprint(:,:,1));
       rf2 = squeeze(footprint(:,:,2));
       RFoverlap = [RFoverlap;sum((rf1(:) & rf2(:)))/min(sum(rf1(:)),sum(rf2(:)))];
       
       temp = zeros(1,numChans);
       for jj=1:numChans
           otherCenter = centerPos(-jj+3,:);
           currentCenter = centerPos(jj,:);
           currentMin = minPos(jj,:);
           temp(jj) = (meanIms(otherCenter(1),otherCenter(2),jj)-...
               meanIms(currentMin(1),currentMin(2),jj))./(...
               meanIms(currentCenter(1),currentCenter(2),jj)-...
               meanIms(currentMin(1),currentMin(2),jj));
       end
       VEPmagAtOtherCenter = [VEPmagAtOtherCenter;temp];
   end
   
   cd('~/CloudStation/ByronExp/RestrictSRP/');
end

numAnimals = length(AnimalNames);
strobeCodes = cell(numAnimals,4);
VEPs = cell(numAnimals,4);
screenCenter = cell(numAnimals,1);
electrodePositions = cell(numAnimals,1);
targetChannels = zeros(numAnimals,1);
distBetweenChannels = zeros(numAnimals,1);
channelLocations = zeros(numAnimals,numChans,2);
for ii=1:numAnimals
   Day1 = sprintf('RestrictSRPDataDay1_%d.mat',AnimalNames(ii));
   [ChanData,timeStamps,tsevs,svStrobed,sampleFreq] = ExtractSignal(Day1,numChans);
   [Response,eventNames] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed,numChans,numStimuli);
   VEPs{ii,1} = Response;
   strobeCodes{ii,1} = eventNames;
   
   Day2 = sprintf('RestrictSRPDataDay2_%d.mat',AnimalNames(ii));
   [ChanData,timeStamps,tsevs,svStrobed,sampleFreq] = ExtractSignal(Day2,numChans);
   [Response,eventNames] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed,numChans,numStimuli);
   VEPs{ii,2} = Response;
   strobeCodes{ii,2} = eventNames;
   
   Day3 = sprintf('RestrictSRPDataDay3_%d.mat',AnimalNames(ii));
   [ChanData,timeStamps,tsevs,svStrobed,sampleFreq] = ExtractSignal(Day3,numChans);
   [Response,eventNames] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed,numChans,numStimuli);
   VEPs{ii,3} = Response;
   strobeCodes{ii,3} = eventNames;
   
   Day4 = sprintf('RestrictSRPDataDay4_%d.mat',AnimalNames(ii));
   [ChanData,timeStamps,tsevs,svStrobed,sampleFreq] = ExtractSignal(Day4,numChans);
   [Response,eventNames] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed,numChans,numStimuli);
   VEPs{ii,4} = Response;
   strobeCodes{ii,4} = eventNames;
   
   load(sprintf('RestrictSRPStimDay1_%d.mat',AnimalNames(ii)),'w_pixels','h_pixels',...
       'targetChan','centerPositions','DistToScreen');
   conv_factor = 0.2363;
   screenCenter{ii} = [w_pixels/2,h_pixels/4];
   electrodePositions{ii} = centerPositions;
   targetChannels(ii) = targetChan;
   
   positionDegrees = zeros(numChans,2);
   for jj=1:numChans
      xPos = centerPositions(jj,1);
      yPos = centerPositions(jj,2);
      
      relativePos = screenCenter{ii}-[xPos,yPos];
      positionDegrees(jj,:) = 2*atand((relativePos.*conv_factor)./(DistToScreen*10*2));
   end
   temp = sqrt((centerPositions(1,1)-centerPositions(2,1))^2+(centerPositions(1,2)-centerPositions(2,2))^2);
   distBetweenChannels(ii) = 2*atand((temp.*conv_factor)./(DistToScreen*10*2));
%    distBetweenChannels(ii) = sqrt((positionDegrees(1,1)-positionDegrees(2,1)).^2+...
%        (positionDegrees(1,2)-positionDegrees(2,2)).^2);
   channelLocations(ii,:,:) = positionDegrees;
end

% the STROBE CODES ARE AS FOLLOWS
% 0 - grey screen onset, 30 seconds 
% 1 and 2 - flip/flop for trained position and trained orientation ...
%  trained position is the target electrode
% 3 and 4 - flip/flop for trained position and orthogonal orientation
% 5 and 6 - flip/flop for novel position (i.e. at site of non-target
%  electrode) and trained orientation
% 7 and 8 - flip/flop for novel position and orientation

negativityWindow = 60:160;
positivityWindow = 120:250;
compactVEPs = cell(numAnimals,4);
for jj=1:numAnimals
%     figure();
    for mm=1:4
        result = zeros(numChans,1);
        eventNames = strobeCodes{jj,mm};
        for kk=1:numChans
            allVEPs = squeeze(VEPs{jj,mm}(kk,:,:,:));
            temp = squeeze(mean(allVEPs(eventNames==1 | eventNames==2,:,:),1));
            meanVEP = squeeze(mean(temp,1));
            
%             subplot(2,1,kk);
%             plot(meanVEP,'LineWidth',2);axis([0 500 -200 200]);hold on;

            result(kk) = max(meanVEP(positivityWindow))-min(meanVEP(negativityWindow));
        end
        if mm<4
            compactVEPs{jj,mm} = result;
        end
    end
%     subplot(2,1,1);
%     title(sprintf('%d - %3.2f - %3.2f - %d',AnimalNames(jj),distBetweenChannels(jj),...
%         RFoverlap(jj),targetChannels(jj)));
    
end

for jj=1:numAnimals
   result = zeros(numChans,4);
   eventNames = strobeCodes{jj,4};
   for kk=1:numChans
      allVEPs = squeeze(VEPs{jj,4}(kk,:,:,:));
      
      temp1 = squeeze(mean(allVEPs(eventNames==1 | eventNames==2,:,:),1));
      meanVEP = squeeze(mean(temp1,1));
      result(kk,1) = max(meanVEP(positivityWindow))-min(meanVEP(negativityWindow));
      
      temp2 = squeeze(mean(allVEPs(eventNames==3 | eventNames==4,:,:),1));
      meanVEP = squeeze(mean(temp2,1));
      result(kk,2) = max(meanVEP(positivityWindow))-min(meanVEP(negativityWindow));
      
      temp3 = squeeze(mean(allVEPs(eventNames==5 | eventNames==6,:,:),1));
      meanVEP = squeeze(mean(temp3,1));
      result(kk,3) = max(meanVEP(positivityWindow))-min(meanVEP(negativityWindow));
      
      temp4 = squeeze(mean(allVEPs(eventNames==7 | eventNames==8,:,:),1));
      meanVEP = squeeze(mean(temp4,1));
      result(kk,4) = max(meanVEP(positivityWindow))-min(meanVEP(negativityWindow));
   end
   compactVEPs{jj,4} = result;
end

save('RestrictSRPData_AllAnimals.mat','compactVEPs','distBetweenChannels','RFoverlap',...
'channelLocations','targetChannels','AnimalNames','electrodePositions','screenCenter',...
'numAnimals','numChans','VEPmagAtOtherCenter');

% QUESTION 1: When a channel is targeted, is it's VEP bigger? ... on day 1,
%  regardless of orientation

% cannot do normalization based on day 1 VEP ... because there is no day 1
%  VEP for the non-target channel
%  column 1 means the channel was targeted ... column 2 means the channel
%  was not targeted
VEP_when_targeted1 = zeros(numAnimals,2);
distToChan1 = zeros(numAnimals,1);
distToScreenCenter1 = zeros(numAnimals,1);
for ii=1:numAnimals
    targetChan = targetChannels(ii);
    distToChan1(ii) = distBetweenChannels(ii);
    %distToChan1(ii) = VEPmagAtOtherCenter(ii,-targetChan+3);
    distToScreenCenter1(ii) = abs(channelLocations(ii,targetChan,1));

    VEPs = compactVEPs{ii,1};
    for jj=1:numChans
        if jj==targetChan
            VEP_when_targeted1(ii,1) = VEPs(jj);
        else
            VEP_when_targeted1(ii,2) = VEPs(jj);
        end
    end
end
figure();histogram(VEP_when_targeted1(:,1),0:50:400);hold on;
histogram(VEP_when_targeted1(:,2),0:50:400);
title('VEP Magnitudes For Target and Off-Target Channels: Day 1');
xlabel('VEP Magnitude (\muV)');ylabel('Count');
legend('Channel Stimulated','Channel NOT Stimulated');

[~,p] = ttest(VEP_when_targeted1(:,1),VEP_when_targeted1(:,2));
fprintf('\np-value t-test for Target Channel: %3.2e\n',p);
[p2,~] = ranksum(VEP_when_targeted1(:,1),VEP_when_targeted1(:,2));
fprintf('p-value for Wilcoxon rank sum test: %3.2e\n\n',p2);

% QUESTION 2: Does the magnitude of the VEP on the non-target channel
%  depend on the distance between the target and non-target channels (when
%  the target channel is actively stimulated)?
figure();scatter(distToChan1,VEP_when_targeted1(:,1),'LineWidth',2);
hold on;scatter(distToChan1,VEP_when_targeted1(:,2));

for ii=1:numAnimals
   y = linspace(VEP_when_targeted1(ii,1),...
       VEP_when_targeted1(ii,2),50);
   x = ones(50,1).*distToChan1(ii); 
   plot(x,y,'k');
end
title('VEP Magnitudes: Day 1');
xlabel('Distance from Target to Off-Target Channel (DOA)');
ylabel('VEP Magnitude (\muV)');
legend('Channel Stimulated','Channel NOT Stimulated');

figure();scatter(RFoverlap,VEP_when_targeted1(:,2),'LineWidth',2);
title('VEP Magnitudes: Day1');
xlabel('Retinotopic Response Region Overlap');
ylabel('VEP Magnitude (\muV)');

[r3,p3] = corrcoef(distToChan1,VEP_when_targeted1(:,2));
fprintf('p-value for correlation between VEP magnitude and channel distance: %3.2e\n\n',p3(1,2));


VEP_when_targeted2 = zeros(numAnimals*2,2);
distToChan2 = zeros(numAnimals*2,1);
distToScreenCenter2 = zeros(numAnimals*2,1);
for ii=1:numAnimals
    index = (ii-1)*2+1;
    targetChan = targetChannels(ii);
    distToChan2(index:index+1) = distBetweenChannels(ii);

    VEPs = compactVEPs{ii,4};
    for jj=1:numChans
        if jj==targetChan
            VEP_when_targeted2(index,1) = VEPs(jj,2);
            VEP_when_targeted2(index,2) = VEPs(jj,4);
            
            distToScreenCenter2(index) = abs(channelLocations(ii,jj,1));
        else
            VEP_when_targeted2(index+1,1) = VEPs(jj,4);
            VEP_when_targeted2(index+1,2) = VEPs(jj,2);
            distToScreenCenter2(index+1) = abs(channelLocations(ii,jj,1));
        end
    end
end

figure();scatter(distToChan2,VEP_when_targeted2(:,1),'LineWidth',2);
hold on;scatter(distToChan2,VEP_when_targeted2(:,2));
title('VEP Magnitudes: Day 4 Novel');
xlabel('Distance from Target to Off-Target Channel (DOA)');
ylabel('VEP Magnitude (\muV)');
legend('Channel Stimulated','Channel NOT Stimulated');

% QUESTION 3: Does the trained orientation potentiate?
%  i.e. is the VEP at the target electrode for trained position and trained
%  orientation bigger on day 4 than on day 1? 

day1vsDay4fam = zeros(numAnimals,2);
for ii=1:numAnimals
    targetChan = targetChannels(ii);
    day1vsDay4fam(ii,1) = compactVEPs{ii,1}(targetChan);
    day1vsDay4fam(ii,2) = compactVEPs{ii,4}(targetChan,1);
end
% [f,x] = ecdf(day1vsDay4fam(:,2)./day1vsDay4fam(:,1));
% figure();plot(x,f);title('VEP Magnitude Ratio eCDF at Target: Day 4 / Day 1');
figure();histogram(day1vsDay4fam(:,2)./day1vsDay4fam(:,1),0:0.25:3,...
    'FaceColor','c');
title('VEP Magnitude Ratio at Target Channel: Day 4 / Day 1');
xlabel('VEP Magnitude Ratio [Day 4 / Day 1]');
ylabel('Count');

[~,p4] = ttest(day1vsDay4fam(:,2)-day1vsDay4fam(:,1));
forBinoTest = day1vsDay4fam(:,2)./day1vsDay4fam(:,1)>1;
if sum(forBinoTest)<length(forBinoTest)/2
    forBinoTest = ~forBinoTest;
end
p5 = binocdf(sum(forBinoTest),length(forBinoTest),0.5,'upper')*2;
fprintf('p-value t-test Day 1 vs Day 4 ratio: %3.2e\n',p4);
fprintf('p-value binomial test Day 1 vs Day 4: %3.2e\n\n',p5);

% QUESTION 4: At the trained site (the target electrode site), is the VEP
%  for the trained orientation and trained position greater than the VEP 
%  for the novel orientation and trained position? 

trainVsNovelTarget = zeros(numAnimals,2);
for ii=1:numAnimals
    targetChan = targetChannels(ii);
    trainVsNovelTarget(ii,1) = compactVEPs{ii,4}(targetChan,1);
    trainVsNovelTarget(ii,2) = compactVEPs{ii,4}(targetChan,2);
end

figure();subplot(2,1,1);
histogram(trainVsNovelTarget(:,1)./trainVsNovelTarget(:,2),0:0.25:3,...
    'FaceColor','c');
title('VEP Magnitude Ratio at Target Channel: Trained Position, Trained / Novel Orientation');
xlabel('VEP Magnitude Ratio [Trained/Novel]');
ylabel('Count');

[~,p6] = ttest(trainVsNovelTarget(:,1)-trainVsNovelTarget(:,2));
forBinoTest2 = trainVsNovelTarget(:,1)./trainVsNovelTarget(:,2)>1;
p7 = binocdf(sum(forBinoTest2),length(forBinoTest2),0.5,'upper')*2;
fprintf('p-value t-test Train vs Novel at Target: %3.2e\n',p6);
fprintf('p-value binomial test Train vs Novel at Target: %3.2e\n\n',p7);

% QUESTION 5: The same as 4, but for the non-target site (novel position)

trainVsNovelOffTarget = zeros(numAnimals,2);
NovelOnVsNovelOff_OffTarget = zeros(numAnimals,1);
for ii=1:numAnimals
    targetChan = targetChannels(ii);
    trainVsNovelOffTarget(ii,1) = compactVEPs{ii,4}(-targetChan+3,3);
    trainVsNovelOffTarget(ii,2) = compactVEPs{ii,4}(-targetChan+3,4);
    
    NovelOnVsNovelOff_OffTarget(ii) = compactVEPs{ii,4}(-targetChan+3,4)/compactVEPs{ii,4}(-targetChan+3,2);
end

subplot(2,1,2);
histogram(trainVsNovelOffTarget(:,1)./trainVsNovelOffTarget(:,2),0:0.25:3,...
    'FaceColor','m');
title('VEP Magnitude Ratio at Off-Target Channel: Novel Position, Trained / Novel Orientation');
xlabel('VEP Magnitude Ratio [Trained/Novel]');
ylabel('Count');

[~,p8] = ttest(trainVsNovelOffTarget(:,1)-trainVsNovelOffTarget(:,2));
forBinoTest3 = trainVsNovelOffTarget(:,1)./trainVsNovelOffTarget(:,2)<1;
if sum(forBinoTest3)<length(forBinoTest3)/2
    forBinoTest3 = ~forBinoTest3;
end
p9 = binocdf(sum(forBinoTest3),length(forBinoTest3),0.5,'upper')*2;
fprintf('p-value t-test Train vs Novel at Off-Target: %3.2e\n',p8);
fprintf('p-value binomial test Train vs Novel at Off-Target: %3.2e\n\n',p9);

% QUESTION 6: Does difference in VEP magnitude on the non-target electrode
%  between 'novel position and trained orientation' and 'novel position and
%  novel orientation' vary as a function of the distance between the
%  channels?
figure();scatter(maskOverlap,...
    trainVsNovelTarget(:,1)./trainVsNovelTarget(:,2),[],'c','LineWidth',2);
hold on;
scatter(maskOverlap,...
    trainVsNovelOffTarget(:,1)./trainVsNovelOffTarget(:,2),[],'m','^','LineWidth',2);
%legend('Target Channel','Off-Target Channel');
hold on;
for ii=1:numAnimals
   y = linspace(trainVsNovelOffTarget(ii,1)./trainVsNovelOffTarget(ii,2),...
       trainVsNovelTarget(ii,1)./trainVsNovelTarget(ii,2),50);
   x = ones(50,1).*maskOverlap(ii); 
   plot(x,y,'k');
end
x = 0:0.1:(max(maskOverlap)+0.1);
y = ones(length(x),1);
hold on;plot(x,y,'k','LineWidth',2);
title('Trained Divided by Novel VEP Magnitude vs. RR Overlap');

xlabel('Stimulus - Retinotopic Response Region Overlap');
ylabel('VEP Magnitude Ratio [Trained/Novel]');

figure();x = [NovelOnVsNovelOff_OffTarget;NovelOnVsNovelOff_OffTarget];
y = [trainVsNovelOffTarget(:,1)./trainVsNovelOffTarget(:,2);trainVsNovelTarget(:,1)./trainVsNovelTarget(:,2)];
group = cell(numAnimals*2,1);

for ii=1:numAnimals
    group{ii} = 'Off-Target';
end

for ii=numAnimals+1:numAnimals*2
    group{ii} = 'Target';
end
scatterhist(x,y,'Group',group,...
    'Kernel','On','Location','SouthEast','Direction','out',...
    'Color','mc','LineStyle',{'-','-.'},'LineWidth',[3,3],'Marker','od','MarkerSize',[10,10])

% as a function of relative VEP magnitude on off-target channel 
%  due to stimulation at target channel
vepMagMetric = zeros(numAnimals,1);
for ii=1:numAnimals
   targetChan = targetChannels(ii);
   vepMagMetric(ii) = VEPmagAtOtherCenter(ii,targetChan);
   %vepMagMetric(ii) = compactVEPs{ii,1}(-targetChan+3);
end
figure();scatter(vepMagMetric,...
    trainVsNovelTarget(:,1)./trainVsNovelTarget(:,2),[],'c','LineWidth',2);
hold on;
scatter(vepMagMetric,...
    trainVsNovelOffTarget(:,1)./trainVsNovelOffTarget(:,2),[],'m','^','LineWidth',2);
%legend('Target Channel','Off-Target Channel');
hold on;
for ii=1:numAnimals
   y = linspace(trainVsNovelOffTarget(ii,1)./trainVsNovelOffTarget(ii,2),...
       trainVsNovelTarget(ii,1)./trainVsNovelTarget(ii,2),50);
   x = ones(50,1).*vepMagMetric(ii); 
   plot(x,y,'k');
end
x = 0:0.1:(max(vepMagMetric)+0.1);
y = ones(length(x),1);
hold on;plot(x,y,'k','LineWidth',2);
title('Trained Divided by Novel VEP Magnitude (Plasticity Metric)');

xlabel('Relative VEP Magnitude at Off-Target Channel');
ylabel('VEP Magnitude Ratio [Trained/Novel]');

% QUESTION 7: Maybe the distance betwen the channels is not the best
% metric, given that even at a large distance, there may still be some
% overlap between the retinotopic response regions on each electrode. The
% same plot as six but as a function of the magnitude of the potentiation
% on the target channel? As a function of the day 1 VEP magnitude on the
% off-target channel (which would be a proxy for overlap)?

figure();scatter(trainVsNovelTarget(:,1)./trainVsNovelTarget(:,2),...
    trainVsNovelOffTarget(:,1)./trainVsNovelOffTarget(:,2));
title('VEP Magnitude Ratio on Target Channel vs Off-Target Channel');
xlabel('Target Channel');ylabel('Off-Target Channel');

[~,p10] = corrcoef(trainVsNovelTarget(:,1)./trainVsNovelTarget(:,2),...
    trainVsNovelOffTarget(:,1)./trainVsNovelOffTarget(:,2));
fprintf('p-value for correlation between Potentiation on Target and Off-Target: %3.2e\n\n',p10(1,2));


% EXTRA FIGURE ... potentiation across days on each channel
numDays = 4;
data_target = zeros(numAnimals,numDays+1);
data_offtarget = zeros(numAnimals,numDays+1);
for ii=1:numDays
    for jj=1:numAnimals
        target = targetChannels(jj);
        data_target(jj,ii) = compactVEPs{jj,ii}(target,1);
        data_offtarget(jj,ii) = compactVEPs{jj,ii}(-target+3,1);
    end
end

for jj=1:numAnimals
   target = targetChannels(jj);
   data_target(jj,numDays+1) = compactVEPs{jj,numDays}(target,2);
   data_offtarget(jj,numDays+1) = compactVEPs{jj,numDays}(-target+3,2);
end

figure();subplot(2,1,1);hold on;
meanTarget = mean(data_target,1)';
semTarget = std(data_target,[],1)'./sqrt(numAnimals);
meanOffTarget = mean(data_offtarget,1)';
semOffTarget = std(data_offtarget,[],1)'./sqrt(numAnimals);
b = bar([mean(data_target,1)',mean(data_offtarget,1)']);
b(1).FaceColor = 'c';
b(2).FaceColor = 'm';
errorbar((1:numDays+1)-0.15,meanTarget,semTarget,'ko','LineWidth',2);hold on;
errorbar((1:numDays+1)+0.15,meanOffTarget,semOffTarget,'ko','LineWidth',2);
title('VEP Magnitude Across Days');
ylabel('VEP Magnitude (\muV) [mean +/- sem]');
xlabel('Experimental Day (5 is Response on Day 4 to Novel)');
legend('Target','Off-Target');

subplot(2,1,2);hold on;
numConditions = 4;
data_target = zeros(numAnimals,numConditions);
data_offtarget = zeros(numAnimals,numConditions);
for ii=1:numConditions
    for jj=1:numAnimals
        target = targetChannels(jj);
        data_target(jj,ii) = compactVEPs{jj,numDays}(target,ii);
        data_offtarget(jj,ii) = compactVEPs{jj,numDays}(-target+3,ii);
    end
end

meanTarget = mean(data_target,1)';
semTarget = std(data_target,[],1)'./sqrt(numAnimals);
meanOffTarget = mean(data_offtarget,1)';
semOffTarget = std(data_offtarget,[],1)'./sqrt(numAnimals);
b = bar([mean(data_target,1)',mean(data_offtarget,1)']);
b(1).FaceColor = 'c';
b(2).FaceColor = 'm';
errorbar((1:numConditions)-0.15,meanTarget,semTarget,'ko','LineWidth',2);hold on;
errorbar((1:numConditions)+0.15,meanOffTarget,semOffTarget,'ko','LineWidth',2);
title('VEP Magnitude Across Conditions (Test Day)');
ylabel('VEP Magnitude (\muV) [mean +/- sem]');
xlabel('Experimental Condition');
legend('Target','Off-Target');

end

function [ChanData,timeStamps,tsevs,svStrobed,sampleFreq] = ExtractSignal(EphysFileName,numChans)
    % Extract LFP signals from allad, filter, get timestamps
    % read in the .plx file

    load(EphysFileName)

    
    sampleFreq = adfreq;

    Chans = find(~cellfun(@isempty,allad));
    
    dataLength = length(allad{1,Chans(1)});
    ChanData = zeros(dataLength,numChans);
    % lowpass filter the data
    
    preAmpGain = 1;
    for ii=1:numChans
        voltage = 1000.*((allad{1,Chans(ii)}).*SlowPeakV)./(0.5*(2^SlowADResBits)*adgains(Chans(ii))*preAmpGain);
        n = 100;
        lowpass = 100/(sampleFreq/2); % fraction of Nyquist frequency
        blo = fir1(n,lowpass,'low',hamming(n+1));
        temp = filter(blo,1,voltage);
        
        notch = 60/(sampleFreq/2);
        bw = notch/n;
        [b,a] = iirnotch(notch,bw);
        ChanData(:,ii) = filter(b,a,temp);
    end
    
    
    timeStamps = 1/sampleFreq:1/sampleFreq:(dataLength/sampleFreq);

    if length(timeStamps) ~= dataLength
        fprintf('Error: Review allad cell array and timing')
        return;
    end
    
end

function [Response,eventNames] = CollectVEPS(ChanData,timeStamps,tsevs,svStrobed,numChans,numStimuli)
    strobeStart = 33;
    strobeTimes = tsevs{1,strobeStart};
    % COLLECT DATA IN THE PRESENCE OF VISUAL STIMULI
   
    eventNames = unique(svStrobed);
    numEvents = length(eventNames);
    
    stimLen = 500;
    Response = zeros(numChans,numEvents,numStimuli/2,stimLen);

    for ii=1:numChans
        for kk=1:numEvents
            stimStrobes = strobeTimes(svStrobed == eventNames(kk));
            for jj=1:length(stimStrobes)
                stimOnset = stimStrobes(jj);
                [~,index] = min(abs(timeStamps-stimOnset));
                temp = ChanData(index:index+stimLen-1,ii);
                Response(ii,kk,jj,:) = temp;
            end
        end
    end
end