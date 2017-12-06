function [] = GetReduce(AnimalName,Day)
% GetReduce.m
%  Smaller representation of each dataset

filename = sprintf('RestrictSRPDataDay%d_%d.mat',Day,AnimalName);
load(filename);

stimFile = sprintf('RestrictSRPStimDay%d_%d.mat',Day,AnimalName);
load(stimFile,'DayType','DistToScreen','targetChan','Radius','numStimuli');

numChans = 2;
sampleFreq = adfreq;

Chans = find(~cellfun(@isempty,allad));

dataLength = length(allad{1,Chans(1)});
ChanData = zeros(dataLength,numChans);
% lowpass filter the data

preAmpGain = 1;
for ii=1:numChans
    voltage = 1000.*((allad{1,Chans(ii)}).*SlowPeakV)./(0.5*(2^SlowADResBits)*adgains(Chans(ii))*preAmpGain);
    n = 3;
    
    notch = 60/(sampleFreq/2);
    bw = notch/n;
    [b,a] = iirnotch(notch,bw);
    ChanData(:,ii) = filtfilt(b,a,voltage);
end


timeStamps = 1/sampleFreq:1/sampleFreq:(dataLength/sampleFreq);

if strcmp(DayType,'train') == 1

    eventTimes = tsevs{33};
    eventTimes = eventTimes(svStrobed==1 | svStrobed==2);
    timeIndex = zeros(numStimuli,2);
    blockIndex = ones(numStimuli,2);
    
    currentBlock = 1;
    currentBlockIndex = 1;
    for ii=1:numStimuli
        [~,timeIndex(ii,1)] = min(abs(eventTimes(ii)-timeStamps));
        timeIndex(ii,2) = 1;
        if ii>1
            difference = timeIndex(ii,1)-timeIndex(ii-1,1);
            if difference>sampleFreq
                currentBlock = currentBlock+1;
                currentBlockIndex = 1;
            end
        end
        blockIndex(ii,1) = currentBlock; % 1 to 4 (which block of stimuli)
        blockIndex(ii,2) = currentBlockIndex; % 1 to 50, which presentation within block
        currentBlockIndex = currentBlockIndex+1;
    end
    numARParams = 10;numStimParams = 250;
    
    reduceData = cell(numChans,1);
    for ii=1:numChans
        temp = zeros(numStimuli,numStimParams+numARParams);
        for jj=1:numStimuli
            temp(jj,:) = ChanData(timeIndex(jj,1)-numARParams...
                    :timeIndex(jj,1)+numStimParams-1,ii);
        end
        reduceData{ii} = temp;
    end
    
elseif strcmp(DayType,'test') == 1
    
    eventTimes = tsevs{33};
    eventTimes1 = eventTimes(svStrobed==1 | svStrobed==2);
    eventTimes2 = eventTimes(svStrobed==3 | svStrobed==4);
    eventTimes3 = eventTimes(svStrobed==5 | svStrobed==6);
    eventTimes4 = eventTimes(svStrobed==7 | svStrobed==8);
    
    eventTimes = [eventTimes1(:),eventTimes2(:),eventTimes3(:),eventTimes4(:)];
    
    timeIndex = zeros(4,numStimuli,2);
    blockIndex = ones(4,numStimuli,2);
    
    for jj=1:4
        currentBlock = 1;
        currentBlockIndex = 1;
        for ii=1:numStimuli
            [~,timeIndex(jj,ii,1)] = min(abs(eventTimes(ii,jj)-timeStamps));
            timeIndex(jj,ii,2) = 1;
            if ii>1
                difference = timeIndex(jj,ii,1)-timeIndex(jj,ii-1,1);
                if difference>sampleFreq
                    currentBlock = currentBlock+1;
                    currentBlockIndex = 1;
                end
            end
            blockIndex(jj,ii,1) = currentBlock; % 1 to 4 (which block of stimuli)
            blockIndex(jj,ii,2) = currentBlockIndex; % 1 to 50, which presentation within block
            currentBlockIndex = currentBlockIndex+1;
        end
    end
    numARParams = 10;numStimParams = 250;
    
    reduceData = cell(numChans,4);
    for kk=1:4
        for ii=1:numChans
            temp = zeros(numStimuli,numStimParams+numARParams);
            for jj=1:numStimuli
                temp(jj,:) = ChanData(timeIndex(kk,jj,1)-numARParams...
                    :timeIndex(kk,jj,1)+numStimParams-1,ii);
            end
            reduceData{ii,kk} = temp;
        end
    end
end

filename = sprintf('RestrictSRPReduce-Day%d_%d.mat',Day,AnimalName);
save(filename,'reduceData','numARParams','numStimParams',...
    'timeIndex','blockIndex','DayType','DistToScreen','targetChan',...
    'Radius','numStimuli');

end