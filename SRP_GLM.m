function [] = SRP_GLM(AnimalName,Day)
% SRP_GLM.m
%  2017/12/04
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
    numARParams = 10;numStimParams = 250;numStimBases = 25;
    
    stimBasisFuns = zeros(numStimParams,numStimBases);
    stdev = (numStimParams/numStimBases)/2;
    centerPoints = linspace(0,numStimParams,numStimBases);
    for ii=1:numStimBases
        temp = exp(-(linspace(0,numStimParams-1,numStimParams)-centerPoints(ii)).^2./(2*stdev*stdev));
        stimBasisFuns(:,ii) = temp;
    end
    
%     numBlockBases = 4;blockMax = 50;
%     numBlockBases = 5;blockMax = max(blockIndex(:,2));
%     blockBasisFuns = zeros(blockMax,numBlockBases);
%     stdev = (blockMax/numBlockBases)/2;
%     centerPoints = linspace(0,blockMax,numBlockBases);
%     for ii=1:numBlockBases
%         temp = exp(-(linspace(0,blockMax-1,blockMax)-centerPoints(ii)).^2./(2*stdev*stdev));
%         blockBasisFuns(:,ii) = temp;
%     end
%     count = 1;stepper = floor(blockMax/numBlockBases);
%     for ii=1:numBlockBases
%         temp = zeros(blockMax,1);
%         temp(count:count+stepper-1) = 1;
%         blockBasisFuns(:,ii) = temp;
%         count = count+stepper;
%     end
    
    basisInds = cell(2,1);
    basisInds{1} = 1:numARParams;
    basisInds{2} = max(basisInds{1})+1:max(basisInds{1})+1+numStimBases-1;
%     basisInds{3} = max(basisInds{2})+1:max(basisInds{2})+1+numBlockBases-1;
    
    [result,classicVEPs] = GetTraining(ChanData,numARParams,numStimParams,...
        stimBasisFuns,timeIndex,numChans,numStimuli,basisInds);
    
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
    numARParams = 10;numStimParams = 250;numStimBases = 25;
    
    stimBasisFuns = zeros(numStimParams,numStimBases);
    stdev = (numStimParams/numStimBases)/2;
    centerPoints = linspace(0,numStimParams,numStimBases);
    for ii=1:numStimBases
        temp = exp(-(linspace(0,numStimParams-1,numStimParams)-centerPoints(ii)).^2./(2*stdev*stdev));
        stimBasisFuns(:,ii) = temp;
    end
    
%     numBlockBases = 6;blockMax = max(squeeze(blockIndex(1,:,2)));
%     blockBasisFuns = zeros(blockMax,numBlockBases);
%     stdev = (blockMax/numBlockBases)/2;
%     centerPoints = linspace(0,blockMax,numBlockBases);
%     for ii=1:numBlockBases
%         temp = exp(-(linspace(0,blockMax-1,blockMax)-centerPoints(ii)).^2./(2*stdev*stdev));
%         blockBasisFuns(:,ii) = temp;
%     end
    
    basisInds = cell(5,1);
    basisInds{1} = 1:numARParams;
    
    basisInds{2} = max(basisInds{1})+1:max(basisInds{1})+1+numStimBases-1;
    basisInds{3} = max(basisInds{2})+1:max(basisInds{2})+1+numStimBases-1;
    basisInds{4} = max(basisInds{3})+1:max(basisInds{3})+1+numStimBases-1;
    basisInds{5} = max(basisInds{4})+1:max(basisInds{4})+1+numStimBases-1;
    
%     basisInds{6} = max(basisInds{5})+1:max(basisInds{5})+1+numBlockBases-1;
%     basisInds{7} = max(basisInds{6})+1:max(basisInds{6})+1+numBlockBases-1;
%     basisInds{8} = max(basisInds{7})+1:max(basisInds{7})+1+numBlockBases-1;
%     basisInds{9} = max(basisInds{8})+1:max(basisInds{8})+1+numBlockBases-1;
%     
%     basisInds{10} = max(basisInds{9})+1;
%     basisInds{11} = max(basisInds{10})+1;
%     basisInds{12} = max(basisInds{11})+1;
%     basisInds{13} = max(basisInds{12})+1;
    
    [result,classicVEPs] = GetTesting(ChanData,numARParams,numStimParams,...
        stimBasisFuns,timeIndex,numChans,numStimuli,basisInds);
end
% figure();plot(classicVEP(:,1),'LineWidth',2);hold on;
% plot(classicVEP(:,2),'LineWidth',2);
% legend('LH','RH');

filename = sprintf('RestrictSRPResults-Day%d_%d.mat',Day,AnimalName);
save(filename,'result','classicVEPs','stimBasisFuns',...
    'timeIndex','blockIndex','basisInds','DayType','DistToScreen','targetChan',...
    'Radius','numStimuli');
end

function [result,classicVEP] = GetTraining(ChanData,numARParams,numStimParams,...
    stimBasisFuns,timeIndex,numChans,numStimuli,basisInds)
    timeIndex(:,1) = timeIndex(:,1)-numARParams;
    result = struct('b',cell(numChans,3),'dev',zeros(numChans,3),...
        'se',cell(numChans,3),'Fstat',cell(numChans,3),'Fpval',cell(numChans,3));
    numStimBases = size(stimBasisFuns,2);
    N = numStimuli*numStimParams+numARParams;
    for ii=1:numChans
        y = ChanData(numARParams+1:end,ii);
        totParams = numARParams+numStimBases;
        design = zeros(length(y),totParams);
        count = 1;
        for jj=numARParams:-1:1
            design(:,count) = ChanData(jj:end-count,ii);
            count = count+1;
        end
        
        stimBasisDesign = zeros(length(y),numStimParams);
        % blockNum = zeros(length(y),1);
%         blockBasisDesign = zeros(length(y),numBlockBases);
        for kk=1:numStimuli
            for jj=1:numStimParams
                stimBasisDesign(timeIndex(kk,1)+jj,jj) = 1;
                %blockNum(timeIndex(kk,1)+jj) = blockIndex(kk,1);
%                 blockBasisDesign(timeIndex(kk,1)+jj,blockIndex(kk,1)) = 1;
            end
        end
        
        design(:,basisInds{2}) = stimBasisDesign*stimBasisFuns;
%         design(:,basisInds{3}) = blockBasisDesign;
        %design(:,basisInds{4}) = blockNum;
        [b,dev,stats] = glmfit(design,y,'normal');

        result(ii,1).b = b;
        result(ii,1).dev = dev;
        result(ii,1).se = stats.se;
        
        allInds = 1:totParams;
        [~,indsToExclude,~] = intersect(allInds,basisInds{1});
        newInds = allInds;
        newInds(indsToExclude) = [];
        [b,dev,stats] = glmfit(design(:,newInds),y,'normal');
        
        F = ((dev-result(ii,1).dev)/(length(indsToExclude)))/(result(ii,1).dev/(length(y)-totParams-2));
        result(ii,2).b = b;
        result(ii,2).dev = dev;
        result(ii,2).se = stats.se;
        result(ii,2).Fstat = F;
        result(ii,2).Fpval = fcdf(F,length(indsToExclude),length(y)-(totParams+1),'upper');
        
        for jj=2
            [~,indsToExclude,~] = intersect(allInds,basisInds{jj});
            newInds = allInds;
            newInds(indsToExclude) = [];
            [b,dev,stats] = glmfit(design(:,newInds),y,'normal');
            
            F = ((dev-result(ii,1).dev)/(length(indsToExclude)))/(result(ii,1).dev/(N-totParams-2));
            result(ii,jj+1).b = b;
            result(ii,jj+1).dev = dev;
            result(ii,jj+1).se = stats.se;
            result(ii,jj+1).Fstat = F;
            result(ii,jj+1).Fpval = fcdf(F,length(indsToExclude),N-(totParams+1),'upper');
        end
    end
    
    classicVEP = zeros(numStimParams,numChans);
    for jj=1:numChans
        for ii=1:numStimuli
            classicVEP(:,jj) = classicVEP(:,jj)+ChanData(timeIndex(ii,1)+numARParams...
                :timeIndex(ii,1)+numARParams+numStimParams-1,jj)./numStimuli;
        end
    end


end


function [result,classicVEP] = GetTesting(ChanData,numARParams,numStimParams,...
    stimBasisFuns,timeIndex,numChans,numStimuli,basisInds)

    timeIndex(:,:,1) = timeIndex(:,:,1)-numARParams;
    result = struct('b',cell(numChans,6),'dev',zeros(numChans,6),...
        'se',cell(numChans,6),'Fstat',cell(numChans,6),'Fpval',cell(numChans,6));
    numStimBases = size(stimBasisFuns,2);
    N = numStimuli*numStimParams+numARParams;
    for ii=1:numChans
        y = ChanData(numARParams+1:end,ii);
        totParams = numARParams+4*numStimBases;
        design = zeros(length(y),totParams);
        count = 1;
        for jj=numARParams:-1:1
            design(:,count) = ChanData(jj:end-count,ii);
            count = count+1;
        end
        
        for ll=1:4
            stimBasisDesign = zeros(length(y),numStimParams);
%             blockNum = zeros(length(y),1);
%             blockBasisDesign = zeros(length(y),blockMax);
            for kk=1:numStimuli
                for jj=1:numStimParams
                    stimBasisDesign(timeIndex(ll,kk,1)+jj,jj) = 1;
%                     blockNum(timeIndex(ll,kk,1)+jj) = blockIndex(ll,kk,1);
%                     blockBasisDesign(timeIndex(ll,kk,1)+jj,blockIndex(ll,kk,2)) = 1;
                end
            end
            
            design(:,basisInds{1+ll}) = stimBasisDesign*stimBasisFuns;
%             design(:,basisInds{5+ll}) = blockBasisDesign*blockBasisFuns;
%             design(:,basisInds{9+ll}) = blockNum;
        end
        
        [b,dev,stats] = glmfit(design,y,'normal');

        result(ii,1).b = b;
        result(ii,1).dev = dev;
        result(ii,1).se = stats.se;
        
        allInds = 1:totParams;
        [~,indsToExclude,~] = intersect(allInds,basisInds{1});
        newInds = allInds;
        newInds(indsToExclude) = [];
        [b,dev,stats] = glmfit(design(:,newInds),y,'normal');
        
        F = ((dev-result(ii,1).dev)/(length(indsToExclude)))/(result(ii,1).dev/(length(y)-totParams-2));
        result(ii,2).b = b;
        result(ii,2).dev = dev;
        result(ii,2).se = stats.se;
        result(ii,2).Fstat = F;
        result(ii,2).Fpval = fcdf(F,length(indsToExclude),length(y)-(totParams+1),'upper');
        
        for jj=2:5
            [~,indsToExclude,~] = intersect(allInds,basisInds{jj});
            newInds = allInds;
            newInds(indsToExclude) = [];
            [b,dev,stats] = glmfit(design(:,newInds),y,'normal');
            
            F = ((dev-result(ii,1).dev)/(length(indsToExclude)))/(result(ii,1).dev/(N-totParams-2));
            result(ii,jj+1).b = b;
            result(ii,jj+1).dev = dev;
            result(ii,jj+1).se = stats.se;
            result(ii,jj+1).Fstat = F;
            result(ii,jj+1).Fpval = fcdf(F,length(indsToExclude),N-(totParams+1),'upper');
        end
    end
    
    classicVEP = cell(4,numChans);
    for ll=1:4
        for jj=1:numChans
            classicVEP{ll,jj} = zeros(numStimParams,1);
            for ii=1:numStimuli
                classicVEP{ll,jj} = classicVEP{ll,jj}+ChanData(timeIndex(ll,ii,1)+numARParams...
                    :timeIndex(ll,ii,1)+numARParams+numStimParams-1,jj)./numStimuli;
            end
        end
    end
end

% % BAYESIAN APPROACH
% X = design;clear design;
% N = size(X,1);
% 
% numParams = numARParams+numStimParams+4;
% numIter = 1e6;burnIn = 1e5;skipRate = 500;
% params = zeros(numParams,numIter);
% posteriorProb = zeros(numIter,1);
% 
% priorMu = 0;
% alpha = 1;beta = 0.9;delta = 1.1;gama = 1;
% abprior1=1e-3;abprior2=1e-3;
% 
% params(1:numARParams+numStimParams+1,1) = ...
%     mvnrnd(zeros(numARParams+numStimParams+1,1),eye(numARParams+numStimParams+1))';
% params(end-3,1) = log(alpha);
% params(end-2,1) = log(beta);
% params(end-1,1) = log(gama);
% params(end,1) = log(delta);
% 
% logNormalPDF = @(y,mu,variance) -0.5*log(variance)-(y-mu).^2./(2*variance);
% logGammaPDF = @(x,a,b) a*log(b)-log(gamma(a))+(a-1).*log(x)-x.*b;
% 
% tempMu = X*params(1:end-4,1);
% variance = (1/N)*sum((tempMu-y).^2);
% prevLogLikelihood = sum(logNormalPDF(y,tempMu,variance));
% 
% smoothPrior = del2(params(numARParams+2:end-4,1));
% prevLogPrior = ((numARParams)*0.5)*log(beta)-0.5*beta*...
%     (params(2:numARParams+1,1)'*params(2:numARParams+1,1))+...
%             0.5*log(alpha)-0.5*alpha*(params(1,1)-priorMu)^2+...
%             ((numStimParams)*0.5)*log(gama)-0.5*gama*...
%             (params(numARParams+2:end-4,1)'*params(numARParams+2:end-4,1))+...
%             ((numStimParams)*0.5)*log(delta)-0.5*delta*(smoothPrior'*smoothPrior)+...
%             sum(logGammaPDF([alpha,beta,gama,delta],abprior1,abprior2));
% 
% posteriorProb(1) = prevLogLikelihood+prevLogPrior;
% 
% % FOR AUTOMATIC CREATION OF UPDATE MATRIX
% updateParam = logspace(-0.3,-2,burnIn);
% loglambda = ones(numParams,1).*log(2.38^2);
% updateMu = zeros(numParams,1);
% updateMu(1:numARParams+numStimParams+1) = ...
%     mvnrnd(zeros(numARParams+numStimParams+1,1),eye(numARParams+numStimParams+1))';
% updateMu(end-3) = log(1.5);updateMu(end-2) = log(1.4);
% updateMu(end-1) = log(1.3);updateMu(end) = log(1.2);
% optimalAccept = 0.234;
% 
% q = numParams;qIdentity = eye(q);identity = eye(numParams);
% W = normrnd(0,1,[numParams,q]);
% 
% % sigmasquare = 1.5;expectedMean = 0.7;sigmamean = 1;
% M = W'*W;%+sigmasquare.*qIdentity;
% eigenvals = M(1:q+1:end)';
% % halfCov = cholcov(M);
% halfSigma = cholcov(eye(numParams));
% 
% p = ones(q,1)./q;
% % figure(2);scatter(1,posteriorProb(1));hold on;pause(0.1);
% for ii=2:burnIn
%     index = find(mnrnd(1,p)==1);
%     lambda = loglambda(index);
%     stdev = sqrt(exp(lambda).*eigenvals(index));
%     pStar = params(:,ii-1)+W(:,index)*normrnd(0,stdev);
%     
% %     if sum(pStar(end-1:end)<=-100) == 0
%         
%         tempMu = X*pStar(1:end-4);
%         variance = (1/N)*sum((tempMu-y).^2);
%         pStarLogLikelihood = sum(logNormalPDF(y,tempMu,variance));
%         
%         smoothPrior = del2(pStar(numARParams+2:end-4));
%         alpha = exp(pStar(end-3));beta = exp(pStar(end-2));gama = exp(pStar(end-1));
%         delta = exp(pStar(end));
%         pStarLogPrior = ((numARParams)*0.5)*log(beta)-0.5*beta*...
%             (pStar(2:numARParams+1)'*pStar(2:numARParams+1))+...
%             0.5*log(alpha)-0.5*alpha*(pStar(1)-priorMu)^2+...
%             ((numStimParams)*0.5)*log(gama)-0.5*gama*...
%             (pStar(numARParams+2:end-4)'*pStar(numARParams+2:end-4))+...
%             ((numStimParams)*0.5)*log(delta)-0.5*delta*(smoothPrior'*smoothPrior)+...
%             sum(logGammaPDF([alpha,beta,gama,delta],abprior1,abprior2));
% 
%         logA = (pStarLogLikelihood+pStarLogPrior)-posteriorProb(ii-1);
%         
%         if log(rand) < logA
%             params(:,ii) = pStar;
%             posteriorProb(ii) = pStarLogLikelihood+pStarLogPrior;
%         else
%             params(:,ii) = params(:,ii-1);
%             posteriorProb(ii) = posteriorProb(ii-1);
%         end
%         
%         if mod(ii,150) == 0
%             meanSubtract = params(:,ii)-updateMu;
%             updateMu = updateMu+updateParam(ii).*meanSubtract;
%             %         sigma2 = sigma2+updateParam(ii).*(meanSubtract*meanSubtract'-sigma2);
%             halfSigma = halfSigma+updateParam(ii).*(triu((inv(halfSigma))*(halfSigma'*halfSigma+meanSubtract*...
%                 meanSubtract')*((inv(halfSigma))')-identity)-halfSigma);
%             sigma = halfSigma'*halfSigma;
%             
%             Z = inv(tril(W'*W)')'*W';
%             W = sigma*Z'*inv(triu(Z*sigma*Z'));
%             W = normc(W);
%             eigenvals = diag(W'*sigma*W);
%             lambda = lambda+updateParam(ii).*(exp(min(0,logA))-optimalAccept);
%             
%         end
%         
%     loglambda(index) = lambda;
% %     scatter(ii,posteriorProb(ii));hold on;pause(0.01);
% %     error(ii) = mean(abs([log(baseRate);historyB]-updateMu));
% %     scatter(ii,error);hold on;pause(0.01);
% end
% 
% [V,D] = eig(cov(params(:,1e4:10:burnIn)'));
% W = V*sqrtm(D);
% eigenvals = diag(W'*W);
% 
% tempW = [];
% tempEigs = [];
% for jj=1:numParams
%    if eigenvals(jj) > 1e-6
%        tempW = [tempW,W(:,jj)];
%        tempEigs = [tempEigs,eigenvals(jj)];
%    end
% end
% 
% W = fliplr(tempW);
% eigenvals = fliplr(tempEigs);
% q = length(eigenvals);
% p = ones(q,1)./q;
% updateParam = 1e-2;
% 
% acceptRate = 0;
% for ii=burnIn+1:numIter
%     index = find(mnrnd(1,p)==1);
%     lambda = loglambda(index);
%     stdev = sqrt(exp(lambda).*eigenvals(index));
%     pStar = params(:,ii-1)+W(:,index)*normrnd(0,stdev);
%     
% %     if sum(pStar(end-1:end)<=-100) == 0
%         tempMu = X*pStar(1:end-4);
%         variance = (1/N)*sum((tempMu-y).^2);
%         pStarLogLikelihood = sum(logNormalPDF(y,tempMu,variance));
%         
%         smoothPrior = del2(pStar(numARParams+2:end-4));
%         alpha = exp(pStar(end-3));beta = exp(pStar(end-2));gama = exp(pStar(end-1));
%         delta = exp(pStar(end));
%         pStarLogPrior = ((numARParams)*0.5)*log(beta)-0.5*beta*...
%             (pStar(2:numARParams+1)'*pStar(2:numARParams+1))+...
%             0.5*log(alpha)-0.5*alpha*(pStar(1)-priorMu)^2+...
%             ((numStimParams)*0.5)*log(gama)-0.5*gama*...
%             (pStar(numARParams+2:end-4)'*pStar(numARParams+2:end-4))+...
%             ((numStimParams)*0.5)*log(delta)-0.5*delta*(smoothPrior'*smoothPrior)+...
%             sum(logGammaPDF([alpha,beta,gama,delta],abprior1,abprior2));
% 
%         logA = (pStarLogLikelihood+pStarLogPrior)-posteriorProb(ii-1);
%         
%         
%         if log(rand) < logA
%             params(:,ii) = pStar;
%             posteriorProb(ii) = pStarLogLikelihood+pStarLogPrior;
%             acceptRate = acceptRate+1;
%         else
%             params(:,ii) = params(:,ii-1);
%             posteriorProb(ii) = posteriorProb(ii-1);
%         end
%         lambda = lambda+updateParam.*(exp(min(0,logA))-optimalAccept);
% %     else
% %         params(:,ii) = params(:,ii-1);
% %         posteriorProb(ii) = posteriorProb(ii-1);
% %         lambda = lambda+updateParam.*(-optimalAccept);
% %     end
%     loglambda(index) = lambda;
% end
% 
% % figure();plot(error);
% fprintf('Final Acceptance Rate: %3.2f\n',acceptRate/(numIter-burnIn-1));
% posteriorSamples = params(:,burnIn+1:skipRate:end);
% % figure();histogram(posteriorSamples(2,:));
% 
% [~,ind] = max(posteriorProb);
% MAP = params(:,ind);
% posteriorMean = mean(posteriorSamples,2);
% posteriorMedian = median(posteriorSamples,2);