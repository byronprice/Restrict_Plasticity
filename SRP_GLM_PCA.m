% SRP_GLM_PCA.m
% numChans = 2;numARParams = 10;numStimParams = 250;
% histData = [];stimData = [];
% files = dir('RestrictSRPResults*.mat');
% for ii=1:length(files)
% filename = files(ii).name;
% load(filename);
%    for jj=1:numChans
%        b = result(jj,1).b;
%        histData = [histData,b(basisInds{1}+1)];
%        for kk=3:size(result,2)
%            if result(jj,kk).Fpval <= 0.05
%                stimData = [stimData,stimBasisFuns*b(basisInds{kk-1}+1)];
%            end
%        end
%    end
% end
% 
% numIter = 1000;
% bootStrapSamples = zeros(numIter,numStimParams);
% 
% for ii=1:1000
%     ind = randperm(size(stimData,2),1);
%     bootStrapSamples(ii,:) = stimData(:,ind)'; 
% end
% 
% q = quantile(bootStrapSamples,[0.05/2,0.5,1-0.05/2],1);
% 
% figure;boundedline(1:250,q(2,:)',[q(2,:)'-q(1,:)',q(3,:)'-q(2,:)']);
% title('Mean Visually-Evoked Potential with 95% Confidence Bound');
% 
% [Whist,eigsHist,sigmaHist,muHist,qHist] = BayesianPCA(histData,5);
% [Wstim,eigsStim,sigmaStim,muStim,qStim] = BayesianPCA(stimData,50);
% 
% fprintf('History Dimensions: %d\n',qHist);
% fprintf('Stimulus Dimensions: %d\n',qStim);
% 
% save('RestrictSRP-FinalPCA.mat','Whist','eigsHist','sigmaHist','muHist','qHist',...
%     'Wstim','eigsStim','sigmaStim','muStim','qStim','histData','stimData');
% 
% for ii=1:qHist
%    figure;plot(Whist(:,ii)); 
%    title(sprintf('History-Dependence PC-%d',ii));
% end
% 
% for ii=1:qStim
%    figure;plot(Wstim(:,ii)); 
%    title(sprintf('Stimulus PC-%d',ii));
% end

load('RestrictSRP-FinalPCA.mat','Whist','qHist','Wstim','qStim');
for ii=1:qHist
    Whist(:,ii) = Whist(:,ii)./norm(Whist(:,ii));
end

for ii=1:qStim
   Wstim(:,ii) = Wstim(:,ii)./norm(Wstim(:,ii)); 
end

% qStim = 6;  ... do we need PCs that have an effect before 50ms latency to
%                   V1??
% Wstim = Wstim(:,1:qStim);

files = dir('RestrictSRPReduce-Day4*.mat');

numConditions = 4;numParams = qStim+1+qHist;
data_Target = zeros(numConditions,17,numParams);
data_Off = zeros(numConditions,17,numParams);

numChans = 2;
for ii=1:17
    index = regexp(files(ii).name,'_');
    Day = str2double(files(ii).name(index-1));
    
    load(files(ii).name,'numARParams','numStimParams','numStimuli','reduceData','targetChan')
    totParams = 1+qHist+qStim;
    histInds = 2:qHist+1;
    stimInds = histInds(end)+1:totParams;
    for jj=1:numChans
        for zz=1:numConditions
            data = reduceData{jj,zz};
            temp = data(:,numARParams+1:end)';
            y = temp(:);N = length(y);
            design = zeros(numStimuli*numStimParams,totParams);
            design(:,1) = 1;
            
            histDesign = zeros(numStimuli*numStimParams,numARParams);
            count = 1;
            for kk=1:numStimuli
                for ll=1:numStimParams
                    for mm=1:numARParams
                        histDesign(count,mm) = data(kk,numARParams+ll-mm);
                    end
                    count = count+1;
                end
            end
            
            stimDesign = zeros(length(y),numStimParams);
            count = 1;
            for kk=1:numStimuli
                for ll=1:numStimParams
                    stimDesign(count,ll) = 1;
                    count = count+1;
                end
            end
            
            design(:,histInds) = histDesign*Whist;
            design(:,stimInds) = stimDesign*Wstim;
            
            [b,dev,~] = glmfit(design,y,'normal','constant','off');
            
            [~,dev2,~] = glmfit(design(:,1:max(histInds)),y,'normal','constant','off');
            
            F = ((dev2-dev)/qStim)/(dev/(N-totParams-1));
            Fpval = fcdf(F,qStim,N-totParams,'upper');
            fprintf('F: %3.2e\n',F);
            fprintf('F-test p-value: %3.2e\n',Fpval);
            
            if jj==targetChan
                data_Target(zz,ii,:) = b(:);
            else
                data_Off(zz,ii,:) = b(:);
            end
        end
    end
end

save('RestrictSRP_Day4PCAResults.mat','Whist','Wstim','qHist','qStim','histInds',...
    'stimInds','data_Target','data_Off');