function optOut = bayesOptFromC(N,resultsPath, workingPath, outPath)
%get the next N parallel evaluation points from the model

%% Define acceptable parameters
nFeat = optimizableVariable('nFeat',[5 100],'Type','integer');
nFolds = optimizableVariable('nFolds',[2 10],'type','integer');
nLambda = optimizableVariable('nLambda',[5 50],'type','integer');
alpha = optimizableVariable('alpha',[0 1],'type','real');
nRepeats = optimizableVariable('nRepeats',[5 20],'type','integer');
maxDepth = optimizableVariable('maxDepth',[2 8],'type','integer');
minSize = optimizableVariable('minSize',[5 100],'type','integer','transform','log');

%constrained to treeStopCount <= minTrees <= maxTrees,
maxTrees = optimizableVariable('maxTrees',[25 100],'type','integer','transform','log');
minTrees = optimizableVariable('minTrees',[20 100],'type','integer','transform','log');
treeStopCount = optimizableVariable('treeStopCount',[10 50],'type','integer');

treeStopThresh = optimizableVariable('treeStopThresh',[.01 .5],'type','real','transform','log');

ensembleSize = optimizableVariable('ensembleSize',[32 100],'type','integer');

%these variables are constrained to sum <= 1
probNeg = optimizableVariable('probNeg',[.1 .9],'type','real');
probPos = optimizableVariable('probPos',[.1 .9],'type','real');

allParams = [nFeat nFolds nLambda alpha nRepeats maxDepth minSize maxTrees minTrees treeStopCount treeStopThresh ensembleSize probNeg probPos ];
allParamsStrings = {'nFeat','nFolds','nLambda','alpha','nRepeats','maxDepth','minSize','maxTrees','minTrees','treeStopCount','treeStopThresh','ensembleSize','probNeg','probPos'};

%allParamsStrings = {'nFeat','nRepeats','nLambda','alpha','nFolds','maxDepth','minSize','maxTrees','minTrees','treeStopCount','treeStopThresh','ensembleSize','probNeg','probPos'};

%% Load results
k = ls(sprintf('%s\\*.in', resultsPath));
s = struct();
trainTime = zeros(size(k,1),1);
wloss = zeros(size(k,1),27);
for l = 1:14
    s.(allParamsStrings{l}) = 0;
end

for kk = 1:size(k,1)
    t = struct();
    f = fopen(sprintf('%s\\%s',resultsPath,k(kk,:)),'r');
    for l = 1:14
        t.(allParamsStrings{l}) = str2num(fgets(f));
    end
    for l = 1:t.ensembleSize+1
        fgets(f); %discard the ecoc matrix
    end
    trainTime(kk) = str2num(fgets(f));
    for l = 1:9
        fgets(f); %discard the transform indicators
    end
    for l=1:27
        temp = split(fgets(f),',');
        wloss(kk,l) = str2num(temp{4});
    end
    
    fclose(f);
    s(kk) = t;
end
wloss= reshape(wloss,size(k,1),3,3,3);
loss = mean(wloss(:,[6,8,12,16,20,22]),2);

%% Load working data
k2 = ls(sprintf('%s\\*.in',workingPath));
for kk = 1:size(k2,1)
    t = struct();
    f = fopen(sprintf('%s\\%s',workingPath,k2(kk,:)),'r');
    for l = 1:14
        t.(allParamsStrings{l}) = str2num(fgets(f));
    end
    fclose(f);
    s(kk+size(k,1)) = t;
end

initX = struct2table(s(1:size(k,1)));
fprintf('Completed %d paramsets so far of %d assigned.\n', size(k,1), size(k,1)+ size(k2,1));
%% Determine next samples
optOut = bayesopt(@(x) x, allParams', 'initialX',initX, 'initialObjective', loss, 'initialObjectiveEvaluationTimes', trainTime, 'maxObjectiveEvaluations',1, 'xconstraintfcn',@doConstraints, 'verbose',0, 'plotfcn',[]);

%we will assume that our prediction of the in-progress points is accurate,
%subject to the contraint it must fall above the minimum observed objective
working = struct2table(s(size(k,1)+1 : end));
pred = optOut.predictObjective(working);

minObs = min(optOut.ObjectiveTrace(optOut.FeasibilityTrace));
pred(pred<minObs) = minObs;
% pred(pred<optOut.MinEstimatedObjective) = optOut.MinEstimatedObjective;

predTime = optOut.predictObjectiveEvaluationTime(working);

if isempty(k2)
    if isempty(k)
        allFileInds = [];
    else
        allFileInds = str2num(k(:,1:4));
    end
else
    allFileInds = cat(1,str2num(k(:,1:4)),str2num(k2(:,1:4)));
end
loss = cat(1,loss,pred);
trainTime = cat(1,trainTime,predTime);
s = struct2table(s);
%add results to model and resume

if isempty(allFileInds)
    optOut = bayesopt(@(x) x, allParams', 'maxObjectiveEvaluations', N, 'NumSeedPoints', N, 'xconstraintfcn', @doConstraints, 'verbose', 0, 'plotfcn', []);
    for n=1:N
        a = makeparamsfile(table2array(optOut.XTrace(n,:)), allFileInds, outPath);
        allFileInds = cat(1, allFileInds, a);
    end
    return
end

for n=1:N
    %re-run the model with the assumed data
    
    optOut = bayesopt(@(x) x, allParams', 'initialX',s, 'initialObjective', loss, 'initialObjectiveEvaluationTimes', trainTime, 'maxObjectiveEvaluations',1, 'xconstraintfcn',@doConstraints, 'verbose',0, 'plotfcn',[]);
    
    %get the next evaluation point from the model
    nextObj = optOut.predictObjective(optOut.NextPoint);
    nextObj = max(nextObj, minObs);
    
    nextTime = optOut.predictObjectiveEvaluationTime(optOut.NextPoint);
    a = makeparamsfile(table2array(optOut.NextPoint),allFileInds,outPath);
    
    %append this point to the model
    allFileInds = cat(1,allFileInds,a);
    s = cat(1,s,optOut.NextPoint);
    loss = cat(1,loss,nextObj);
    trainTime = cat(1,trainTime,nextTime);
end

end

function tf = doConstraints(x)
tf = ((x.probNeg + x.probPos) <= 1) & (x.treeStopCount <= x.minTrees) & (x.minTrees <= x.maxTrees);
end
