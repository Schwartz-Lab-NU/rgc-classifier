function P = readparamsfile(filename)

P=struct();
f = fopen(filename,"r");

%nFeat, nFolds, nLambda, alpha
P.nFeatures = fscanf(f,"%d\n",1);
P.nFolds = fscanf(f,"%d\n",1);
P.nLambda = fscanf(f,"%d\n",1);
P.alpha = fscanf(f,"%f\n",1);

%nRepeats maxDepth minSize
P.nRepeats = fscanf(f,"%d\n",1);
P.maxDepth = fscanf(f,"%d\n",1);
P.minSize = fscanf(f,"%d\n",1);

%maxTrees, minTrees, stopCount, stopThresh
P.maxTrees = fscanf(f,"%d\n",1);
P.minTrees = fscanf(f,"%d\n",1);
P.stopCount = fscanf(f,"%d\n",1);
P.stopThreshold = fscanf(f,"%f\n",1);

%ensembleSize, probNeg, probPos, nLabels
P.ensembleSize = fscanf(f,"%d\n",1);
P.negativeProbability = fscanf(f,"%f\n",1);
P.positiveProbability = fscanf(f,"%f\n",1);
P.nLabels = fscanf(f,"%d",1);

%error-correcting output code
P.ECOC = fscanf(f,"%d",[P.nLabels, P.ensembleSize]);
fclose(f);

end