function N=makeparamsfile(trace,k,path) 

% rng shuffle; %seeds from current time
if isempty(k)
    k=10000;
end

N = k(1);
while any(k==N)
    N = randi(9999);
end
f = fopen(sprintf("%s%04dparams.in",path,N),"w");

%nFeat, nFolds, nLambda, alpha
fprintf(f,"%d\n",trace(1));
fprintf(f,"%d\n",trace(2));
fprintf(f,"%d\n",trace(3));
fprintf(f,"%f\n",trace(4));

%nRepeats maxDepth minSize
fprintf(f,"%d\n",trace(5));
fprintf(f,"%d\n",trace(6));
fprintf(f,"%d\n",trace(7));

%maxTrees, minTrees, stopCount, stopThresh
fprintf(f,"%d\n",trace(8));
fprintf(f,"%d\n",trace(9));
fprintf(f,"%d\n",trace(10));
fprintf(f,"%f\n",trace(11));

%ensembleSize, probNeg, probPos, nLabels
fprintf(f,"%d\n",trace(12));
fprintf(f,"%f\n",trace(13));
fprintf(f,"%f\n",trace(14));
fprintf(f,"%d",31);

%ecoc has been moved to cpp
% M = designecocN(31, trace(12), [trace(13) 1-trace(13)-trace(14) trace(14)]);
% 
% for m=1:trace(12)
% for n=1:31
% fprintf(f, '%d', M(n,m));
% if n<31
%     fprintf(f," ");
% elseif m<trace(12)
% 	fprintf(f,"\n");
% end
% end
% end
fclose(f);
