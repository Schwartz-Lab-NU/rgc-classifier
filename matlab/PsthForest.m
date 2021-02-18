classdef PsthForest < handle
    properties
        trees
        transform
    end
    
    properties (Access = private)
        ensemble
    end
    
    properties (Dependent)
        size
    end
    
    methods
        function obj = PsthForest(forestDir, ensemble)
            obj.ensemble = ensemble;
            nTrees = numel(dir(forestDir)) - 5;
            obj.trees = PsthTree.empty(nTrees,0);
            obj.transform = struct('score',cell(3,1),'probability',cell(3,1));
            for tree = 1:nTrees
                obj.trees(tree) = PsthTree(sprintf('%s\\tree%d.out',forestDir,tree), ensemble);
            end
            for fold = 1:3
                f = fopen(sprintf('%s\\fold%d.LUT',forestDir,fold));
                fread(f,1,'ushort');
                transform = reshape(fread(f,Inf,'double'),[],2);
                [~,uniqueInd,~] = unique(transform(:,1));
                obj.transform(fold).score = transform(uniqueInd,1);
                obj.transform(fold).probability = transform(uniqueInd,2);
                fclose(f);
            end
        end
        
        function p = classify(self, fold)
            nTrees = numel(self.trees);
            s = zeros(self.ensemble.lastCellCount, 1);
            for t = 1:nTrees
                self.trees(t) = self.trees(t).classify();
                s = s + self.trees(t).lastScores;
            end
            p = interp1(self.transform(fold).score, self.transform(fold).probability, s, 'nearest', 'extrap');
        end
        
        function s = get.size(self)
            s = 0;
            for fold = 1:3
                s = s + 16*numel(self.transform(fold).score);
            end
            for tree = self.trees
                s = s + tree.size;
            end
        end
    end
end