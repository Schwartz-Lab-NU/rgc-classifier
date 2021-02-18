classdef PsthTree < PsthNode & handle
    properties
        sampleWeights
        treeWeight
        lastScores
    end
    
    properties (Access = private)
        ensemble
    end
    
    properties (Dependent)
        size
    end
    
    methods
        function obj = PsthTree(treeName, ensemble)
            f = fopen(treeName);
            N = fread(f,1,'int');
            fread(f,88,'char'); %paramset info
            fread(f,2*5016,'char'); %seed info

            sampleWeights = fread(f,N,'double'); %corresponds to training fold
            treeWeight = fread(f,1,'double');
            fread(f,1,'char');
            
            obj = obj@PsthNode(f, [], ensemble);
            obj.ensemble = ensemble;
            obj.sampleWeights = sampleWeights;
            obj.treeWeight = treeWeight;
            fclose(f);
        end
        
        function self = classify(self)
            self.lastScores = zeros(self.ensemble.lastCellCount,1);
            classify@PsthNode(self, 1:self.ensemble.lastCellCount);
            self.lastScores = self.treeWeight * self.lastScores;
        end
        
        function s = get.size(self)
            nw = numel(self.timeIndex);
            
            %matlab doesn't allow overriding base class getters
            s = 64 + nw*64 + 8*numel(self.sampleWeights) + 8*numel(self.lastScores) + self.left.getSize() + self.right.getSize();
        end
    end
    
    
end