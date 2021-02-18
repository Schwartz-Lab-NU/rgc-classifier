classdef PsthLeaf
    properties
        positiveCount
        negativeCount
        score
        totalWeight
        size = 40
    end
    
    properties (Access=private)
        tree
    end
    
    methods
        function obj = PsthLeaf(f, tree)
            obj.tree = tree;
            obj.positiveCount = fread(f,1,'uint');
            obj.negativeCount = fread(f,1,'uint');
            obj.score = fread(f,1,'double');
            obj.totalWeight = fread(f,1,'double');
        end
        
        function classify(self, indices)
            self.tree.lastScores(indices) = self.score;
        end
        
        function s = getSize(self)
            s = self.size;
        end
    end
    
end