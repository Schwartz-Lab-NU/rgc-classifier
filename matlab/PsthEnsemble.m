classdef PsthEnsemble < handle
    properties
        forests
        ecoc %k-by-l
        lastPsth
        lastCellCount
    end
    
    properties (Dependent)
        size
    end
    
    methods 
        function obj = PsthEnsemble(ensembleDir, ecoc)
            obj.ecoc = ecoc;
            nForests = numel(dir(ensembleDir)) - 3;
            obj.forests = PsthForest.empty(nForests,0);
            for forest= 1:nForests
                forestDir = sprintf('%s\\forest%d', ensembleDir, forest);
                obj.forests(forest) = PsthForest(forestDir, obj);
            end
        end
        
        function p = classify(self, psth, normalizationFold)
            self.lastPsth = psth;
            self.lastCellCount = size(psth,3); %#ok<CPROPLC>
            
            nForests = numel(self.forests);
            s = zeros(self.lastCellCount, nForests);
            for f = 1:nForests
                s(:,f) = self.forests(f).classify(normalizationFold);
            end
            %NOTE: n is the number of filled rows in the ecoc matrix, because
            %all the samples were weighted by class
            p = ECOC_coupling(self.ecoc,s,sum(abs(self.ecoc),1)');
        end
    
        function s = get.size(self)
            %returns the size of the ensemble in bytes
            s = 8*(numel(self.ecoc) + numel(self.lastPsth) + 1);
            for forest = self.forests
                s = s + forest.size;
            end
        end
    end
end