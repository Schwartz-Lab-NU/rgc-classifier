classdef PsthNode < handle
    properties
        N
        timeIndex
        timeIndexQuality
        spotSizeIndex
        spotSizeIndexQuality
        indexQuality
        beta0
        beta1
        mu
        istd
        lambda
        deviance
        uncertainty
        support
        left
        right
    end
    
    properties (Dependent, Access = private)
        size
    end
    
    properties (Access = private)
        tree
        ensemble
    end
    
    methods
        function obj = PsthNode(f, tree, ensemble)
            obj.N = fread(f,1,'uint');
            nw = fread(f,1,'uint');
            obj.timeIndex = fread(f,nw,'double') + 1;
            obj.spotSizeIndex = fread(f,nw,'double') + 1;
            obj.beta1 = fread(f,nw,'double');
            obj.mu = fread(f,nw,'double');
            obj.istd = fread(f,nw,'double');
            obj.timeIndexQuality = fread(f,nw,'double');
            obj.spotSizeIndexQuality = fread(f,nw,'double');
            obj.indexQuality = fread(f,nw,'double');
            
            obj.beta0 = fread(f,1,'double');
            obj.lambda = fread(f,1,'double');
            obj.deviance = fread(f,1,'double');
            obj.uncertainty = fread(f,1,'double');
            obj.support = fread(f,1,'double');
            
            if isempty(tree)
               obj.tree = obj;
            else
                obj.tree = tree;
            end
            obj.ensemble = ensemble;
            %save node
            left = fread(f,1,'char');
            if left==1
                obj.left = PsthNode(f, obj.tree, ensemble);
            else
                obj.left = PsthLeaf(f, obj.tree);
            end
            right = fread(f,1,'char');
            if right==1
                obj.right = PsthNode(f, obj.tree, ensemble);
            else
                obj.right = PsthLeaf(f, obj.tree);
            end
        end
        
        function classify(self, indices)
            % psth is (nTimes -by- nSpots -by- nCells)
            nCells = numel(indices);
            li = sub2ind([351,1601],self.timeIndex,self.spotSizeIndex);
            li = sub2ind([351*1601, self.ensemble.lastCellCount], repmat(li,1,nCells), repmat(indices,numel(li), 1));
            
            decision = sum(self.beta1.*self.istd.*(self.ensemble.lastPsth(li) - self.mu),1) <= self.beta0;
            if nnz(decision)
                self.left.classify(indices(decision));
            end
            
            decision = ~decision;
            
            if nnz(decision)
                self.right.classify(indices(decision));
            end
            
        end
        
        function s = get.size(self)
            nw = numel(self.timeIndex);
            s = 56 + nw*64 + self.left.getSize() + self.right.getSize();
        end
        
        function s = getSize(self)
            s = self.size;
        end
    end    
end