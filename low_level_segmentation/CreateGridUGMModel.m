function [edgePot,edgeStruct]=CreateGridUGMModel(nRows, nCols, K, lambda)
%
%
% nRows, nCols: image dimension
% K: number of states
% lambda: smoothing factor

tic

nNodes = nRows*nCols;

adj = sparse(nNodes,nNodes);

% Add Down Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols); % No Down edge for last row
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1)) = 1;

% Add Right Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows])); % No right edge for last column
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+nRows)) = 1;

% Add Up/Left Edges
adj = adj+adj';
edgeStruct = UGM_makeEdgeStruct(adj,K);

edgePot = zeros(K, K);
for i = 1:K
    for j = 1:K
        if i == j
            edgePot(i,j) = exp(lambda(2));
        else
            edgePot(i,j) = exp(lambda(1));
        end
    end
end
edgePot = repmat(edgePot, [1 1 edgeStruct.nEdges]);

toc;
