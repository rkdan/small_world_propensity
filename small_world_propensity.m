function [SWP,delta_C,delta_L] = small_world_propensity(A, varargin)

% a function for calculating the small world propensity of
% a given network - assumes that matrix is undirected (symmeteric) and if
% not, creates a symmetric matrix which is used for the calculations

%NOTE:  This code requires the Bioinformatics Toolbox to be installed
%        (uses graphallshortestpaths.m)

%Inputs:
%   A           the connectivity matrix, weighted or binary
%   varargin    a string corresponding to the method of clustering
%               to be used, where 'O' is Onnela, 'Z' is Zhang,
%               'B' is Barrat, 'bin' is binary (default is Onnela).
%               If the user specifies binary analysis, a
%               weighted matrix will be converted to a binary matrix
%               before proceeding.
%

%Outputs:
%   SWP         the small world propensity of the matrix
%   delta_C     the fractional deviation from the expected culstering coefficient of a
%                   random network
%   delta_L     the fractional deviation from the expected path length of a
%                   random network

%written by Eric Bridgeford and modified by Sarah F. Muldoon

% Reference: Muldoon, Bridgeford, and Bassett (2015) "Small-World Propensity in Weighted,
%               Real-World Networks" http://arxiv.org/abs/1505.02194

if isempty(varargin)
    varargin{1} = 'O';
end

if sum(sum(A)) > 0

bin_matrix = 0;
if strcmp(varargin{1},'bin') == 1
   bin_matrix = 1;
   A = A > 0;
end

%check to see if matrix is symmeteric
symcheck=abs(A-A');
if sum(sum(symcheck)) > 0
    % adjust the input matrix to symmeterize
    disp('Input matrix is not symmetric. Symmetrizing.')
    W = symm_matrix(A, bin_matrix);
else
    W=A;
end

%calculate the number of nodes
n = length(W);
%compute the weighted density of the network
dens_net = sum(sum(W))/(max(max(W))*n*(n-1));

%compute the average degree of the unweighted network, to give
%the approximate radius
numb_connections = length(find(W>0));
avg_deg_unw = numb_connections/n;
avg_rad_unw = avg_deg_unw/2;
avg_rad_eff = ceil(avg_rad_unw);


%compute the regular and random matrix for the network W
W_reg = regular_matrix_generator(W, avg_rad_eff);
W_rand = randomize_matrix(W);
%compute all path length calculations for the network
reg_path = avg_path_matrix(1./W_reg);      %path of the regular network
rand_path = avg_path_matrix(1./W_rand);    %path of the random netowork
net_path = avg_path_matrix(1./W);          %path of the network

A = (net_path - rand_path);
if A < 0
    A = 0;
end
diff_path =  A/ (reg_path - rand_path);
if net_path == Inf || rand_path == Inf || reg_path == Inf
    diff_path = 1;
end
if diff_path > 1
    diff_path = 1;
end


%compute all clustering calculations for the network
reg_clus = avg_clus_matrix(W_reg,varargin{1});
rand_clus = avg_clus_matrix(W_rand,varargin{1});
net_clus = avg_clus_matrix(W,varargin{1});

B = (reg_clus - net_clus);
if B < 0
    B = 0;
end

diff_clus = B / (reg_clus - rand_clus);
if isnan(reg_clus) || isnan(rand_clus) || isnan(net_clus)
    diff_clus = 1;
end
if diff_clus > 1
    diff_clus = 1;
end

%calculate small world value, the root sum of the squares of
%diff path and diff clus
SWP = 1 - (sqrt(diff_clus^2 + diff_path^2)/sqrt(2));
delta_C=diff_clus;
delta_L=diff_path;
end
end


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%the above code requires the following funcitons

function [Clus] = avg_clus_matrix(W, met)
%a function to compute the average clusteirng coefficient for a
%input matrix M

%Inputs:
%   W     a matrix, weighted or unweighted
%   met   a string, to represent the method to be used for computing
%         the clustering coefficient
%         possible strings: 'O' (Onnela), 'Z' (Zhang), 'B' (Barrat),
%         'bin' (binary)
%         default if none is chosen is Onnela

%Outputs:
%   Clus  the average clustering coefficient

%written by Eric Bridgeford

n = length(W);
[C] = clustering_coef_matrix(W, met);
Clus = nanmean(C);
end


function [Len] = avg_path_matrix(M)

%a function to compute the average path length of a given matrix
%using the graphallshortestpaths built-in matlab function

%written by Eric Bridgeford

n = length(M);
M = sparse(M);
D = graphallshortestpaths(M);

%checks if a node is disconnected from the system, and replaces
%its value with 0
for i = 1:n
    for j = 1:n
        if isinf(D(i,j)) == 1
            D(i,j) = 0;
        end
    end
end

Len = mean(mean(D));
end


function [C] = clustering_coef_matrix(W, met)

%a modification of the clustering coefficient function provided
%in the brain connectivity toolbox

%improved definition of Onnela Clustering Coefficient, as well as
%implementation of function for Zhang and Barrat clustering values

%Reference:
%   Onnela et al., Phys. Rev. E71, 065103(R)(2005)
%   B.Zhang and S. Horvath, Stat. App. Genet. Mol. Biol.4, 17(2005)
%   Barrat et al., Proc. Natl. Acad. Sci. U.S.A.101, 3747(2004)
%   Watts and Strogatz (1998) Nature 393:440-442

%Inputs:
%   W    the weighted or unweighted connectivity matrix
%   met   a string, to represent the method to be used for computing
%         the clustering coefficient
%         possible strings: 'O' (Onnela), 'Z' (Zhang), 'B' (Barrat), 'bin'
%         (binary)
%         default if none is chosen is Onnela


%code originally written by Mika Rubinov, UNSW, 2007-2010
%modified/written by Eric Bridgeford

if met == 'O'
    K=sum(W~=0,2);
    W = double(W);
    W2 = W/max(max(W));
    cyc3=diag(W2.^(1/3)^3);
    K(cyc3==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
    C=cyc3./(K.*(K-1));
end
if met == 'bin'
    G = double(W>0);
    n=length(G);
    C=zeros(n,1);
    for u=1:n
        V=find(G(u,:));
        k=length(V);
        if k>=2;                %degree must be at least 2
            S=G(V,V);
            C(u)=sum(S(:))/(k^2-k);
        end
    end
end

if met == 'Z'
    K=sum(W~=0,2);
    W = double(W);
    W2 = W/max(max((W)));
    cyc3=diag((W2)^3);
    denom = zeros(length(W),1);
    for i = 1:length(W)
        denom(i) = (sum(W2(i,:))^2-sum(W2(i,:).^2));
    end
    C = cyc3./denom;
end

if met == 'B'
    A = double(W>0);
    C = zeros(length(W),1);
    for i = 1:length(W)
        sum1 = 0;
        for j = 1:length(W)
            for k = 1:length(W)
                sum1 = ((W(i,j)+W(i,k))/2)*A(i,j)*A(j,k)*A(i,k)+sum1;
            end
        end
        C(i) = 1/(sum(W(i,:))*(sum(A(i,:))-1))*sum1;
    end
end


end


function A_rand=randomize_matrix(A);

%This code creates a random undirected network from the connectivity
%distribution of an undirected adjacency matrix, ie, the intital matrix
%must be symmetric.

% INPUTS:
%   A: an undirected adjacency matrix (symmetric) with no self connections

% OUTPUTS:
%   A_rand: a comparable random network with same number of nodes and
%       connectivity distribution

% written by Sarah F. Muldoon

num_nodes=length(A);
A_rand=zeros(num_nodes);
mask=triu(ones(num_nodes),1);
grab_indices=find(mask > 0);

orig_edges=A(grab_indices);
num_edges=length(orig_edges);

rand_index=randperm(num_edges);
randomized_edges=orig_edges(rand_index);

edge=1;
for i=1:num_nodes-1
    for j=i+1:num_nodes
        A_rand(i,j)=randomized_edges(edge);
        A_rand(j,i)=randomized_edges(edge);
        edge=edge+1;
    end
end
end


function M = regular_matrix_generator(G,r)
%generates a regular matrix, with weights obtained form the
%original adjacency matrix representation of the network

% note that all inputs should be symmeterized prior to forming a regular
% matrix, since otherwise half of the connnections will be trashed. This
% can be accomplished with the built in symm_matrix function, however,
% the function is not performed here so that users can use their own
% symmeterization procedure.

%Inputs:
%   G    the adjacency matrix for the given network; must be symmmeterized
%   r    the approximate radius of the regular network

%Outputs:
%   M    the regular matrix for the given network, where all
%        weights are sorted such that the inner radius has the
%        highest weights randomly distributed across the nodes,
%        and so on

%written by Eric W. Bridgeford

n = length(G);
G = triu(G);
%reshape the matrix G into an array, B
B = reshape(G,[length(G)^2,1]);
%sorts the array in descending order
B = sort(B,'descend');
%computes the number of connections and adds zeros if
%numel(G) < 2*n*r
num_els =ceil(numel(G)/(2*n));
num_zeros = 2*n*num_els - numel(G);
%adds zeros to the remaineder of the list, so length(B) = 2*n*r
B = cat(1,B,zeros(num_zeros,1));
%reshapes B into a matrix, where the values descend top to
%bottom, as well as left to right. The greatest value in each
%column is less than the smallest value of the column to its left.
B = reshape(B,[n],[]);

M = zeros(length(G));

%distributes the connections into a regular network, M, where
%the innermost radius represents the highest values
for i = 1:length(G)
    for z = 1:r
        a = randi([1,n]);

        %random integer chosen to take a value from B
        while (B(a,z) == 0 && z ~= r) || (B(a,z) == 0 && z == r && ~isempty(find(B(:,r),1)))
            a = randi([1,n]);
        end
        %finds the two nodes a distance of z from the origin node
        %and places the entries from the matrix B
        y_coor_1 = mod(i+z-1,length(G))+1;
        %y_coor_2 = mod(i-z-1,length(G))+1;
        M(i,y_coor_1) = B(a,z);
        M(y_coor_1,i) = B(a,z);
        %removes the weights from the matrix B so they cannot be
        %reused
        B(a,z) = 0;

    end
end
end


function [W] = symm_matrix(A, bin_key)

% a function to symmetrize an input matrix. The procedure
% by which items are symmetrized such that:
%   in the binary case:
%       if a(i,j) || a(j,i) == 1 for i,j in A
%           w(i,j) && w(j,i) == 1
%   in  the weighted case:
%       if (a(i,j) || a(j,i) > 0 for i,j in A
%           w(i,j) && w(j,i) == (a(i,j) + a(j,i) )/ 2

% Inputs:
%   A:          The binary or weighted input matrix
%   bin_key:    the key to indicate whether weighted or binary analysis
%               will take place
%               1 indicates binarized, 0 indicates weighted

% Outputs
%   W:          The symmeterized matrix

% if binary analysis is specified, let binary symmeterization take place

% written by Eric W. Bridgeford

W = zeros(length(A));

if bin_key == 1
    A = A > 0; % verify that the input matrix is binary
    for i = 1:length(A)
        for j = i:length(A)
            if A(i,j) || A(j,i)
                W(i,j) = 1;
                W(j,i) = 1;
            end
        end
    end
else
    for i = 1:length(A)
        for j = i:length(A)
            if A(i,j) || A(j,i)
                val = (A(i,j) + A(j,i)) / 2;
                W(i,j) = val;
                W(j,i) = val;
            end
        end
    end
end
end
