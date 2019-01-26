%save('/Users/Jacopo/dfg4grn/Arabidopsis/Jacopo.mat','KNO3_gen')


% GRN_CreateCombinations
%   Create replicate-based combinations of micro-array data
%
% Syntax:
%   [nComb, KNO3, KCl, ratio] = ...
%     GRN_CreateCombinations(geneNames, use_log, dataKNO31, dataKNO32, ...
%                            [dataKCl1, dataKCl2])
%
% Inputs:
%   <geneNames>  : cell of size <nGenes> with gene names
%   <use_log>    : use log normalization or not
%   <dataKNO31>  : matrix of size <nGenes> x <nTimes> containing
%                  micro-array data for experiments on the reaction to KN03
%                  for replicate 1
%   <dataKNO32>  : matrix of size <nGenes> x <nTimes> containing
%                  micro-array data for experiments on the reaction to KN03
%                  for replicate 2
%   <dataKCl1>   : matrix of size <nGenes> x <nTimes> containing
%                  micro-array data for experiments on the reaction to KCl
%                  for replicate 1 (optional)
%   <dataKCl2>   : matrix of size <nGenes> x <nTimes> containing
%                  micro-array data for experiments on the reaction to KCl
%                  for replicate 2 (optional)
%
% Return:
%   <nComb>  : number of combinations of replicates 
%              (for 2 replicates A and B, there are 4 combinations for
%              a time series of 5 time points:
%              AAAAA, ABABA, BBBBB, BABAB
%   <KNO3>   : cell array with <nComb> cells, one for each combination 
%              of replicates, each containing a matrix of size 
%              <nGenes> x <nTimes> with micro-array data for KNO3
%   <KCL>    : cell array with <nComb> cells, one for each combination 
%              of replicates, each containing a matrix of size 
%              <nGenes> x <nTimes> with micro-array data for KCl (optional)
%   <ratio>  : cell array with <nComb> cells, one for each combination 
%              of replicates, each containing a matrix of size 
%              <nGenes> x <nTimes> with micro-array data for KNO3/KCl (optional)
%
% Copyright (C) 2009 Piotr Mirowski
%                    Courant Institute of Mathematical Sciences
%                    New York University

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% Version 1.0, New York, 5 July 2009
% (c) 2009, Piotr Mirowski,
%     Ph.D. candidate at the Courant Institute of Mathematical Sciences
%     Computer Science Department
%     New York University
%     719 Broadway, 12th Floor, New York, NY 10003, USA.
%     email: mirowski [AT] cs [DOT] nyu [DOT] edu

function [nComb, KNO3, KCl] = ...
  GRN_CreateCombinations_shoot(geneNames, use_log, dataKNO31, dataKNO32, dataKNO33, ...
                         dataKCl1, dataKCl2, dataKCl3)

nGenes = length(geneNames);

if (use_log)
  dataKNO3{1} = log2(dataKNO31);
  dataKNO3{2} = log2(dataKNO32);
  dataKNO3{3} = log2(dataKNO33);
  if nargin == 8
    dataKCl{1} = log2(dataKCl1);
    dataKCl{2} = log2(dataKCl2);
    dataKCl{3} = log2(dataKCl3);
  end
else
  dataKNO3{1} = dataKNO31;
  dataKNO3{2} = dataKNO32;
  dataKNO3{3} = dataKNO33;
  if nargin == 8
    dataKCl{1} = dataKCl1;
    dataKCl{2} = dataKCl2;
    dataKCl{3} = dataKCl3;
  end
end

% Create new series by combining replicates
lenComb = size(dataKNO31, 2);
switch lenComb
  case 5
    comb = ['11111'; '12121'; '21212'; '22222'];
  case 7
    comb = ['1111111'; '1212121'; '2121212'; '2222222'];
  case 8
    comb = ['11111111'; '12121212'; '21212121'; '22222222'; '33333333'; '13131313'; '23232323'; '32323232'; '31313131'];
  case 10
    comb = ['1111111111'; '1212121212'; '2121212121'; '2222222222'; '3333333333'; '1313131313'; '2323232323'; '3232323232'; '3131313131'];
  otherwise 
    error('Not implemented');
end
nComb = size(comb, 1);

% Initialize output data
KNO3 = cell(nComb, 1);
if nargin == 8
  KCl = cell(nComb, 1);
  ratio = cell(nComb, 1);
end

% Loop over all combinations
for k = 1:nComb
  KNO3{k} = zeros(nGenes, lenComb);
  if nargin == 8
    KCl{k} = zeros(nGenes, lenComb);
  end
  c = comb(k, :);
  for j = 1:lenComb
    if (c(j) == '1')
      KNO3{k}(:, j) = dataKNO3{1}(:, j);
      if nargin == 8
        KCl{k}(:, j) = dataKCl{1}(:, j);
        ratio{k}(:, j) = dataKNO3{1}(:, j) - dataKCl{1}(:, j);
      end
    elseif (c(j) == '2')
      KNO3{k}(:, j) = dataKNO3{2}(:, j);
      if nargin == 8
        KCl{k}(:, j) = dataKCl{2}(:, j);
        ratio{k}(:, j) = dataKNO3{2}(:, j) - dataKCl{2}(:, j);
      end
    elseif (c(j) == '3')
      KNO3{k}(:, j) = dataKNO3{3}(:, j);
      if nargin == 8
        KCl{k}(:, j) = dataKCl{3}(:, j);
        ratio{k}(:, j) = dataKNO3{3}(:, j) - dataKCl{3}(:, j);
      end
    end
  end
end
