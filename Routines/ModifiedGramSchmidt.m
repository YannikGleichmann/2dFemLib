function [Q,r] = ModifiedGramSchmidt(Q,B,l)

% This method implements the modified Gram-Schmidt orthonormalization with
% respect to an arbitrary bilinear form B(v, w) = v'*B*w with a matrix B.
% This algorithm overwrites the input matrix Q to the new orthogonal Matrix
% Q.
% 
% Reference:
% [1] Van Loan, Charles F., and Gene H. Golub. | Matrix computations. p.232
% 
% Usage:
%           Q = ModifiedGramSchmidt(Q);
%           Q = ModifiedGramSchmidt(Q, B);
%           [Q, r] = ModifiedGramSchmidt(__);
% 
% Input:
% * Q = [q_1, ..., q_n] - n x m (n > m) matrix with ordered vectors q_i 
% * B                   - the matrix for the bilinear form B(v, w)
%                           default: I (euclidean norm)
% 
% Output:
% * Q = [q_1, ..., q_n] - the (new) orthonormalized vectors
% * r                   - diagonal of R in the QR decomposition
%                         a small value in r(i) indicates that the i-th
%                         column of Q may be linearly dependent 
% 
% 
% Author: Yannik Gleichmann
% Date:   29.01.2024

[m,n] = size(Q);
if exist('B','var') == 0 || isempty(B); B = speye(m,m); end
if exist('l','var') == 0 || isempty(l); l = 1; end

r = zeros(n,1);
for i = l:n
  q_i = Q(:,i);
  r(i) = sqrt(q_i'*(B*q_i));
  if r(i) < eps*1e1
    Q(:,i) = zeros(size(q_i));
  else
    Q(:,i) = q_i/r(i);
    q_i = Q(:,i);
    q_j = Q(:,i+1:n);
    Q(:,i+1:n) = q_j - q_i*(q_i'*(B*q_j));
  end
end
