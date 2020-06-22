% The mapping from SO(n + 1) to SE(n)
% 
% Author: Jin Wu
% e-mail: jin_wu_uestc@hotmail.com
% website: www.jinwu.science
%          www.ram-lab.com
%
% References: 
%
%  [1] Wu, J., Liu, M. (2019) 
%               Simultaneous SO(n) Solutions to Hand-eye Calibration Problems,
%               IEEE Transactions on Automation Science and Engineering
%
%  [2] Wu, J., Sun, Y., Wang, M., Liu, M. (2019) 
%               Hand-eye Calibration: 4D Procrustes Analysis Approach,
%               IEEE Transactions on Instrumentation and Measurement




function [R, t, X] = SOnp1_SEn(A, d)
    s = size(A);
    dim = s(1);
    
    t = A(1 : dim - 1, dim) * d;
    R = A(1 : dim - 1, 1 : dim - 1);
    X = [R, t;
         zeros(1, dim - 1), 1];

end