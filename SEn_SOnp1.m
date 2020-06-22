% The mapping from SE(n) to SO(n + 1)
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


function SA = SEn_SOnp1(A, dd)
s = size(A);
dim = s(1);
RA = A(1 : dim - 1, 1 : dim - 1);
tA = A(1 : dim - 1, dim);
SA = [RA, tA / dd;
      - tA' * RA / dd, 1];
end