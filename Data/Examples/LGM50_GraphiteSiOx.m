function OCP = LGM50_GraphiteSiOx
% An example negative electrode OCP from equation (9) of this paper:
% C.-H. Chen et al., Journal of The Electrochemical Society, 167:080534,
% 2020. doi.org/10.1149/1945-7111/ab9050.

% Negative electrode potential
OCP = @(Sn) ...
        1.9793*exp(-39.3631*Sn) + 0.2482 ...
      - 0.0909*tanh(29.8538*(Sn-0.1234)) ...
      - 0.04478*tanh(14.9159*(Sn-0.2769)) ...
      - 0.0205*tanh(30.4444*(Sn-0.6103));

end
