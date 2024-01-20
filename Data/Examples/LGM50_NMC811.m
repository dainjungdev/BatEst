function OCP = LGM50_NMC811
% An example positive electrode OCP from equation (8) of this paper:
% C.-H. Chen et al., Journal of The Electrochemical Society, 167:080534,
% 2020. doi.org/10.1149/1945-7111/ab9050.

% Positive electrode potential
OCP = @(Sp) ...
    - 0.8090*Sp + 4.4875 ...
    - 0.0428*tanh(18.5138*(Sp-0.5542)) ...
    - 17.7326*tanh(15.7890*(Sp-0.3117)) ...
    + 17.5842*tanh(15.9308*(Sp-0.3120));

end