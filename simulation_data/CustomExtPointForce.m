%Function to calculate the custom external point force. Make sure to edit
%CustomExtForce.m file too
%Last modified by Aditya Prakash 5/20/2026

function Fext=CustomExtPointForce(Tr,q,g,J,t,qd,eta,Jdot)

global ext_forces
n    = Tr.nsig-Tr.N; %everything except rigid joints
Fext = zeros(6*n,1);


for i = 1:size(ext_forces, 1)
    i_sig = ext_forces{i, 1};
    g_here = g((i_sig-1)*4+1:i_sig*4,:);
    Ad_g_here_inv = dinamico_Adjoint(ginv(g_here));
    F_vec = ext_forces{i, 2};
    i_sig_nj = i_sig-1;
    Fext((i_sig_nj-1)*6+1:i_sig_nj*6) = [0; 0; 0; F_vec]; % this has been verified. It is correct.
end
end