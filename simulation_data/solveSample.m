% Worker-side helper for parallel dataset generation.
% force_pos : (N_FORCES x 1)  normalised arc-length per force
% force_vec : (3 x N_FORCES)  [Fx;Fy;Fz] per force column, body frame

function [base_T, tip_T, meas_strain, gt_transforms, force_vec_global, disc_positions] = ...
    solveSample(manipulator, Xs, node_sig_idx, strain_nodes, ...
               n_strain_meas, N_GT, nNodes, force_pos, force_vec, ...
               N_FORCES, tendon_tension, disc_i_sig)

global ext_forces u_cust %#ok<GVMIS>

u_cust = tendon_tension;

ext_forces    = cell(N_FORCES, 2);
force_sig_idx = zeros(N_FORCES, 1);
for f = 1:N_FORCES
    [~, idx]         = min(abs(force_pos(f) - Xs));
    force_sig_idx(f) = idx;
    ext_forces{f, 1} = idx;
    ext_forces{f, 2} = force_vec(:, f);
end

q = manipulator.statics(zeros(manipulator.ndof, 1));

g_full  = manipulator.FwdKinematics(q);
xi_full = ScrewStrain(manipulator, q);

i_sig_base = node_sig_idx(1);
base_T     = g_full((i_sig_base-1)*4+1:i_sig_base*4, :);

i_sig_tip = node_sig_idx(nNodes);
tip_T     = g_full((i_sig_tip-1)*4+1:i_sig_tip*4, :);

meas_strain = zeros(3, n_strain_meas);
for k = 1:n_strain_meas
    i_sig            = node_sig_idx(strain_nodes(k) + 1);
    xi_here          = xi_full((i_sig-1)*6+1:i_sig*6);
    meas_strain(:,k) = xi_here(1:3);
end

gt_transforms = zeros(4, 4*N_GT);
for i_gt = 1:N_GT
    s      = (i_gt - 1) / (N_GT - 1);
    [~, i] = min(abs(s - Xs));
    gt_transforms(:, (i_gt-1)*4+1:i_gt*4) = g_full((i-1)*4+1:i*4, :);
end

force_vec_global = zeros(3, N_FORCES);
for f = 1:N_FORCES
    idx    = force_sig_idx(f);
    g_here = g_full((idx-1)*4+1:idx*4, :);
    R_here = g_here(1:3, 1:3);
    force_vec_global(:, f) = R_here * force_vec(:, f);
end

disc_positions = zeros(11, 3);
for k = 1:numel(disc_i_sig)
    gh = g_full((disc_i_sig(k)-1)*4+1:disc_i_sig(k)*4, :);
    disc_positions(k, :) = gh(1:3, 4)';
end

end
