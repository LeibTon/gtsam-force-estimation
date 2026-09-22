% Generates static-equilibrium samples for the experimental-geometry rod
% (manipulator_new.mat) and writes them to HDF5, in the layout GTSAM's
% estimator reads directly (see GTSAM_SOROSIM_NODE_MAPPING.md).

clc; clear;

N_FORCES   = 2;
N_SAMPLES  = 4000;
F_MIN      = 0.1;      % N, body-frame magnitude; For single force: 0.15 N; for two forces: 0.1 N
F_MAX      = 0.25;      % N, body-frame magnitude; For single force: 0.3 N; For two forces: 0.25 N
S_MIN      = 0.37;      % normalised arc-length
S_MAX      = 0.95;
N_GT       = 101;
nPieces    = 20;
outFile    = sprintf('dataset_%dforces.h5', N_FORCES);

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'Force Estimation Paper', 'SoRoSim V6.3')));
load('manipulator_new.mat');
manipulator_act = manipulator_new;

manipulator_act.CEFP = true;
manipulator_act.CAS  = false;   % no tendon actuation on this rod

nNodes  = nPieces + 1;
g0      = manipulator_act.FwdKinematics(zeros(manipulator_act.ndof, 1));
Xs_phys = zeros(manipulator_act.nsig, 1);
for i_sig = 1:manipulator_act.nsig
    g_here         = g0((i_sig-1)*4+1:i_sig*4, :);
    Xs_phys(i_sig) = g_here(1, 4);
end
Xs = Xs_phys / Xs_phys(end);   % normalised [0, 1]

node_sig_idx = zeros(1, nNodes);
for n = 0:nPieces
    [~, i]            = min(abs(n/nPieces - Xs));
    node_sig_idx(n+1) = i;
end

disc_i_sig = 1:9:91;   % 11 exact disc positions, stride 9 (GTSAM_SOROSIM_NODE_MAPPING.md)

strain_nodes  = 1:2:nPieces-1;
n_strain_meas = length(strain_nodes);
tendon_tension = zeros(0, 1);   % no tendons on this rod

strain_arc_lengths = double(strain_nodes) / nPieces;
shape_arc_lengths  = linspace(0, 1, N_GT);
disc_arc_lengths   = (0:10) / 10;

rng(42);
all_force_pos = S_MIN + (S_MAX - S_MIN) * rand(N_FORCES, N_SAMPLES);
all_dirs      = randn(2, N_FORCES * N_SAMPLES);   % Fy, Fz directions only
all_dirs      = all_dirs ./ vecnorm(all_dirs, 2, 1);
all_mags      = F_MIN + (F_MAX - F_MIN) * rand(N_FORCES, N_SAMPLES);
all_yz        = reshape(all_dirs .* all_mags(:)', 2, N_FORCES, N_SAMPLES);
all_force_vec_body = [zeros(1, N_FORCES, N_SAMPLES); all_yz];   % Fx=0, body frame

res_base_T    = zeros(4, 4, N_SAMPLES);
res_tip_T     = zeros(4, 4, N_SAMPLES);
res_strain    = zeros(3, n_strain_meas, N_SAMPLES);
res_gt_T      = zeros(4, 4*N_GT, N_SAMPLES);
res_force_pos = zeros(N_FORCES, N_SAMPLES);
res_force_vec = zeros(3, N_FORCES, N_SAMPLES);   % global frame
res_disc_pos  = zeros(11, 3, N_SAMPLES);

pool = gcp('nocreate');
if isempty(pool), pool = parpool('local'); end
fprintf('parallel workers: %d\n', pool.NumWorkers);
fprintf('%d samples, %d force(s), F_body in [%.3f %.3f] N, s in [%.2f %.2f]\n', ...
        N_SAMPLES, N_FORCES, F_MIN, F_MAX, S_MIN, S_MAX);

tic;
parfor s_idx = 1:N_SAMPLES
    fp = all_force_pos(:, s_idx);
    fv_body = reshape(all_force_vec_body(:, :, s_idx), 3, N_FORCES);

    [base_T, tip_T, meas_strain, gt_transforms, fv_global, disc_pos] = ...
        solveSample(manipulator_act, Xs, node_sig_idx, strain_nodes, ...
                   n_strain_meas, N_GT, nNodes, fp, fv_body, ...
                   N_FORCES, tendon_tension, disc_i_sig);

    res_base_T(:,:,s_idx)    = base_T;
    res_tip_T(:,:,s_idx)     = tip_T;
    res_strain(:,:,s_idx)    = meas_strain;
    res_gt_T(:,:,s_idx)      = gt_transforms;
    res_force_pos(:,s_idx)   = fp;
    res_force_vec(:,:,s_idx) = fv_global;
    res_disc_pos(:,:,s_idx)  = disc_pos;
end
elapsed = toc;
fprintf('solved %d samples in %.1f s (%.2f s/sample wall)\n', N_SAMPLES, elapsed, elapsed/N_SAMPLES);

if isfile(outFile), delete(outFile); end

h5create(outFile, '/measurements/base_pose/transform',  [4 4 N_SAMPLES]);
h5write( outFile, '/measurements/base_pose/transform',  res_base_T);

h5create(outFile, '/measurements/tip_pose/transform',   [4 4 N_SAMPLES]);
h5write( outFile, '/measurements/tip_pose/transform',   res_tip_T);

h5create(outFile, '/measurements/strain/values',        [3 n_strain_meas N_SAMPLES]);
h5write( outFile, '/measurements/strain/values',        res_strain);

h5create(outFile, '/measurements/strain/arc_lengths',   [1 n_strain_meas]);
h5write( outFile, '/measurements/strain/arc_lengths',   strain_arc_lengths);

n_tendons = numel(tendon_tension);
if n_tendons > 0
    h5create(outFile, '/measurements/tendon_tensions',  [n_tendons N_SAMPLES]);
    h5write( outFile, '/measurements/tendon_tensions',  repmat(tendon_tension, 1, N_SAMPLES));
end

h5create(outFile, '/ground_truth/shape/transforms',     [4 4*N_GT N_SAMPLES]);
h5write( outFile, '/ground_truth/shape/transforms',     res_gt_T);

h5create(outFile, '/ground_truth/shape/arc_lengths',    [1 N_GT]);
h5write( outFile, '/ground_truth/shape/arc_lengths',    shape_arc_lengths);

h5create(outFile, '/ground_truth/disc/positions',       [11 3 N_SAMPLES]);
h5write( outFile, '/ground_truth/disc/positions',       res_disc_pos);

h5create(outFile, '/ground_truth/disc/arc_lengths',     [1 11]);
h5write( outFile, '/ground_truth/disc/arc_lengths',     disc_arc_lengths);

h5create(outFile, '/ground_truth/force/position',       [N_FORCES N_SAMPLES]);
h5write( outFile, '/ground_truth/force/position',       res_force_pos);

h5create(outFile, '/ground_truth/force/vec',            [3 N_FORCES N_SAMPLES]);
h5write( outFile, '/ground_truth/force/vec',            res_force_vec);

h5create(outFile, '/metadata/n_samples',  [1 1]); h5write(outFile, '/metadata/n_samples',  double(N_SAMPLES));
h5create(outFile, '/metadata/n_forces',   [1 1]); h5write(outFile, '/metadata/n_forces',   double(N_FORCES));

fprintf('Done. Saved %d samples (%d forces each) to %s\n', N_SAMPLES, N_FORCES, outFile);
