%% Segmentation with non-overlappping surfaces
mex('mex/mex_surfcut_planesep_qpbo.cpp', 'src/*.cpp','-Iinclude');

%% Segmentation with unconstrained surfaces
mex('mex/mex_surfcut.cpp', 'src/*.cpp','-Iinclude');