%% Make volume containing three overlapping balls
V1 = false(100,100,100);
V1(35,50,35) = true;
V1 = single(bwdist(V1) <= 30);

V2 = false(100,100,100);
V2(70,50,35) = true;
V2 = single(bwdist(V2) < 30);

V3 = false(100,100,100);
V3(50,50,70) = true;
V3 = single(bwdist(V3) < 30);

V = max(V1,V2);
V = max(V,V3);
clear V1 V2 V3;

%% Display volume
figure(1); clf;
patch(isosurface(V,0.5),'EdgeColor','none','FaceColor',[0.8,0.8,0.8]);
axis equal;
axis([1 size(V,2) 1 size(V,1) 1 size(V,3)]);
view(75,15);
lighting phong
camlight right
material dull
title('Data');

%% Segment balls
% Define (very) simple region cost volume
Cost = (1 - V) - V;

% Segmentation parameters
Cen = single([35 50 35;   % Mesh centers
              70 50 35;
              50 50 70]);
R = single([10,10,10]);   % Mesh radii
nsub = 3;                 % Subdiv. lvl. for subdiv. icosahedron
nsamples = 50;            % Number of samples
step = 1;                 % Step size for samples
delta = 4;               % Smoothness parameter
costtype = 1;             % Cost type (1 = region costs)
Conn = ones(3) - eye(3);  % Mesh connectivity
Conn = logical(sparse(Conn));

% Segmentation without overlap constraints
[Fcs_olap,Vtx_olap] = mex_surfcut(Cost,Cen,R,nsub,nsamples,step,...
    delta,costtype);

% Segmentation with overlap constraints
[Fcs_nos,Vtx_nos] = mex_surfcut_planesep_qpbo(Cost,Cen,R,nsub,nsamples,...
    step,delta,Conn,costtype);

%% Display segmentation
figure(2); clf;

subplot(121);
patch(isosurface(V,0.5),'EdgeColor','none','FaceColor',[0.8,0.8,0.8],...
    'FaceAlpha',0.6);
hold on;
Colors = {'r','g','b'};
for i = 1:length(Fcs_olap)
    % Note that vertices need to be adjusted
    patch('Faces',Fcs_olap{i},'Vertices',Vtx_olap{i}(:,[2 1 3])+1,...
        'EdgeColor',Colors{i},'FaceColor','none');
end

axis equal;
axis([1 size(V,1) 1 size(V,1) 1 size(V,1)]);
view(75,15);
lighting phong
camlight right
material dull
title('Without overlap constraints');

subplot(122);
patch(isosurface(V,0.5),'EdgeColor','none','FaceColor',[0.8,0.8,0.8],...
    'FaceAlpha',0.6);
hold on;
Colors = {'r','g','b'};
for i = 1:length(Fcs_nos)
    % Note that vertices need to be adjusted
    patch('Faces',Fcs_nos{i},'Vertices',Vtx_nos{i}(:,[2 1 3])+1,...
        'EdgeColor',Colors{i},'FaceColor','none');
end

axis equal;
axis([1 size(V,1) 1 size(V,1) 1 size(V,1)]);
view(75,15);
lighting phong
camlight right
material dull
title('With overlap constraints');