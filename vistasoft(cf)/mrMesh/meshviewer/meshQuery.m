function mesh = meshQuery(mesh,summaryParams)
%
%   mesh = meshQuery(mesh,summaryParams)
%
%Author: BW
%Purpose:
%   Utility to set mesh parameters
%

if(summaryParams)
    % Get a partial list of key parameters.  This is the default.
    
    prompt = {'Decimate (0-1, 1=none):',...
            'Decimate Iterations:',...
            'Decimate Subiterations:',...
            'Smooth with Windowed Sinc (0|1):',...
            'Smooth Iterations:',...
            'Smooth Relaxation (0-2):',...
            'Do pre-decimation smoothing:',...
            'Curvature Mod Depth (0-1):',...
            'Curvature Color (0|1):'};
    defAns = {num2str(mesh.decimate_reduction),...
            num2str(mesh.decimate_iterations),...
            num2str(mesh.decimate_subiterations),...
            num2str(mesh.smooth_sinc_method),...
            num2str(mesh.smooth_iterations),...
            num2str(mesh.smooth_relaxation),...
            num2str(mesh.smooth_pre),...
            num2str(mesh.curvature_mod_depth),...
            num2str(mesh.curvature_color)};
    
    resp = inputdlg(prompt, 'Set Mesh Build Parameters', 1, defAns);
    
    if(~isempty(resp))
        mesh.decimate_reduction = str2num(resp{1});
        mesh.decimate_iterations = str2num(resp{2});
        mesh.decimate_subiterations = str2num(resp{3});
        mesh.smooth_sinc_method = str2num(resp{4});
        mesh.smooth_iterations = str2num(resp{5});
        mesh.smooth_relaxation = str2num(resp{6});
        mesh.smooth_pre = str2num(resp{7});
        mesh.curvature_mod_depth = str2num(resp{8});
        mesh.curvature_color = str2num(resp{9});
    else
        mesh = [];
        lights = [];
        return;
    end
    
else
    % Get the whole damn list of parameters.
    prompt = {'Decimate (0-1, 1=none):',...
            'Decimate Iterations:',...
            'Decimate Subiterations:',...
            'Decimate Preserve Edges (0|1):',...
            'Decimate Preserve Topology (0|1):',...
            'Decimate Boudary Vertex Delete (0|1):',...
            'Decimate Aspect Ratio (?):',...
            'Decimate Degree (?):',...
            'Smooth with Windowed Sinc (0|1):',...
            'Smooth Iterations:',...
            'Smooth Relaxation (0-1):',...
            'Smooth Feature Angle:',...
            'Smooth Edge Angle:',...
            'Smooth Boundary (0|1):',...
            'Smooth Feature Angle Smoothing On (0|1):',...
            'Do pre-decimation smoothing:',...
            'Curvature Mod Depth (0-1):',...
            'Curvature Color (0|1):'};
    defAns = {num2str(mesh.decimate_reduction),...
            num2str(mesh.decimate_iterations),...
            num2str(mesh.decimate_subiterations),...
            num2str(mesh.decimate_preserve_edges),...
            num2str(mesh.decimate_preserve_topology),...
            num2str(mesh.decimate_boudary_vertex_deletion),...
            num2str(mesh.decimate_aspect_ratio),...
            num2str(mesh.decimate_degree),...
            num2str(mesh.smooth_sinc_method),...
            num2str(mesh.smooth_iterations),...
            num2str(mesh.smooth_relaxation),...
            num2str(mesh.smooth_feature_angle),...
            num2str(mesh.smooth_edge_angle),...
            num2str(mesh.smooth_boundary),...
            num2str(mesh.smooth_feature_angle_smoothing),...
            num2str(mesh.smooth_pre),...
            num2str(mesh.curvature_mod_depth),...
            num2str(mesh.curvature_color)};
    
    resp = inputdlg(prompt, 'Set Mesh Build Parameters', 1, defAns);
    
    if(~isempty(resp))
        mesh.decimate_reduction = str2num(resp{1});
        mesh.decimate_iterations = str2num(resp{2});
        mesh.decimate_subiterations = str2num(resp{3});
        mesh.decimate_preserve_edges = str2num(resp{4});
        mesh.decimate_preserve_topology = str2num(resp{5});
        mesh.decimate_boudary_vertex_deletion = str2num(resp{6});
        mesh.decimate_aspect_ratio = str2num(resp{7});
        mesh.decimate_degree = str2num(resp{8});
        mesh.smooth_sinc_method = str2num(resp{9});
        mesh.smooth_iterations = str2num(resp{10});
        mesh.smooth_relaxation = str2num(resp{11});
        mesh.smooth_feature_angle = str2num(resp{12});
        mesh.smooth_edge_angle = str2num(resp{13});
        mesh.smooth_boundary = str2num(resp{14});
        mesh.smooth_feature_angle_smoothing = str2num(resp{15});
        mesh.smooth_pre = str2num(resp{16});
        mesh.curvature_mod_depth = str2num(resp{17});
        mesh.curvature_color = str2num(resp{18});
    else
        mesh = [];
        lights = [];
        return;
    end
end

% *** FIX THIS- should allow the relaxtion to use a different method, if
% desired.
if(mesh.smooth_sinc_method==1)
    relaxFactor = .0001;
else
    relaxFactor = 1.0;
end

return;