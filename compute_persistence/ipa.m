function ipa(in_file, out_file, max_distance)
% IPA The incremental persistence algorithm performs the following computations 
% on a 3D point cloud dataset: 
%
% 1) Approximates the topology of the input space using a Vietoris-Rips complex 
%
% 2) Computes the persistence of the zeroth and first homology groups using the 
%    algorithms described in "Topological Persistence and Simplification" by
%    Edelsbrunner et al.
%
close all;
pcd_data = loadpcd(in_file);

x = []; y = []; z = [];
r = []; g = []; b = [];

fprintf('Extracting positions and colors ... ');
if size(pcd_data, 3) > 1  % Organized point cloud 
    [m,n] = size(pcd_data(:,:,1)) 
    k = 1;
    for i=1:m
        for j=1:n
            if ~isnan(pcd_data(i,j,1))
                x(k) = pcd_data(i,j,1);
                y(k) = pcd_data(i,j,2);
                z(k) = pcd_data(i,j,3);
                r(k) = pcd_data(i,j,4);
                g(k) = pcd_data(i,j,5);
                b(k) = pcd_data(i,j,6);
                k = k + 1;
            end
        end
    end
else  % Unorganized point cloud 
    for i=1:size(pcd_data, 2)
        x(i) = pcd_data(1,i);
        y(i) = pcd_data(2,i);
        z(i) = pcd_data(3,i);
        r(i) = pcd_data(4,i);
        g(i) = pcd_data(5,i);
        b(i) = pcd_data(6,i);
    end
end
fprintf('done\n');

fprintf('Sorting %d points ... ', size(x,2));
tic; 
data = sortrows([x' y' z' r' g' b'], 1:3);
elapsed_time = toc;
fprintf('done: %.3f secs\n', elapsed_time);

fprintf('Initializing data structure ... ');
tic;
parfor i=1:size(data,1)
    point(i).x                  = data(i,1);
    point(i).y                  = data(i,2);
    point(i).z                  = data(i,3);
    point(i).r                  = data(i,4);
    point(i).g                  = data(i,5);
    point(i).b                  = data(i,6);
    point(i).parent             = i;
    point(i).death_idx          = 0;
    point(i).cycle              = [];
    point(i).neighbors          = [];
    edge_list(i).edges          = [];
    edge_list(i).offset         = 0;
    edge_list(i).n_edges        = 0;
    triangle_list(i).triangles  = [];
    triangle_list(i).n_triangles = 0;
end
elapsed_time = toc;
fprintf('done: %.3f secs\n', elapsed_time);

fprintf('Initializing k-d tree ... ');
tic; 
tree = KDTreeSearcher(data(:,1:3));
elapsed_time = toc;
fprintf('done: %.3f secs\n', elapsed_time);

fprintf('Finding neighbors within distance %.3f ... ', max_distance);
n_neighbors = 0;
tic; 
parfor i=1:size(data, 1)
    query = [point(i).x point(i).y point(i).z];
    idx = rangesearch(tree, query, max_distance);
    k = 1;
    for j=1:size(idx{1}, 2)
        if i < idx{1}(:,j) 
            point(i).neighbors(k) = idx{1}(:,j);
            k = k + 1;
        end
    end
    point(i).neighbors = sort(point(i).neighbors);
    n_neighbors = n_neighbors + size(point(i).neighbors, 2);
end
neighbors_elapsed_time = toc;
fprintf('done: %.1f neighbors per point, %.3f secs\n', n_neighbors / size(data, 1), neighbors_elapsed_time);

fprintf('Computing edges ... ');
total_edges = 0;
tic; 
parfor i=1:size(data, 1)
    n_edges = 0;
    for j=1:size(point(i).neighbors, 2)
        edge_list(i).edges = [edge_list(i).edges; i point(i).neighbors(j)];
        n_edges = n_edges + 1;
        total_edges = total_edges + 1;
    end
    edge_list(i).n_edges = n_edges;
end
edge_count = edge_list(1).n_edges; 
for i=2:size(data, 1)
    edge_list(i).offset = edge_count;
    edge_count = edge_count + edge_list(i).n_edges;
end
edge_list_elapsed_time = toc;
fprintf('done: %d edges computed, %.3f secs\n', total_edges, edge_list_elapsed_time);

fprintf('Constructing edge list ... ');
k = 1;
n_points = size(data, 1);
tic;
for i=1:size(data, 1)
    for j=1:size(edge_list(i).edges, 1)
        edge(k).u = edge_list(i).edges(j, 1);
        edge(k).v = edge_list(i).edges(j, 2);
        edge(k).idx = n_points + k;
        edge(k).cycle = [];
        edge(k).negative = 0;
        edge(k).death_idx = 0;
        k = k + 1;
    end
end
edges_elapsed_time = toc;
fprintf('done: %d edges processed, %.3f secs\n', size(edge, 2), edges_elapsed_time);

fprintf('Computing triangles ... ');
total_triangles = 0;
tic; 
parfor i=1:size(data, 1)
    n_triangles = 0;
    for j=1:edge_list(i).n_edges
        u = edge_list(i).edges(j,1);
        v = edge_list(i).edges(j,2);
        for k=1:size(point(i).neighbors, 2)
            w = point(i).neighbors(k);
            if w > v
                % If (u,w) and (v,w) are edges then we have a triangle
                [uw, uw_offset] = isEdge(w, point(u).neighbors);
                [vw, vw_offset] = isEdge(w, point(v).neighbors);
                if uw && vw 
                    uv_idx = edge_list(u).offset + j;
                    uw_idx = edge_list(u).offset + uw_offset;
                    vw_idx = edge_list(v).offset + vw_offset;
                    triangle_list(i).triangles = [triangle_list(i).triangles; u v w uv_idx uw_idx vw_idx];
                    n_triangles = n_triangles + 1;
                    total_triangles = total_triangles + 1;
                end
            end
        end
    end
    triangle_list(i).n_triangles = n_triangles;
end
triangle_list_elapsed_time = toc;
fprintf('done: %d triangles computed, %.3f secs\n', total_triangles, triangle_list_elapsed_time);

fprintf('Constructing triangle list ... ');
k = 1;
n_points = size(data, 1);
tic;
for i=1:size(data, 1)
    for j=1:size(triangle_list(i).triangles, 1)
        triangle(k).u = triangle_list(i).triangles(j, 1);
        triangle(k).v = triangle_list(i).triangles(j, 2);
        triangle(k).w = triangle_list(i).triangles(j, 3);
        triangle(k).uv_idx = triangle_list(i).triangles(j, 4);
        triangle(k).uw_idx = triangle_list(i).triangles(j, 5);
        triangle(k).vw_idx = triangle_list(i).triangles(j, 6);
        triangle(k).idx = n_points + total_edges + k;
        triangle(k).negative = 0;
        k = k + 1;
    end
end
triangles_elapsed_time = toc;
fprintf('done: %d triangles processed, %.3f secs\n', size(triangle, 2), triangles_elapsed_time);
fprintf('Total time: %.3f secs\n\n', neighbors_elapsed_time + edges_elapsed_time + ...
         edge_list_elapsed_time + triangles_elapsed_time + triangle_list_elapsed_time);

fprintf('Computing persistence of 0-cycles ... ');
negative_edges = 0;
tic;
for i=1:size(edge, 2)
    u = edge(i).u;
    v = edge(i).v;
    lambda = [v u];

    youngest = max(lambda);
    if isempty(point(youngest).cycle)
        point(youngest).cycle = lambda;
        point(youngest).death_idx = edge(i).idx;
        edge(i).negative = 1;
        negative_edges = negative_edges + 1;
    else
        while 1
            lambda_i = point(youngest).cycle;
            lambda = setdiff(union(lambda, lambda_i), intersect(lambda, lambda_i)); 
            if isempty(lambda) break; end
            youngest = max(lambda);
            if isempty(point(youngest).cycle)
                point(youngest).cycle = lambda;
                point(youngest).death_idx = edge(i).idx;
                edge(i).negative = 1;
                negative_edges = negative_edges + 1;
                break;
            end
        end
    end
end
elapsed_time = toc;
betti_0 = n_points - negative_edges;
fprintf('done: Betti 0 = %d, %d deaths, %.3f secs\n', betti_0, negative_edges, elapsed_time);

fprintf('Computing persistence of 1-cycles ... ');
negative_triangles = 0;
tic;
for i=1:size(triangle, 2)
    uv_idx = triangle(i).uv_idx;
    uw_idx = triangle(i).uw_idx;
    vw_idx = triangle(i).vw_idx;
    uv = edge(uv_idx);
    uw = edge(uw_idx);
    vw = edge(vw_idx);
    lambda = [];

    if uv.negative == 0
        lambda = [lambda uv_idx];
    end
    if uw.negative == 0
        lambda = [lambda uw_idx];
    end
    if vw.negative == 0
        lambda = [lambda vw_idx];
    end

    youngest = max(lambda);
    if isempty(edge(youngest).cycle)
        edge(youngest).cycle = lambda;
        edge(youngest).death_idx = triangle(i).idx;
        triangle(i).negative = 1;
        negative_triangles = negative_triangles + 1;
    else
        while 1
            lambda_i = edge(youngest).cycle;
            lambda = setdiff(union(lambda, lambda_i), intersect(lambda, lambda_i)); 
            if isempty(lambda) break; end
            youngest = max(lambda);
            if isempty(edge(youngest).cycle)
                edge(youngest).cycle = lambda;
                edge(youngest).death_idx = triangle(i).idx;
                triangle(i).negative = 1;
                negative_triangles = negative_triangles + 1;
                break;
            end
        end
    end
end
elapsed_time = toc;
positive_edges = size(edge, 2) - negative_edges;
betti_1 = positive_edges - negative_triangles;
fprintf('done: Betti 1 = %d, %d deaths, %.3f secs\n\n', betti_1, negative_triangles, elapsed_time);

file_id = fopen(out_file, 'w');
fprintf('Writing birth and death times to %s ... ', out_file);
for i=1:size(data, 1)
    if (point(i).death_idx == 0)
        fprintf(file_id, '%d inf\n', i);
    else
        fprintf(file_id, '%d %d\n', i, point(i).death_idx);
    end
end
for i=1:size(edge, 2)
    if (edge(i).negative == 0)
        if (edge(i).death_idx == 0)
            fprintf(file_id, '%d inf\n', edge(i).idx);
        else
            fprintf(file_id, '%d %d\n', edge(i).idx, edge(i).death_idx);
        end
    end
end
fclose(file_id);
fprintf('done\n');

end

function [value, offset] = isEdge(p, neighbors)
    value = 0;
    offset = 0;
    for i=1:size(neighbors, 2)
        if neighbors(i) == p
            value = 1;
            offset = i;
            break
        end
    end
end
