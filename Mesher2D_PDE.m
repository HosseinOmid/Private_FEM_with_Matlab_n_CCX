%% Funktion zur Vernetzung einer Geometrie mit 2D-Triangular-Shell Elementen
% Autor: Hossein Omid Beiki
%
% Inputs:
%   inputFile: a jason / txt file in which data are stored in json format
%       inputStruct: matlab struct which is decoded from the json file. it
%       has the following fields:
%          pgon: polyshape of geometry
%          Hmin: Target minimum mesh edge length
%          Hmax: Target maximum mesh edge length
%          GeometricOrder: Element type, 'linear' or 'quadratic'.
%             In general, 'quadratic' elements produce more accurate solutions.
%          doPlot: determines if the mesh should be plotted
%          filename: name of the output file
% outputs:
%    nodesXY: Mesh nodes, returned as a matrix. nodesXY is a 2-by-Nn matrix,
%       where Nn is the number of nodes in the mesh.
%       Each column of Nodes contains the x, y coordinates for that mesh node.
%       2-D meshes have nodes at the mesh triangle corners for linear elements,
%       and at the corners and edge midpoints for 'quadratic' elements.
%
%    elementsQuadratTriangular: Mesh elements, returned as an M-by-Ne matrix,
%       where Ne is the number of elements in the mesh, and M is:
%          3 for 2-D triangles with 'linear' GeometricOrder
%          6 for 2-D triangles with 'quadratic' GeometricOrder
%
%    wasSuccessful: determines if the Mesher2D_PDE was successful

function [nodesXY, elementsQuadratTriangular, wasSuccessful] = Mesher2D_PDE(inputFile)
%% init
wasSuccessful = false;
nodesXY = [];
elementsQuadratTriangular = [];
%% read the input
try
    inputStruct = jsondecode(fileread(inputFile));
    pgon = polyshape(inputStruct.pgon.Vertices);
    Hmin = inputStruct.Hmin;
    Hmax = inputStruct.Hmax;
    GeometricOrder = inputStruct.GeometricOrder;
    doPlot = inputStruct.doPlot;
    filename = inputStruct.filename;
    
catch ME
    disp('error reading the input file')
    disp(ME.message)
    return
end

%% create & plot 2D mesh
try
    tr2d = triangulation(pgon);
    
    % creat a mesh using PDE
    model = createpde('structural','transient-planestress');
    geometryFromMesh(model,tr2d.Points',tr2d.ConnectivityList'); % generate mesh
    
    % refine the mesh using PDE
    % Assigning the min element size
    generateMesh(model,'Hmin',Hmin,'Hmax', Hmax,'GeometricOrder',GeometricOrder');
catch ME
    disp('error meshing with PDE')
    disp(ME.message)
    return
end

try
    if doPlot
        h2 = figure;
        subplot(1,2,1)
        triplot(tr2d); axis equal
        title('Triangulation','interpreter','latex')
        subplot(1,2,2)
        pdemesh(model)
        title('Mesh','interpreter','latex')
        saveas(h2, [filename '_h2.png'],'png')
        savefig(h2, [filename '_h2' '.fig'],'compact');
        close(h2)
    end
catch ME
    disp('error plotting')
    disp(ME.message)
    return
end

%% prepair nodes and elements sets and write the output
try
    fid = fopen(filename, 'w'); % w:overwrite
catch ME
    disp('error creating output file')
    disp(ME.message)
    return
end
if fid == -1
    disp('error creating output file')
    return
end

try
    fprintf(fid,'%s', jsonencode(model));
catch ME
    disp('error writing in the output file')
    disp(ME.message)
end
try
    fclose(fid);
catch ME
    disp('could not close the output file, please consider that the output file could be blocked by matlab and should be closed manually')
    disp(ME.message)
    return
end

%% if successful, then return the output
nodesXY = model.Mesh.Nodes;
elementsQuadratTriangular = model.Mesh.Elements;

wasSuccessful = true;
disp('Mesher was successful')

end