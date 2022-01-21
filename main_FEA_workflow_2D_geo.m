%% Programm zur FEM Simulation von beliebigen 2D-Geometrien mit freiem CCX solver
% Verfasser:        Hossein Omid Beiki

%% initialization
if ~isdeployed
    clc
    clear
    close all
    fclose all;
end

%% input (example)
pointsX = [0, 100, 100, 0];
pointsY = [50, 50, 0, 0];

boundary_pointsX = [0, 10, 10, 0];
boundary_pointsY = [5, 5, 0, 0];

%% set the result path
resultpath = uigetdir(pwd,'Select the directory for the results!');
if resultpath == 0
    disp('Aborted! No directory for the results was selected.')
    return
end
iFileName = 'test';

% temporary files for communicating with the mesher
tmpOutFileName = fullfile(resultpath,'tmpout.json');
tmpInFileName = fullfile(resultpath,'tmpin.json');

%% set mesher and solver parameters and the material constants
youngsModulus  = 2.1e+8; % E-Modul [kg / mm . s^2]
poissonsRatio  = 0.3; % poissonsRatio [-]
massDensity = 7.83E-06; % [kg / mm^3]
gravitationalAcceleration = 9.81e+3; % [mm / s^2]

% --> displacement results are in [mm]

Mesh2D_GeometricOrder = 'quadratic'; % es werden quadratische dreieckige Elemente für das 2D-Netz verwendet
iThickness = 10;
Hmin = 2; % Hmin: Target minimum mesh edge length [mm]
% here the Hmin is set to the thickness of the corresponding part. this
% ensures that the mesh will not be too fine.
Hmax = 2; % Hmax: Target maximum mesh edge length [mm]
% here the Hmax is set to an empirical value related to the root of the part area, this
% ensures that the mesh will not be too rough.

% Solver mode: linear or nonlinear? 
% decide if the solver should consider nonlinear effects of large displacements.
% should account for geometric nonlinearity during the step.
nonlinear = 0;

OMP_NUM_THREADS = 2; % sets the number of threads to use for parallel regions for solver

resolution = 1; % for creating training data: for example: resolution = 1 means one pixel corresponds to 1 mm

doPlot = true;

%% process the input
pgon = polyshape(pointsX,pointsY, 'Simplify', true); % create a polyshape from the geometry
boundary_polyshape = polyshape(boundary_pointsX,boundary_pointsY, 'Simplify', true);

% plot the initial situation
figure;
plot(pgon);hold on;
plot(boundary_polyshape);
axis equal

% generate the mesh
% create & plot 2D mesh
% prepare the input argument of the mesher fuction
inputStruct = struct('pgon', pgon, 'Hmin',Hmin, 'Hmax',Hmax,...
    'GeometricOrder', Mesh2D_GeometricOrder, 'doPlot',false, 'filename', tmpOutFileName);

% creat/open the temp file
fid = fopen(tmpInFileName , 'w'); % w:overwrite
if fid == -1 % if error
end

fprintf(fid,'%s', jsonencode(inputStruct));
fclose(fid);
if exist(tmpOutFileName, 'file')==2
    delete(tmpOutFileName);
end

try
    [~,~,mesherWasSuccessful] = Mesher2D_PDE(tmpInFileName);
    %system([pwd '\Mesher2D_PDE\for_testing\Mesher2D_PDE.exe' ' ' tmpInFileName]);
    model = jsondecode(fileread(tmpOutFileName));
    if isempty(model)
        error('error reading mesher output file or there is no mesher output file')
    end
catch ME
    disp('mesher was NOT successful');
    disp(ME.message)
end
           
% prepair nodes and elements sets
nodesXY =  model.Mesh.Nodes;
elementsQuadTriangular =  model.Mesh.Elements;
nodes = [(1:size(nodesXY,2)) ;nodesXY; zeros(1,size(nodesXY,2))];
elements = [(1:size(elementsQuadTriangular,2)) ;elementsQuadTriangular];

id_boundary = isinterior(boundary_polyshape, nodesXY');

if doPlot
    h3 = figure;
    plot(nodesXY(1,:),nodesXY(2,:),'LineStyle','non',...
        'Marker','.','MarkerEdgeColor', 'b', 'DisplayName','Nodes')
    hold on
    plot(nodesXY(1,id_boundary), nodesXY(2,id_boundary),...
        'LineStyle','non','Marker','.','MarkerEdgeColor','r', 'DisplayName','Boundary')
    axis equal
    legend('show','Location','northeastoutside')
    close(h3)
end

%% write calculix inp data
fileID = fopen([iFileName '.inp'],'w');

% Definition von Knoten, Elementen und Sets
fprintf(fileID,'%s\n','**');
fprintf(fileID,'%s\n','*NODE');
fprintf(fileID,'%d, %f, %f, %f \n',nodes);

fprintf(fileID,'%s\n','**');
fprintf(fileID,'%s\n','*ELEMENT, TYPE=s6, ELSET=ELEMENTS_ALL');
fprintf(fileID,'%d, %d, %d, %d, %d, %d, %d\n',elements);

fprintf(fileID,'%s\n','**');
fprintf(fileID,'%s\n','*NSET, NSET=N_FIXED');
fprintf(fileID,'%d, \n',nodes(1,id_boundary));

% Materialdefinition und -zuweisung
fprintf(fileID,'%s\n','**');
fprintf(fileID,'%s\n','*MATERIAL, NAME=steel');
fprintf(fileID,'%s\n','*ELASTIC');
fprintf(fileID,'%d, %d\n', youngsModulus, poissonsRatio); % E-Modul [N/m^2] & Querkontraktionszahl [-]
fprintf(fileID,'%s\n','*DENSITY');
fprintf(fileID,'%d\n',massDensity); % [kg/m^3]
fprintf(fileID,'%s\n','**');
fprintf(fileID,'%s\n','*SHELL SECTION, ELSET=ELEMENTS_ALL, MATERIAL=steel');
fprintf(fileID,'%d\n',iThickness); % thickness

% STEPS
% static step
fprintf(fileID,'%s\n','**');
if nonlinearityGeom
    fprintf(fileID,'%s\n','*STEP, NLGEOM');
else
    fprintf(fileID,'%s\n','*STEP');
end
fprintf(fileID,'%s\n','*STATIC');

% see: http://web.mit.edu/calculix_v2.7/CalculiX/ccx_2.7/doc/ccx/node179.html
%fprintf(fileID,'%s\n','*CONTROLS,PARAMETERS=FIELD');
%fprintf(fileID,'%s\n','0.5,1,0.01,,0.02,1.e-5,1,1.e-8');

% Fesselungen
fprintf(fileID,'%s\n','**');
fprintf(fileID,'%s\n','*BOUNDARY,OP=NEW');
fprintf(fileID,'%s\n','N_FIXED,1,3');

fprintf(fileID,'%s\n','*DLOAD');
fprintf(fileID,'%s %f %s \n','ELEMENTS_ALL,GRAV,', gravitationalAcceleration, ',0.,0.,-1.'); % [N]

% Speicherung von Ergebnissen in der Ergebnis-Datei
fprintf(fileID,'%s\n','**');
fprintf(fileID,'%s\n','*NODE FILE');
fprintf(fileID,'%s\n','U');

fprintf(fileID,'%s\n','**');
fprintf(fileID,'%s\n','*END STEP');

fclose(fileID);

%% Berechnungsinitiierung mit CalculiX CrunchiX
try
    %Starte CalculiX CrunchiX mit Input-Datei
    system(['set OMP_NUM_THREADS=' num2str(OMP_NUM_THREADS) '&&.\CCX\ccx\ccx .\' iFileName]); % verwende 4 threads
    
    %Lösche nicht-benötigte Ergebnis-Dateien
    delete([iFileName '.cvg']);
    delete([iFileName '.sta']);
    delete([iFileName '.dat']);
    delete('spooles.out');
    
    % system(['del .\' fileName '.dat']);
    % system(['del .\spooles.out']);
    
    source = [fullfile(pwd, iFileName) '.inp'];
    iFemResultFile = [fullfile(femResultpath,iFileName) '.inp'];
    movefile( source, iFemResultFile)
    source = [ iFileName '.frd'];
    iFemResultFile = [fullfile(femResultpath,iFileName) '.frd'];
    movefile( source, iFemResultFile)
    
catch
    disp('error solving or moving')
end

pause(1) % pause to make sure, the data is moved to the target path and is now available

%% lese die output-Datei von ccx und extrahiere die Verschiebungsmatrix
% Öffne CalculiX Ergebnisdatei
fID = fopen([fullfile(femResultpath,iFileName) '.frd']);
if fID ==-1
    disp(['error reading result frd file: ' iFileName]);
end

[dispNodesValue,dispNodesCoo] = parseFEMresultCCX(fID); %interne Funktion, siehe unten
fclose(fID);

%% visualize the result
% to be completed


%% internal functions
function [dispNodesValue,dispNodesCoo] = parseFEMresultCCX(fID)
% some variables
look_for_nodes = '2C';     % Stichwort für den Abschnitt in dem Knoten Koordinaten zu finden sind
look_for_disp = 'DISP';    % Stichwort für den Abschnitt in dem die Knotenverschiebungen zu finden sind

% extract the #node
Nodes_Number_not_found = true;
while Nodes_Number_not_found
    tline = fgetl(fID);
    found = strfind(tline,look_for_nodes);
    if found~=0                                 % Starte, falls look_for gefunden ist
        tmp = textscan(tline(found+2:end),'%f %f', 'Delimiter',' ','MultipleDelimsAsOne',1);
        nodesNo = cell2mat(tmp(1));
        Nodes_Number_not_found = false;
    end
end

% extract the node coordinates and numbers
dispNodesCoo = zeros(nodesNo,5);
fomatSpec = '%f %f %f %f %f'; % five numbers as: -1 nodeId X Y Z
for j = 1 : nodesNo
    tline = fgetl(fID);
    tmp = textscan(tline,fomatSpec, 'Delimiter',' ','MultipleDelimsAsOne',1);
    dispNodesCoo(j,:) = cell2mat(tmp);
end

% extract the node displacement
dispNodesValue = zeros(nodesNo,5);
tline = 'dummy';                                % Speichert aktuell gescannte Zeile
while ischar(tline)                             % Starte while-Schleife, solange tline char ist
    tline = fgetl(fID);                         % Aktuelle Zeile auf tline speichern
    found = strfind(tline,look_for_disp);       % Suche in tline nach look_for
    if found~=0                                 % Starte, falls look_for gefunden ist
        for i = 1:4                             % Überlese die ersten vier Zeilen, nachdem Displacement Abschnitt gefunden wurde
            tline = fgetl(fID);
        end
        for i = 1 : nodesNo                     % Lese Daten für alle Knoten aus
            tline = fgetl(fID);                 % Speicher aktuelle Zeile
            tmp = textscan(tline,fomatSpec, 'Delimiter',' ','MultipleDelimsAsOne',1);
            dispNodesValue(i,:) = cell2mat(tmp);
        end
        tline=-1;                               % Beende while-Schleife
    end
end
end