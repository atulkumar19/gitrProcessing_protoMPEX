% cleanup

%check yields_ntv
%check Ftridyn file

tic %time meshing - for high resolution this is a large cost
planes = planes_from_stl('helicon_cylinder_4p14m.stl',0.024);
refinementnum=4;
toc

nFaces = length(planes);
abcd = zeros(length(planes),4);
area = zeros(nFaces,1);
centroid = zeros(nFaces,3);
plane_norm = zeros(nFaces,1);
BCxBA = zeros(nFaces,1);
CAxCB = zeros(nFaces,1);
density = zeros(nFaces,1);

for i=1:length(planes)
    A = planes(i,1:3);
    B = planes(i,4:6);
    C = planes(i,7:9);
    
    AB = B-A;
    AC = C-A;
    BC = C-B;
    BA = A-B;
    CA = -AC;
    CB = -BC;
    
    norm1 = norm(AB);
    norm2 = norm(BC);
    norm3 = norm(AC);
    
    s = (norm1+norm2+norm3)/2;
    area(i) = sqrt(s*(s-norm1)*(s-norm2)*(s-norm3));
    normalVec = cross(AB,AC);
    normalVec = normalVec./norm(normalVec);
   
    d = -(dot(normalVec,A));
    
    abcd(i,:) = [normalVec,d];
    plane_norm(i) = norm(normalVec);
    
    BCxBA(i) = sign(dot(cross(BC,BA),normalVec));
    CAxCB(i) = sign(dot(cross(CA,CB),normalVec));
    centroid(i,:) = [1/3*(planes(i,1)+planes(i,4)+planes(i,7)), ...
                1/3*(planes(i,2)+planes(i,5)+planes(i,8)), ...
                1/3*(planes(i,3)+planes(i,6)+planes(i,9))];
              
end

%% Find a region you'd like to refine

helcenter=1.7462;% helicon
% helcenter=3.5; % ICH

helregion1=helcenter-0.2;
Upstream=0.5;
helregion2=helcenter+0.2;
Downstream=4.14;
%centroid=mean([planes(:,3),planes(:,6),planes(:,9)]');

% helicon region
y_faces = find(centroid(:,3) > helregion1 & centroid(:,3) < helregion2);

ret = plot_panes(planes,y_faces)

planes_refined = refine_planes(planes(y_faces,:),refinementnum);

ret = plot_panes(planes_refined,1:1:length(planes_refined))

%% Remove old planes and add new ones
planes(y_faces,:) = [];
planes = [planes;planes_refined];
ret = plot_panes(planes,1:1:length(planes))


%% end cap refinement

for i=1:length(planes)
    A = planes(i,1:3);
    B = planes(i,4:6);
    C = planes(i,7:9);
    
    AB = B-A;
    AC = C-A;
    BC = C-B;
    BA = A-B;
    CA = -AC;
    CB = -BC;
    
    norm1 = norm(AB);
    norm2 = norm(BC);
    norm3 = norm(AC);
    
    s = (norm1+norm2+norm3)/2;
    area(i) = sqrt(s*(s-norm1)*(s-norm2)*(s-norm3));
    normalVec = cross(AB,AC);
    normalVec = normalVec./norm(normalVec);
   
    d = -(dot(normalVec,A));
    
    abcd(i,:) = [normalVec,d];
    plane_norm(i) = norm(normalVec);
    
    BCxBA(i) = sign(dot(cross(BC,BA),normalVec));
    CAxCB(i) = sign(dot(cross(CA,CB),normalVec));
    centroid(i,:) = [1/3*(planes(i,1)+planes(i,4)+planes(i,7)), ...
                1/3*(planes(i,2)+planes(i,5)+planes(i,8)), ...
                1/3*(planes(i,3)+planes(i,6)+planes(i,9))];
              
end

y_faces = find(abs(abcd(:,3)) > 0.9);
%y_faces = find(centroid(:,3) > Downstream & centroid(:,3) < Upstream);

ret = plot_panes(planes,y_faces)

planes_refined = refine_planes(planes(y_faces,:),refinementnum);

ret = plot_panes(planes_refined,1:1:length(planes_refined))

% Remove old planes and add new ones
planes(y_faces,:) = [];
planes = [planes;planes_refined];
ret = plot_panes(planes,1:1:length(planes))

%% recalculate plane norms, etc,

for i=1:length(planes)
    A = planes(i,1:3);
    B = planes(i,4:6);
    C = planes(i,7:9);
    
    AB = B-A;
    AC = C-A;
    BC = C-B;
    BA = A-B;
    CA = -AC;
    CB = -BC;
    
    norm1 = norm(AB);
    norm2 = norm(BC);
    norm3 = norm(AC);
    
    s = (norm1+norm2+norm3)/2;
    area(i) = sqrt(s*(s-norm1)*(s-norm2)*(s-norm3));
    normalVec = cross(AB,AC);
    normalVec = normalVec./norm(normalVec);
   
    d = -(dot(normalVec,A));
    
    abcd(i,:) = [normalVec,d];
    plane_norm(i) = norm(normalVec);
    
    BCxBA(i) = sign(dot(cross(BC,BA),normalVec));
    CAxCB(i) = sign(dot(cross(CA,CB),normalVec));
    centroid(i,:) = [1/3*(planes(i,1)+planes(i,4)+planes(i,7)), ...
                1/3*(planes(i,2)+planes(i,5)+planes(i,8)), ...
                1/3*(planes(i,3)+planes(i,6)+planes(i,9))];
              
end

nFaces=length(planes);
%% Write to cfg file





%%
%Plasma input files
% plasma_file = 'assets/Beers_Helicon_3D_100kW_LowDensity_Plasma.txt';
% plasma_file = '../assets/Beers_Helicon_3D_100kW_HighDensity_Plasma.txt';
% plasma_file = 'assets/Beers_Helicon_3D_100kW_LowDensityHighTe_Plasma.txt';
%  plasma_file = '../assets/Beers_Helicon_3D_100kW_HighDensityHighTe_Plasma.txt';

plasma_file = '../assets/Beers_Helicon_3D_ModelPaperDensity_Plasma.txt';
%plasma_file='assets/Beers_Helicon_3D_HigherDensity_SheathEdge.txt';


%Model Paper voltage files
%surface_file = 'assets/Beers_Helicon_3D_ModelPaperDensity_SheathCenter_NonMag10.xlsx';
%surface_file = 'C:\Users\cxe\Documents\Proto-MPEX\COMSOLmodel\Helicon\Inputs\Beers_Inputs\ModelPaperDensity\Beers_Helicon_3D_ModelPaperDensity_SheathCenter_Mag10.xlsx';
 
%Surface voltage files
% surface_file = 'assets/Beers_Helicon_3D_100kW_LowDensity_HeliconData_LayerMiddle9.xlsx';
% surface_file = 'assets/Beers_Helicon_3D_100kW_LowDensity_HeliconData_LayerMiddle7_MagCase.xlsx';
% surface_file = 'assets/Beers_Helicon_3D_100kW_HighDensity_HeliconData_LayerMiddle6.xlsx';
 surface_file = '../assets/Beers_Helicon_3D_100kW_HighDensity_HeliconData_LayerMiddle11_MagCase.xlsx';
%surface_file='assets/Beers_Helicon_3D_HigherDensity_SheathCenter_Mag13.xlsx';

nP = 1e6;
radius = 0.06256;


materialZ = zeros(length(planes),1);
surfs = zeros(length(planes),1);

HeliconRefined=find(centroid(:,3) > helregion1 & centroid(:,3) < helregion2);
materialZ(HeliconRefined) = 74;
surfs(HeliconRefined)=1;

%%
plotSet=1:length(planes);
plotSet(find(surfs))=[];
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];

% [X,Y,Z] = refineXYZ(X,Y,Z,6)

figure
patch(transpose(X),transpose(Y),transpose(Z),'g','FaceAlpha',.3,'EdgeColor','none')%'none')
title({'GITR Geometry'})
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
hold on


plotSet=find(surfs);
X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
patch(transpose(X),transpose(Y),transpose(Z),'b','FaceAlpha',.3,'EdgeColor','none')%'none')

%%
inDir = ones(1,length(planes));
figure(10)
hold on

for i=1:length(surfs)
    ind = i;
    normal = -abcd(ind,1:3)/plane_norm(ind);
    l_normal = 0.01;
    normal = l_normal*normal;
%     centroid = [1/3*(planes(ind,1)+planes(ind,4)+planes(ind,7)), ...
%                 1/3*(planes(ind,2)+planes(ind,5)+planes(ind,8)), ...
%                 1/3*(planes(ind,3)+planes(ind,6)+planes(ind,9))];
    dot_product = dot(normal,centroid(i,:));
    
    if(sign(dot_product) == 1)
        inDir(ind) = -1;
    end
    
    if(abs(normal(3)) > 0.005)
        surfs(ind) = 0;
        materialZ(ind) = 0;
    end

    normal = inDir(ind)*normal;
%     quiver3(centroid(1),centroid(2),centroid(3),normal(1),normal(2),normal(3))
end

ProtoMPEX_data
helicon_data


[xx yy zz emag ne_surf] = helicon_read_xls(surface_file);

zz=zz+1.7462;
% zz=zz+3.5;

for i=1:length(centroid)
    distance = sqrt((centroid(i,1) - xx).^2 + (centroid(i,2) - yy).^2 + (centroid(i,3) - zz).^2);
    
    [M I] = min(distance);
    e_value(i) = emag(I);
    dens_value(i) = ne_surf(I);
    if  abs(centroid(i,3))> 1.7462+0.2 % 3.5+0.2 
        e_value(i) = 0;
        dens_value(i) = 0;
    end
end

fileID = fopen('gitrGeometryPointPlane3d.cfg','w');
fprintf(fileID,'geom = \n{ \n   x1 = [');
fprintf(fileID,'%5e',planes(1,1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,1));
end
fprintf(fileID,' ] \n   y1 = [');
fprintf(fileID,'%5e',planes(1,2));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,2));
end
fprintf(fileID,' ] \n   z1 = [');
fprintf(fileID,'%5e',planes(1,3));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,3));
end
fprintf(fileID,' ] \n   x2 = [');
fprintf(fileID,'%5e',planes(1,4));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,4));
end
fprintf(fileID,' ] \n   y2 = [');
fprintf(fileID,'%5e',planes(1,5));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,5));
end
fprintf(fileID,' ] \n   z2 = [');
fprintf(fileID,'%5e',planes(1,6));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,6));
end
fprintf(fileID,' ] \n   x3 = [');
fprintf(fileID,'%5e',planes(1,7));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,7));
end
fprintf(fileID,' ] \n   y3 = [');
fprintf(fileID,'%5e',planes(1,8));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,8));
end
fprintf(fileID,' ] \n   z3 = [');
fprintf(fileID,'%5e',planes(1,9));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',planes(i,9));
end
fprintf(fileID,' ] \n   a = [');
fprintf(fileID,'%5e',abcd(1,1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',abcd(i,1));
end

fprintf(fileID,' ] \n   b = [');
fprintf(fileID,'%5e',abcd(1,2));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',abcd(i,2));
end

fprintf(fileID,' ] \n   c = [');
fprintf(fileID,'%5e',abcd(1,3));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',abcd(i,3));
end
fprintf(fileID,' ] \n   d = [');
fprintf(fileID,'%5e',abcd(1,4));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',abcd(i,4));
end
fprintf(fileID,' ] \n   plane_norm = [');
fprintf(fileID,'%5e',plane_norm(1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',plane_norm(i));
end
% fprintf(fileID,' ] \n   ABxAC = [');
% fprintf(fileID,'%5e',ABxAC(1))
% for i=2:nFaces
% fprintf(fileID, ',')
% fprintf(fileID,'%5e',ABxAC(i))
% end
fprintf(fileID,' ] \n   BCxBA = [');
fprintf(fileID,'%5e',BCxBA(1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',BCxBA(i));
end
fprintf(fileID,' ] \n   CAxCB = [');
fprintf(fileID,'%5e',CAxCB(1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',CAxCB(i));
end
fprintf(fileID,' ] \n   area = [');
fprintf(fileID,'%5e',area(1,1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%5e',area(i));
end
fprintf(fileID,' ] \n   Z = [');
fprintf(fileID,'%f',materialZ(1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%f',materialZ(i));
end
fprintf(fileID,' ] \n   surface = [');
fprintf(fileID,'%i',surfs(1));
for i=2:nFaces
fprintf(fileID, ',');
fprintf(fileID,'%i',surfs(i));
end
fprintf(fileID,' ] \n   inDir = [');
fprintf(fileID,'%i',inDir(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,'%i',inDir(i));
end
fprintf(fileID,' ] \n   potential = [');
fprintf(fileID,'%5e',e_value(1));
for i=2:nFaces
    fprintf(fileID, ',');
    fprintf(fileID,'%5e',e_value(i));
end
            fprintf(fileID,' ] \n');
            fprintf(fileID,'periodic = 0;\n');
            fprintf(fileID,'theta0 = 0.0;\n');
            fprintf(fileID,'theta1 = 0.0\n');
            fprintf(fileID,'periodic_bc_x0 = 0.0;\n');
            fprintf(fileID,'periodic_bc_x1 = 0.0;\n');
            fprintf(fileID,'periodic_bc_x = 0;}\n');
            fclose(fileID)
            


e_value = e_value; %Conversion of E-field to potential
plotSet = 1:1:length(planes);
plotSet = find(surfs>0);

X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];

%%
THETA = atan2(Y,X);
ss=sum(sign(THETA),2);
swit = find((abs(ss) < 3) & (abs(THETA(:,1)) > 1));
THETA(swit,:) = ss(swit).*abs(THETA(swit,:) );
figure(101)
% patch(transpose(X),transpose(Y),transpose(Z),e_value,'FaceAlpha',.3,'EdgeColor','k')
patch(transpose(THETA),transpose(Z),0*transpose(Z),e_value(plotSet),'FaceAlpha',1,'EdgeColor','k')
colorbar
title({'Surface Potential [V]','Low Flux Case'})
xlabel('Theta [radian]') % x-axis label
ylabel('z [m]') % y-axis label
set(gca,'fontsize',16)

% Potential figure

figure
patch(transpose(Z),transpose(THETA),0*transpose(Z),e_value(plotSet),'FaceAlpha',1,'EdgeColor','k')
c=colorbar;
ylabel(c, 'Surface Potential [V]','FontSize',15)
title('Unmagnetized Case Surface Potential [V]')
xlabel('Z [m]')
ylabel('Angle along window (deg.)')
ylim([0 360])
xlim([-0.15 0.15])
set(gca,'fontsize',15)

%%

file = '../assets/ftridynBackground.nc';
ncid = netcdf.open(file,'NC_NOWRITE');
[dimname, nE] = netcdf.inqDim(ncid,0);
[dimname, nA] = netcdf.inqDim(ncid,1);
if strcmp(file,'../assets/ftridynBackground.nc')
[dimname, nS] = netcdf.inqDim(ncid,2);
else
    nS = 1;
end
energy = ncread(file,'E');
angle = ncread(file,'A');
spyld = ncread(file,'spyld');
rfyld = ncread(file,'rfyld');
cosxDist = ncread(file,'cosXDist');
cosxDistRef = ncread(file,'cosXDistRef');
cosyDist = ncread(file,'cosYDist');
% coszDist = ncread(file,'cosZDist');
eDist = ncread(file,'energyDist');
eDistRef = ncread(file,'energyDistRef');
eDistEgrid = ncread(file,'eDistEgrid');
eDistEgridRef = ncread(file,'eDistEgridRef');
phiGrid = ncread(file,'phiGrid');
thetaGrid = ncread(file,'thetaGrid');
thisEdistRef = reshape(eDistRef(:,1,:),length(eDistEgridRef),[]);
% figure(100)
% plot(eDistEgridRef,thisEdistRef)

spyld=reshape(spyld(:,:,1),length(angle),length(energy));

figure(113)
h = pcolor(energy,angle,spyld)
h.EdgeColor = 'none';
colorbar
set(gca,'ColorScale','log')
% set(gca, 'YDir', 'normal')
 set(gca, 'XScale', 'log')
title({'Sputtering Yield D on Al','As a Function of Energy and Angle'})
xlabel('E [eV]') % x-axis label
ylabel('Angle [degrees]') % y-axis label
set(gca,'fontsize',16)

% r_grid=rGrid;
% z_grid=zGrid;
% br_interp=Brnew;
% bz_interp=Bznew;

br_edge = interpn(r_grid,z_grid,br_interp',0*z_grid + radius, z_grid);
bz_edge = interpn(r_grid,z_grid,bz_interp',0*z_grid + radius, z_grid);
bfield_angle_z =acosd(br_edge./sqrt(br_edge.^2 + bz_edge.^2));
bfield_angle_z(find(bfield_angle_z>1))=65;
element_angle = interpn(z_grid,bfield_angle_z,centroid(:,3));
element_angle(find(element_angle >90)) = 180 - element_angle(find(element_angle >90));
%element_angle(find(element_angle >=89.9)) = 89;

%element_angle(:,:)=75;

load('yields_ntv_DonAl_extendedV.mat');

ion_dens_wall = dens_value;

%Y0=interpn(energy,angle,spyld',e_value,element_angle','pchip',0);

Y0 = interpn(narray,Tarray,Varray,yields,ion_dens_wall,0*ion_dens_wall+8,e_value,'pchip',0);

plotSet = 1:1:length(planes);

X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
figure(102)
patch(transpose(X),transpose(Y),transpose(Z),e_value,'FaceAlpha',1,'EdgeColor','none')
colorbar
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

ion_temp_wall = interpn(y,z,ion_temp,0.06256,0);
ion_dens_wall = interpn(y,z,dens,0.06256,0);
ion_temp_wall=ion_temp(1,1);
ion_dens_wall = dens_value;

k=1.38e-23*11604;
c_bar = sqrt(8*k*ion_temp_wall/pi/4/1.66e-27);
flux = 0.25*ion_dens_wall*c_bar;

plotSet = find(surfs>0);

X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];
% THETA = atan2(Y,X);
% ss=sum(sign(THETA),2);
% swit = find((abs(ss) < 3) & (abs(THETA(:,1)) > 1));
% THETA(swit,:) = ss(swit).*abs(THETA(swit,:) );

%%
THETA = 0*Z;
for ii=1:length(THETA)
    %     if ii==2457
%     if ii==4859
%         ii
%     end
    THETA(ii,:) = atan2d(Y(ii,:),X(ii,:));
    diff21 = abs(THETA(ii,2) - THETA(ii,1));
    diff31 = abs(THETA(ii,3) - THETA(ii,1));
    diff32 = abs(THETA(ii,3) - THETA(ii,2));
    
    maxdiff = max([diff21,diff31,diff32]);
    if maxdiff > 90
        if THETA(ii,1)>0
            % THETA(ii,1)=THETA(ii,1)+360;
            %     end
            if THETA(ii,2)<0
                THETA(ii,2)=THETA(ii,2)+360;
            end
            if THETA(ii,3)<0
                THETA(ii,3)=THETA(ii,3)+360;
            end
        end
        if THETA(ii,1)<0
            % THETA(ii,1)=THETA(ii,1)+360;
            %     end
            if THETA(ii,2)>0
                THETA(ii,2)=THETA(ii,2)-360;
            end
            if THETA(ii,3)>0
                THETA(ii,3)=THETA(ii,3)-360;
            end
        end
    end
end
THETA=THETA+180;

%%
figure(103)
patch(transpose(THETA),transpose(Z),0*transpose(Z),flux(plotSet).*Y0(plotSet),'FaceAlpha',1,'EdgeColor','k')
% patch(transpose(X),transpose(Y),transpose(Z),flux(plotSet).*Y0(plotSet),'FaceAlpha',1,'EdgeColor','k')

colorbar
title({'D eroded Al Flux [m^{-2}s^{-1}]','Low Flux Case'})
xlabel('Theta [radian]') % x-axis label
ylabel('z [m]') % y-axis label
set(gca,'fontsize',16)

xx = -0.06;
yy = 0.0;
zz = 1.7462;
distance = sqrt((centroid(:,1) - xx).^2 + (centroid(:,2) - yy).^2 + (centroid(:,3) - zz).^2);
[v i] = min(distance);


erosion = flux.*Y0.*area';

erosion_inds = find(erosion);
erosion_sub = erosion(erosion_inds);
erosion_sub_cdf = cumsum(erosion_sub);
erosion_rate=erosion_sub_cdf(end)
erosion_sub_cdf = erosion_sub_cdf./erosion_sub_cdf(end);

% plot(erosion_sub_cdf)

rand1 = rand(nP,1);

element = interp1([0, erosion_sub_cdf],0:1:length(erosion_sub_cdf),rand1);

element_ceil = ceil(element);
x_sample = [];
y_sample = [];
z_sample = [];
vx_sample = [];
vy_sample = [];
vz_sample = [];
m = 27;
% meanE = 4.0;
% v = sqrt(2*meanE*1.602e-19/m/1.66e-27);
nPoints = 200;
maxE = 20;
Eb = 3.39;%8.79;
a = 5;
E = linspace(0,maxE,nPoints);
dE = E(2);
thompson2 = a*(a-1)*E.*Eb^(a-1)./(E+Eb).^(a+1);
figure(234)
plot(E,thompson2)
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
title('Sample Particle Positions')

m=27;
ecdf = cumsum(thompson2);
ecdf = ecdf./ecdf(end);
rand1 = rand(1,nP);
randTheta = 2*pi*rand(1,nP);
randPhi = 0.5*pi*rand(1,nP);
Esamp = interp1(ecdf,E,rand1);
v = sqrt(2*Esamp*1.602e-19/m/1.66e-27)';
    vx = v'.*cos(randTheta).*sin(randPhi);
    vy = v'.*sin(randTheta).*sin(randPhi);
    vz = v'.*cos(randPhi);
buffer = 1e-5;
plotSet = 1:length(planes);

X = [planes((plotSet),1),planes((plotSet),4),planes((plotSet),7)];
Y = [planes((plotSet),2),planes((plotSet),5),planes((plotSet),8)];
Z = [planes((plotSet),3),planes((plotSet),6),planes((plotSet),9)];

particle_counts = histcounts(erosion_inds(element_ceil),0.5:1:(length(plane_norm)+0.5));
nP0 = 0;

for i=1:length(particle_counts)
    if particle_counts(i) > 0
% i
% ind = erosion_inds(element_ceil);
 x_tri = X(i,:);
 y_tri = Y(i,:);
 z_tri = Z(i,:);
    parVec = [x_tri(2) - x_tri(1), y_tri(2) - y_tri(1) , z_tri(2) - z_tri(1)];
    parVec = parVec./norm(parVec);
samples = sample_triangle(x_tri,y_tri,z_tri,particle_counts(i));

normal = inDir(i)*(-abcd(i,1:3)/plane_norm(i));

v_inds = nP0+1:nP0+particle_counts(i);

x_sample(v_inds) = samples(:,1) + buffer*normal(1);
y_sample(v_inds) = samples(:,2) + buffer*normal(2);
z_sample(v_inds) = samples(:,3) + buffer*normal(3);

parVec2 = cross(parVec,normal);

newV = vx(v_inds)'.*parVec + vy(v_inds)'.*parVec2 + vz(v_inds)'.*normal;
% newV = v(i)*normal;
vx_sample(v_inds) = newV(:,1);
vy_sample(v_inds) = newV(:,2);
vz_sample(v_inds) = newV(:,3);

nP0 = nP0+particle_counts(i)
    end
end

index_array = 1:1:nP;
index_array(randperm(length(index_array)));
x_sample = x_sample(index_array);
y_sample = y_sample(index_array);
z_sample = z_sample(index_array);
vx_sample = vx_sample(index_array);
vy_sample = vy_sample(index_array);
vz_sample = vz_sample(index_array);
%
figure(104)
hold on
theta_sample = atan2(y_sample,x_sample);
scatter3(x_sample,y_sample,z_sample,'k')
% scatter3(theta_sample,z_sample,0*z_sample,'r')
% quiver3(x_sample,y_sample,z_sample,vx_sample./10./v,vy_sample./10./v,vz_sample./10./v)
% quiver3(0*vx',0*vx',0*vx',vx'./v,vy'./v,vz'./v)

xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
title('Sample Particle Positions')

%%

ncid = netcdf.create(['./particle_source_helicon.nc'],'NC_WRITE')
 
dimP = netcdf.defDim(ncid,'nP',nP);

xVar = netcdf.defVar(ncid,'x','double',[dimP]);
yVar = netcdf.defVar(ncid,'y','double',[dimP]);
zVar = netcdf.defVar(ncid,'z','double',[dimP]);
vxVar = netcdf.defVar(ncid,'vx','double',[dimP]);
vyVar = netcdf.defVar(ncid,'vy','double',[dimP]);
vzVar = netcdf.defVar(ncid,'vz','double',[dimP]);

netcdf.endDef(ncid);
 
netcdf.putVar(ncid, xVar, x_sample);
netcdf.putVar(ncid, yVar, y_sample);
netcdf.putVar(ncid, zVar, z_sample);
netcdf.putVar(ncid, vxVar, vx_sample);
netcdf.putVar(ncid, vyVar, vy_sample);
netcdf.putVar(ncid, vzVar, vz_sample);

netcdf.close(ncid);


%% Functions
function ret = plot_panes(planes0,plotSet)

X = [planes0((plotSet),1),planes0((plotSet),4),planes0((plotSet),7)];
Y = [planes0((plotSet),2),planes0((plotSet),5),planes0((plotSet),8)];
Z = [planes0((plotSet),3),planes0((plotSet),6),planes0((plotSet),9)];

% [X,Y,Z] = refineXYZ(X,Y,Z,6)
planes0 = [X(:,1) Y(:,1) Z(:,1) X(:,2) Y(:,2) Z(:,2) X(:,3) Y(:,3) Z(:,3)];
figure
patch(transpose(X),transpose(Y),transpose(Z),'g','FaceAlpha',.3,'EdgeColor','k')%'none')


title({'GITR Geometry'})
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
ret = 1;
end

function planes = planes_from_stl(filename,hmax)
model = createpde;
importGeometry(model,filename);% Import STL file
figure(2)
pdegplot(model,'FaceLabels','on') %Plot stl 

tic %time meshing - for high resolution this is a large cost
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',hmax);% Options Hmax and Hmin can be set, linear order can also be used
figure
pdeplot3D(model,'FaceAlpha',0.5)

[p,e,t] = meshToPet(mesh);

nPoints = length(p);
nTets = length(t);

tess = transpose(t(1:4,:));%sort(t(1:4,:),1);
% all faces
faces=[tess(:,[1 2 3]);tess(:,[1 2 4]); ...
       tess(:,[1 3 4]);tess(:,[2 3 4])];


faces = sort(faces,2);
faces = sortrows(faces);
Y = diff(faces);
zeroRow = [0,0,0];
k = ismember(Y,zeroRow,'rows');
k2 = find(k~=0);

faces([k2;k2+1],:) = [];

C = faces;

planes = zeros(length(C), 9);

planes(1:length(C),:) = [transpose(p(1:3,C(:,1))),transpose(p(1:3,C(:,2))),transpose(p(1:3,C(:,3)))];

end

function planes_refined = refine_planes(planes,n)
X = [planes(:,1),planes(:,4),planes(:,7)];
Y = [planes(:,2),planes(:,5),planes(:,8)];
Z = [planes(:,3),planes(:,6),planes(:,9)];

[Xrefined, Yrefined, Zrefined] = refineXYZ(X,Y,Z,n);

planes_refined = [Xrefined(:,1) Yrefined(:,1) Zrefined(:,1) Xrefined(:,2) Yrefined(:,2) Zrefined(:,2) Xrefined(:,3) Yrefined(:,3) Zrefined(:,3)];
end
function [Xrefined, Yrefined, Zrefined] = refineXYZ(X,Y,Z,n)

for j=1:n
    Xrefined = zeros(2*length(X),3);
    Yrefined = zeros(2*length(X),3);
    Zrefined = zeros(2*length(X),3);
    for i=1:length(X)
        A = [X(i,1) Y(i,1) Z(i,1)];
        B = [X(i,2) Y(i,2) Z(i,2)];
        C = [X(i,3) Y(i,3) Z(i,3)];
        
        AB = B-A;
        AC = C-A;
        BC = C-B;
        
        norms =[norm(AB) norm(BC) norm(AC)];
        [maxVal maxInd] = max(norms);
        if maxInd ==1
            midPtAB = A + 0.5*AB;
            Xrefined(2*i-1,1) = A(1);
            Xrefined(2*i,1) = midPtAB(1);
            Xrefined(2*i-1,2) = midPtAB(1);
            Xrefined(2*i,2) = B(1);
            Xrefined(2*i-1,3) = C(1);
            Xrefined(2*i,3) = C(1);
            Yrefined(2*i-1,1) = A(2);
            Yrefined(2*i,1) = midPtAB(2);
            Yrefined(2*i-1,2) = midPtAB(2);
            Yrefined(2*i,2) = B(2);
            Yrefined(2*i-1,3) = C(2);
            Yrefined(2*i,3) = C(2);
            Zrefined(2*i-1,1) = A(3);
            Zrefined(2*i,1) = midPtAB(3);
            Zrefined(2*i-1,2) = midPtAB(3);
            Zrefined(2*i,2) = B(3);
            Zrefined(2*i-1,3) = C(3);
            Zrefined(2*i,3) = C(3);
        elseif maxInd ==2
            midptBC = B + 0.5*BC;
            Xrefined(2*i-1,1) = A(1);
            Xrefined(2*i,1) = A(1);
            Xrefined(2*i-1,2) = B(1);
            Xrefined(2*i,2) = midptBC(1);
            Xrefined(2*i-1,3) = midptBC(1);
            Xrefined(2*i,3) = C(1);
            Yrefined(2*i-1,1) = A(2);
            Yrefined(2*i,1) = A(2);
            Yrefined(2*i-1,2) = B(2);
            Yrefined(2*i,2) = midptBC(2);
            Yrefined(2*i-1,3) = midptBC(2);
            Yrefined(2*i,3) = C(2);
            Zrefined(2*i-1,1) = A(3);
            Zrefined(2*i,1) = A(3);
            Zrefined(2*i-1,2) = B(3);
            Zrefined(2*i,2) = midptBC(3);
            Zrefined(2*i-1,3) = midptBC(3);
            Zrefined(2*i,3) = C(3);
        elseif maxInd ==3
            midptAC = A + 0.5*AC;
            Xrefined(2*i-1,1) = A(1);
            Xrefined(2*i,1) = midptAC(1);
            Xrefined(2*i-1,2) = B(1);
            Xrefined(2*i,2) = B(1);
            Xrefined(2*i-1,3) = midptAC(1);
            Xrefined(2*i,3) = C(1);
            Yrefined(2*i-1,1) = A(2);
            Yrefined(2*i,1) = midptAC(2);
            Yrefined(2*i-1,2) = B(2);
            Yrefined(2*i,2) = B(2);
            Yrefined(2*i-1,3) = midptAC(2);
            Yrefined(2*i,3) = C(2);
            Zrefined(2*i-1,1) = A(3);
            Zrefined(2*i,1) = midptAC(3);
            Zrefined(2*i-1,2) = B(3);
            Zrefined(2*i,2) = B(3);
            Zrefined(2*i-1,3) = midptAC(3);
            Zrefined(2*i,3) = C(3);
        end
    end
    X = Xrefined;
    Y = Yrefined;
    Z = Zrefined;
end
end

function samples = sample_triangle(x,y,z,nP)
x_transform = x - x(1);
y_transform = y - y(1);
z_transform = z- z(1);

% figure(2)
% plot3([x_transform x_transform(1)],[y_transform y_transform(1)],[z_transform z_transform(1)])

v1 = [x_transform(2) y_transform(2) z_transform(2)];
v2 = [x_transform(3) y_transform(3) z_transform(3)];
v12 = v2 - v1;
normalVec = cross(v1,v2);

a1 = rand(nP,1);
a2 = rand(nP,1);

samples = a1.*v1 + a2.*v2;
% hold on
% scatter3(samples(:,1),samples(:,2),samples(:,3))
samples2x = samples(:,1) - v2(1);
samples2y = samples(:,2) - v2(2);
samples2z = samples(:,3) - v2(3);
samples12x = samples(:,1) - v1(1);
samples12y = samples(:,2) - v1(2);
samples12z = samples(:,3) - v1(3);
v1Cross = [(v1(2).*samples(:,3) - v1(3).*samples(:,2)) (v1(3).*samples(:,1) - v1(1).*samples(:,3)) (v1(1).*samples(:,2) - v1(2).*samples(:,1))];
v2 = -v2;
v2Cross = [(v2(2).*samples2z - v2(3).*samples2y) (v2(3).*samples2x - v2(1).*samples2z) (v2(1).*samples2y - v2(2).*samples2x)];
v12Cross = [(v12(2).*samples12z - v12(3).*samples12y) (v12(3).*samples12x - v12(1).*samples12z) (v12(1).*samples12y - v12(2).*samples12x)];

v1CD = normalVec(1)*v1Cross(:,1) + normalVec(2)*v1Cross(:,2) + normalVec(3)*v1Cross(:,3);
v2CD = normalVec(1)*v2Cross(:,1) + normalVec(2)*v2Cross(:,2) + normalVec(3)*v2Cross(:,3);
v12CD = normalVec(1)*v12Cross(:,1) + normalVec(2)*v12Cross(:,2) + normalVec(3)*v12Cross(:,3);

inside = abs(sign(v1CD) + sign(v2CD) + sign(v12CD));
insideInd = find(inside ==3);
notInsideInd = find(inside ~=3);
% scatter3(samples(insideInd,1),samples(insideInd,2),samples(insideInd,3))

 v2 = -v2;
dAlongV1 = v1(1).*samples(notInsideInd,1) + v1(2).*samples(notInsideInd,2) + v1(3).*samples(notInsideInd,3);
dAlongV2 = v2(1).*samples(notInsideInd,1) + v2(2).*samples(notInsideInd,2) + v2(3).*samples(notInsideInd,3);

dV1 = norm(v1);
dV2 = norm(v2);
halfdV1 = 0.5*dV1;
halfdV2 = 0.5*dV2;

samples(notInsideInd,:) = [-(samples(notInsideInd,1) - 0.5*v1(1))+0.5*v1(1) ...
    -(samples(notInsideInd,2) - 0.5*v1(2))+0.5*v1(2) ...
    -(samples(notInsideInd,3) - 0.5*v1(3))+0.5*v1(3)];
% samples(notInsideInd,:) = [-(samples(notInsideInd,1) - 0.5*v2(1))+0.5*v2(1) ...
%     -(samples(notInsideInd,2) - 0.5*v2(2))+0.5*v2(2) ...
%     -(samples(notInsideInd,3) - 0.5*v2(3))+0.5*v2(3)];
samples(notInsideInd,:) = [(samples(notInsideInd,1) + v2(1)) ...
    (samples(notInsideInd,2) +v2(2)) ...
    (samples(notInsideInd,3) + v2(3))];
% figure(4)
% plot3([x_transform x_transform(1)],[y_transform y_transform(1)],[z_transform z_transform(1)])
% hold on
% scatter3(samples(:,1),samples(:,2),samples(:,3))

samples(:,1) = samples(:,1)+ x(1);
samples(:,2) = samples(:,2)+ y(1);
samples(:,3) = samples(:,3)+ z(1);

% figure(5)
% plot3([x x(1)],[y y(1)],[z z(1)])
% hold on
% scatter3(samples(:,1),samples(:,2),samples(:,3))
end
