% This script is attached to the companion paper Bernard et al. (2021):
% "The impact of lithology on fjord morphology"
%
% It extracts fjord widths and depths, and drainage area within the
% Scoresby sund area from 68°N to 75°N.
% The morphological analyses are performed on the BedMachine v3 DEM
% (Morlighem et al., 2017), which can be downloaded here:
% https://sites.uci.edu/morlighem/dataproducts/bedmachine-greenland/
%
% The 1:2 500 000e geological map (Henriksen, 2000), can be downloaded from
% the GEUS website: https://frisbee.geus.dk/webshop/?customer=nanoq&lang=en
%
% Entry parameters for the script are:
%    - Reference_elevation. To choose the reference elevation from which
%    the fjord width and depth are computed: 'Mean_Topography' or
%    'Sea_level'
%
% Along the script, results can be plotted. to do so, just uncomment the
% relevant line.
%
% Authors: Maxime Bernard, Philippe Steer - University of Rennes 1
clear all; close all;

%% INPUT PARAMETERS

%Define reference elevation from where to compute fjord width and depths
%Mean_topography or Sea-level
Reference_elevation = 'Mean_topography';


%% LOAD DATA

filename = 'BedMachineGreenland-2017-09-20.nc';
% ncdisp(filename)
BedrockError = ncread(filename,'errbed',[1 1],[Inf Inf],[1 1]);
x = ncread(filename,'x',1, Inf,3); %Decrease resolution to 450m
y = ncread(filename,'y',1, Inf, 3);
Z = ncread(filename,'bed',[1 1],[Inf Inf],[3 3]);
xsave=x;ysave=y;y=xsave;x=ysave;

%% MAKE COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zthreshold=0;%imagesc(x./1000,y./1000,Z);axis equal tight;caxis([Zthreshold max(Z(:))]);

% Define a map that filter long-wavelength topography to set a threshold
% for fjord width that does not depend on mean topography
% 45 km window
if strcmp(Reference_elevation,'Mean_topography')
    fun = @(x) mean(x(:));Zmean2 = nlfilter(Z,[100 100],fun);%imagesc(x./1e3,y./1e3,Zmean2);axis equal tight;
    Ztemp=Z;Z3=Ztemp;Z3(Z-Zmean2>0)=0;Z3(Z-Zmean2<=0)=1;%imagesc(x./1e3,y./1e3,Z3);axis equal tight;
elseif strcmp(Reference_elevation,'Sea_level') %Topography at sea-level
    Ztemp=Z;Z3=Ztemp;Z3(Z>0)=0;Z3(Z<=0)=1;%imagesc(x./1e3,y./1e3,Z3);axis equal tight;
end

% Crop image to the interesting part
indx=find(x>-2.4e6 & x<-1.5e6);indy=find(y> 0.36e6 & y< 0.865e6); %ScoresbySund area (65°N-75°N)

% save location into variables and define the mesh of the area
xloc=double(x(indx));yloc=double(y(indy));Zloc=Z(indy,indx);Z2loc=Z3(indy,indx);
[XLOC,YLOC]=meshgrid(xloc,yloc);
bed = Z(indy,indx);
Berr_loc = BedrockError(indy,indx);

% Compute drainage area
Zlocfill = fillsinks(Zloc);
DEM = GRIDobj(XLOC,YLOC,Zlocfill);
Ms = FLOWobj(DEM);
A = flowacc(Ms);As=flipud(A.Z);
As = As.*(abs((xloc(2)-xloc(1))).*(abs(yloc(2)-yloc(1))));
%figure;imagesc(xloc./1000,yloc./1000,log10(As.*1e-6));axis equal tight;

% Remove small valleys or ridges
temp1=abs(1-rmUnwantedRegion(abs(1-Z2loc),100));temp2=rmUnwantedRegion(temp1,200);Z2loc=temp2; 
%subplot(2,1,1);imagesc(xloc,yloc,Z2loc);axis equal tight;subplot(2,1,2);imagesc(xloc./1e3,yloc./1e3,temp2+Z2loc); axis equal tight;

% Create a skeleton
if strcmp(Reference_elevation,'Mean_topography')
    skel= bwmorph(Z2loc,'skel',Inf);skel1=bwmorph(skel, 'spur',2000);skel2=bwmorph(skel1, 'clean',500);skel=skel2;
elseif strcmp(Reference_elevation,'Sea_level')
    skel= bwmorph(Z2loc,'skel',inf);skelB=bwmorph(skel, 'clean',500);skel=skelB;
end
% figure
% imagesc(xloc,yloc,Z2loc);axis equal tight; hold on;
% scatter(XLOC(skel==1),YLOC(skel==1),10,'r','fill'); 

% Correct for small branches
B = bwmorph(skel, 'branchpoints'); E = bwmorph(skel, 'endpoints');
scatter(XLOC(B==1),YLOC(B==1),'gr','fill')
scatter(XLOC(E==1),YLOC(E==1),'b','fill'); hold off
[ytemp,xtemp] = find(E); B_loc = find(B);
Dmask = false(size(skel));
for k = 1:numel(xtemp)
    D = bwdistgeodesic(skel,xtemp(k),ytemp(k));
    distanceToBranchPt = min(D(B_loc));
    Dmask(D < distanceToBranchPt) =true;
end
Z2loc_centerline = logical(skel - Dmask);
% figure;
% imagesc(xloc./1e3,yloc./1e3,Z2loc+Z2loc_centerline);axis equal tight;

% Detect the borderline
Z2loc_borderline = bwmorph(Z2loc, 'remove');Z2loc_borderline = bwmorph(Z2loc_borderline, 'clean', inf);%imagesc(xloc,yloc,Z2loc_borderline);axis equal tight;
imagesc(xloc./1e3,yloc./1e3,Zloc);axis equal tight; hold on
%ind = find(double(Z2loc_borderline)==1);
% scatter(XLOC(ind)./1e3,YLOC(ind)./1e3,10,'fill');
% scatter(XLOC(Zloc<10 & Zloc>-10)./1e3,YLOC(Zloc<10 & Zloc>-10)./1e3,10,'fill');hold off;

%Get elevation along borderline
Zborder = Zloc.*double(Z2loc_borderline);

% Now get the distance transform to the centerline.
[Z2loc_dist,idist] = bwdist(Z2loc_borderline);
Z2loc_dist = Z2loc_dist.*abs(xloc(2)-xloc(1))./1000; %imagesc(xloc./1e3,yloc./1e3,Z2loc_dist);axis equal tight;

% Get the distance to the borderline along the centerline
Z2loc_centerline_dist=Z2loc_centerline.*Z2loc_dist;%imagesc(xloc./1e3,yloc./1e3,Z2loc_centerline_dist);axis equal tight;

% Keep only the maximum value in a radius of 2km around the centerpoints
fun = @(x) max(x(:));Z2loc_centerline_dist_max = Z2loc_centerline.*nlfilter(Z2loc_centerline_dist,[5 5],fun);
%imagesc(xloc,yloc,Z2loc_centerline_dist_max);axis equal tight;

% Determine the coastline coordinates 
yc = [4.7e5,5.25e5,5.539e5,6.24e5,6.817e5,7.27e5,7.88e5,8.46e5,8.2e5,7.69e5,7.551e5,7.24e5,7.37e5,7.32e5,7.42e5];xc = [-2.4e6,-2.36e6,-2.34e6,-2.31e6,-2.26e6,-2.22e6,-2.13e6,-1.973e6,-1.9e6,-1.812e6,-1.776e6,-1.735e6,-1.68e6,-1.57e6,-1.5e6];
yc2 = [3.6e5,4.28e5,5.027e5,4.68e5,4.4e5,5.8e5]; xc2 = [-2.34e6,-2.23e6,-2.02e6,-1.858e6,-1.65e6,-1.51e6];

% save coastline coordinates
yc_temp = yc; xc_temp = xc; yc2_temp = yc2; xc2_temp = xc2;

% Build the coastline vector
coast2 = [yc',xc'];coast3=[yc2',xc2'];
ycoastloc = interp1(coast2(:,2),coast2(:,1),xloc,'spline');xcoastloc=xloc;
ycoastloc2 = interp1(coast3(:,2),coast3(:,1),xloc,'spline');xcoastloc2=xloc;

% Create points along the centerline (not taking into account offshore locations)
clear  xc yc zc wc wc_max 
ind=find(Z2loc_centerline_dist>0);
for i=1:numel(ind) 
    [ind2,ind1]=ind2sub(size(Z2loc_centerline_dist_max),ind(i));
    if yloc(ind2)<ycoastloc(ind1) & yloc(ind2)>ycoastloc2(ind1)
        k=k+1;
        xc(k)=xloc(ind1);yc(k)=yloc(ind2);
        wc(k)=2.*Z2loc_centerline_dist(ind(i)).*1000;wc_max(k)=2.*Z2loc_centerline_dist_max(ind(i)).*1000;
        zc(k)=Zloc(ind(i));
    end
end
ind=find(xc==0 & yc==0);xc(ind)=[];yc(ind)=[];wc(ind)=[];wc_max(ind)=[];
zc(ind)=[];

%Clear boundaries
ind = find(xc>-1.533e6|xc<-2.376e6); xc(ind)=[];yc(ind)=[];wc(ind)=[];
wc_max(ind) =[]; zc(ind) = [];

% Plot the results
% imagesc(xloc,yloc,Zloc.*3);caxis([-1 1]),axis equal tight;hold on;
% plot(xcoastloc2,ycoastloc2,'r');xlim([min(xloc) max(xloc)]);ylim([min(yloc) max(yloc)]);
% scatter(xc_temp,yc_temp,'gr','fill'); 
% scatter(xc,yc,10,wc_max,'fill');plot(xcoastloc,ycoastloc,'r');xlim([min(xloc) max(xloc)]);ylim([min(yloc) max(yloc)]);hold off;


%% Load Geologic map

load('Geology2500k.mat');

% Extract geological units
Gloc_unique=unique(unique(G)); 

%Load Legend colors
load('Geol_legend.mat');

k=1;
Gloc2=zeros(size(G));label_vector = cell(numel(Gloc_unique),1);
for i=1:numel(Gloc_unique)
    ind = Gloc_unique(i);indinR = str2double(extractfield(R,'gm_label')); labelind = extractfield(R,'Legend_Tex');
    Gloc2(G==Gloc_unique(i))=k;
    indlabel = find(indinR == ind,1);
    label_vector{i} = labelind(indlabel);
    k=k+1;
end
xsave=xg;ysave=yg;xG=ysave;yG=xsave;Gloc2=Gloc2';[XG,YG]=meshgrid(xG,yG);
Gloc = griddata(XG,YG,double(Gloc2),double(XLOC),double(YLOC),'nearest');
Gloc_unique=unique(unique(Gloc));

%Get geol color map and local lithology
LocalLabel = cell(numel(Gloc_unique),1);
indxlegend= zeros(size(Gloc_unique));
Gloc1=Gloc;
for i=1:numel(Gloc_unique)
    ind = Gloc_unique(i);
    LocalLabel(i) = label_vector(ind);
    indxlegend(i)= ind;
    Gloc1(Gloc1==Gloc_unique(i)) = i;
end

% Create a simplified geologic map
% Classify by periods
Gloc3=zeros(size(Gloc));
ind0 = find(Gloc==81 | Gloc==82); Gloc3(ind0)=1; % not rock
ind10 =find(Gloc==1); Gloc3(ind10)=0; %Water
ind1 = find(Gloc==10 | Gloc==11 | Gloc>=20 & Gloc<=30 | Gloc==44 | Gloc==50 | Gloc==53 | Gloc==54 | Gloc>=57&Gloc<=80);  Gloc3(ind1)=2; %Ante-Devonian
ind2=find(Gloc>=39 & Gloc<=42);    Gloc3(ind2)=5; % metasediments
ind3 = find(Gloc==38); Gloc3(ind3)=6; %Devonian
ind4 = find(Gloc>=18&Gloc<=19 | Gloc==37); Gloc3(ind4)=7; %Permian and Carboniferous
ind5 = find(Gloc==17 | Gloc>=33&Gloc<=36 | Gloc==12 | Gloc==83); Gloc3(ind5)=8; %Juras-Trias
ind6 = find(Gloc==14 | Gloc==16 | Gloc==32 | Gloc==56); Gloc3(ind6)=9; %Cretaceous
ind7 = find(Gloc==51 | Gloc==55); Gloc3(ind7)=4; %Cenozoic intrusive rock
ind8 = find(Gloc==7 | Gloc==8 | Gloc==31 | Gloc==15 | Gloc==46 | Gloc==47 | Gloc==48 | Gloc==49 |Gloc==45); Gloc3(ind8)=3; %Cenozoic basalts
ind9 = find(Gloc>=2&Gloc<=6 | Gloc==43 | Gloc==52); Gloc3(ind9)=2; %Undifferentiated

% Classify by first order lithology
Gloc5=zeros(size(Gloc));
ind0 = find(Gloc==81 | Gloc==82); Gloc5(ind0)=1; % not rock
ind10 = find(Gloc==1); Gloc5(ind10)=0; %Water
ind1 = find(Gloc==10 | Gloc==11 | Gloc==44 | Gloc==50 | Gloc==53 | Gloc==54 | Gloc>=57&Gloc<=80 | Gloc==51 | Gloc==55);  Gloc5(ind1)=2;%Cenozoic intrusive and gneissic basement
ind2=find(Gloc>=39 & Gloc<=42);    Gloc5(ind2)=4; % metasediments
ind3 = find(Gloc==38 | Gloc>=18&Gloc<=28 | Gloc==37 | Gloc==17 | Gloc>=33&Gloc<=36 | Gloc==12 | Gloc==83 ...
    | Gloc==14 | Gloc==16 | Gloc==32 | Gloc==56 | Gloc==15| Gloc==30); Gloc5(ind3)=5; %sediments
ind7 = find(Gloc==7 | Gloc==8 | Gloc==31 | Gloc==46 | Gloc==47 ...
    | Gloc==48 | Gloc==49 |Gloc==45 | Gloc==29); Gloc5(ind7)=3; % basalts 
ind9 = find(Gloc>=2&Gloc<=6 | Gloc==43 | Gloc==52); Gloc5(ind9)=2; %Undifferentiated

% Create Color scales
%For first order lithology
geolmap5 = [ 80, 133, 188;       % water
            255, 255, 255;       % ice
            245, 222, 179;       % hard rocks
            255,  20, 147;       % volcanics 
            52, 136, 60;         % metasediments
            203, 140,  55;       % sediments                                                                                                            
           %125, 125, 125;      % Undifferentiated
            ]./255;   
        
%For lithology classified by periods
geolmap3 = [ 80, 133, 188;       % water
            255, 255, 255;       % ice
            245, 222, 179;       % Ante-Devonian
            255,  20, 147;       % CENOZOIC BASALTS 
            255,   0,   0;       % CENOZOIC INTRUSIVE ROCKS
            52, 136, 60;         % Late proterozoic sediments
            203, 140,  55;       % DEVONIAN 
            210,  64,  40;       % Carbo-Perm
            52, 178, 201;        % Trias-JURASSIC
            127, 198,  78;       % CRETACEOUS                                                                                                      
            %125, 125, 125;      % Undifferentiated
            ]./255;   

%imagesc(xloc,yloc,Gloc3);axis equal tight;colormap(geolmap3(1:end,:));caxis([0 8.9]);

% Determine local relief along the fjord centerline
zc_min = zeros(size(xc));
zc_max = zeros(size(xc));
zc_mean = zeros(size(xc));  
zc_minborder = zeros(size(xc));
rc = zeros(size(xc));
As_maxb = zeros(size(xc)); 
geol3 = zeros(size(xc));
geol5 = geol3;
ind_geolnan = [];
%Loop through centerpoints
for i=1:numel(xc)
    % Compute the distance from a centerpoint
    dist=sqrt((xc(i)-XLOC).^2+(yc(i)-YLOC).^2);
    % Create a mask identifying points of the DEM where
    % distance is lower than the valley width
    mask=zeros(size(dist));mask(dist<wc_max(i))=1;
    % Assign first order lithology to the centerpoint
    Gtemp=Gloc5.*mask;
    if isempty(Gtemp(Gtemp~=0 & Gtemp~=1))
        geol5(i) = NaN;
    else
        geol5(i)=max(Gtemp(Gtemp~=0 & Gtemp~=1));
    end
    % Assign the period to the centerpoint
    Gtemp=Gloc3.*mask;
    if isempty(Gtemp(Gtemp~=0))
        geol3(i) = 0;
        ind_geolnan=[ind_geolnan,i];
    else
        geol3(i)=max(Gtemp(Gtemp~=0));   
    end
    % Compute local relief around the centerpoint (valley width)
    Ztemp=Zloc.*mask; Gtemp=Gloc3.*mask; Zb=Zborder.*mask;
    zc_min(i) =min( min( Ztemp(Ztemp~=0))); %min bed elevation
    zc_max(i) =max( max( Ztemp(Ztemp~=0 &Gtemp==geol3(i)))); %max bed elevation 
    zc_minborder(i) =min( min(Zborder(Zb~=0))); %min elevation on the valley border
    rc(i)= zc_max(i)-zc_min(i); %local relief
    Astemp=As.*mask; % Drainage area around the centerpoint
    As_maxb(i) =max( max( Astemp(Astemp~=0))); %Max drainage area arount the centerpoint
end

% Convert drainage area to km²
As_max = As_maxb.*1e-6; %in km²

% Handle negative relief
ind_rc=find(rc<0);rc_old = rc;rc(ind_rc)=[]; wc_max(ind_rc)=[]; As_max(ind_rc)=[];
geol5(ind_rc)=[];geol3(ind_rc)=[];zc_min(ind_rc)=[];zc_max(ind_rc)=[];
zc_minborder(ind_rc)=[];
xc(ind_rc) = []; yc(ind_rc)=[];zc(ind_rc)=[];

% Handle no lithology
ind_rc=ind_geolnan;rc(ind_rc)=[]; wc_max(ind_rc)=[]; As_max(ind_rc)=[];
geol5(ind_rc)=[];geol3(ind_rc)=[];zc_min(ind_rc)=[];zc_max(ind_rc)=[];
xc(ind_rc) = []; yc(ind_rc)=[];zc(ind_rc)=[];zc_minborder(ind_rc)=[];

% save results in a temporary file to avoid re-computation
save temp

%% PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all; close all;load temp;

% Compute depth of fjords along centerline
% Minimum elevation between the two sides valley borders - Minimum
% elevation found along the valley width
Dv = zc_minborder-zc_min;

%Log drainage area
As_log = log10(As_max); 

% ------- Plot fjord depth along centerlines ----------------
% To plot with valleys above sea-level, uncomment the second 'ind'
figure;load grayC.mat; set(gcf,'units','normalized','outerposition',[0 0 1 1]);
ind = find(Dv>0 & geol5>1 & zc_min<0); %Keep points with depth > 0 and with lithology identified (no ice, or in water), keep only fjords
%ind = find(Dv>0 & geol5>1); %Keep points with depth > 0 and with lithology identified (no ice, or in water), glacial valleys and fjords
% Background topography
cmap = [grayC;jet]; data = Dv(ind)./1e3; colormap(cmap);
cdata=0.99*(Zloc-(min(Zloc(:))))/(max(Zloc(:))-(min(Zloc(:))));
imagesc(xloc./1000,yloc./1000,cdata);axis equal tight;hold on;colorbar;xlim([min(xloc) max(xloc)]./1000);ylim([min(yloc) max(yloc)]./1000);
caxis([0,2]); vt = [0,0.25,0.5,0.75,1];c=colorbar('Position',[0.9,0.1,0.01,0.3]);
set(c,'ylim',[0,1],'ytick',vt,'yticklabel',min(Zloc(:))+vt*(max(Zloc(:))-min(Zloc(:))));
set(get(c,'ylabel'),'string','Elevation (m)');
% Foreground data
mindata =0;maxdata = 2;
cdata=0.99*(data-mindata)/(maxdata-mindata+1e-16)+1;
scatter(xc(ind)./1000,yc(ind)./1000,20,cdata,'fill'); vt = [0,0.25,0.5,0.75,1];c2=colorbar('Position',[0.9,0.5,0.01,0.3]);
set(c2,'ylim',[1,2],'ytick',vt+1,'yticklabel',mindata+vt*(maxdata-mindata));
set(get(c2,'ylabel'),'string','Valley Depth (km)');
set(c,'ylim',[0,1]);set(c2,'ylim',[1,2]); 
xlabel('Northing (km)');ylabel('Easting (km)');
view([-90 90]);


% ------- Plot fjord width along centerlines  ----------------
% To plot with valleys above sea-level, uncomment the second 'ind'
figure; set(gcf,'units','normalized','outerposition',[0 0 1 1]);
ind = find(Dv>0 & geol5>1 & zc_min<0); %Keep points with depth > 0 and with lithology identified (no ice, or in water),keep only fjords
%ind = find(Dv>0 & geol5>1); %Keep points with depth > 0 and with lithology identified (no ice, or in water), glacial valleys and fjords
% Background topography
cmap = [grayC;jet]; data = wc_max(ind)./1e3; colormap(cmap);
cdata=0.99*(Zloc-(min(Zloc(:))))/(max(Zloc(:))-(min(Zloc(:))));
imagesc(xloc./1000,yloc./1000,cdata);axis equal tight;hold on;colorbar;xlim([min(xloc) max(xloc)]./1000);ylim([min(yloc) max(yloc)]./1000);
caxis([0,2]); vt = [0,0.25,0.5,0.75,1];c=colorbar('Position',[0.9,0.1,0.01,0.3]);
set(c,'ylim',[0,1],'ytick',vt,'yticklabel',min(Zloc(:))+vt*(max(Zloc(:))-min(Zloc(:))));
set(get(c,'ylabel'),'string','Elevation (m)');
% Foreground data
mindata =0;maxdata = 30;
cdata=0.99*(data-mindata)/(maxdata-mindata+1e-16)+1;
scatter(xc(ind)./1000,yc(ind)./1000,20,cdata,'fill'); vt = [0,0.25,0.5,0.75,1];c2=colorbar('Position',[0.9,0.5,0.01,0.3]);
set(c2,'ylim',[1,2],'ytick',vt+1,'yticklabel',mindata+vt*(maxdata-mindata));
set(get(c2,'ylabel'),'string','Valley Width (km)');
set(c,'ylim',[0,1]);set(c2,'ylim',[1,2]); 
xlabel('Northing (km)');ylabel('Easting (km)');
view([-90 90]);

% ------- Plot drainage area along centerlines  ----------------
% To plot with valleys above sea-level, uncomment the second 'ind'
figure; set(gcf,'units','normalized','outerposition',[0 0 1 1]);
ind = find(Dv>0 & geol5>1 & zc_min<0); %Keep points with depth > 0 , with lithology identified (no ice, or in water), keep only fjords
%ind = find(Dv>0 & geol5>1); %Keep points with depth > 0 and with lithology identified (no ice, or in water), glacial valleys and fjords
% Background topography
cmap = [grayC;jet]; data = log10(As_max(ind)); colormap(cmap);
cdata=0.99*(Zloc-(min(Zloc(:))))/(max(Zloc(:))-(min(Zloc(:))));
imagesc(xloc./1000,yloc./1000,cdata);axis equal tight;hold on;colorbar;xlim([min(xloc) max(xloc)]./1000);ylim([min(yloc) max(yloc)]./1000);
caxis([0,2]); vt = [0,0.25,0.5,0.75,1];c=colorbar('Position',[0.9,0.1,0.01,0.3]);
set(c,'ylim',[0,1],'ytick',vt,'yticklabel',min(Zloc(:))+vt*(max(Zloc(:))-min(Zloc(:))));
set(get(c,'ylabel'),'string','Elevation (m)');
% Foreground data
mindata =0;maxdata = max(data(:));
cdata=0.99*(data-mindata)/(maxdata-mindata+1e-16)+1;
scatter(xc(ind)./1000,yc(ind)./1000,20,cdata,'fill'); vt = [0,0.25,0.5,0.75,1];c2=colorbar('Position',[0.9,0.5,0.01,0.3]);
set(c2,'ylim',[1,2],'ytick',vt+1,'yticklabel',mindata+vt*(maxdata-mindata));
set(get(c2,'ylabel'),'string','log10 Drainage Area (km²)');
set(c,'ylim',[0,1]);set(c2,'ylim',[1,2]); 
xlabel('Northing (km)');ylabel('Easting (km)');
view([-90 90]);


% ------------------------  Plot Geology ----------------------------------
% Create color scale
geolmap_cb = zeros(256,3); ind=1:42; for i=1:6, geolmap_cb(ind+42*(i-1),:)=repmat(geolmap5(i,:),42,1); end
cmap = [geolmap_cb;geolmap_cb]; 
ind = find(Dv>0 & zc_min<0); %Keep points with depth > 0 , with lithology identified (no ice, or in water), keep only fjords
%ind = find(Dv>0); %Keep points with depth > 0 and with lithology identified (no ice, or in water), glacial valleys and fjords
data =geol5(ind);
% Legend for lithology
label_litho = {'Sea','Ice','hard rocks','volcanics','Metasediments','sediments'};
figure;set(gcf,'units','normalized','outerposition',[0 0 1 1]); colormap(cmap);
%Backgroung geology
cdata=0.97*(Gloc5-(min(Gloc5(:))))/(5-(min(Gloc5(:))));
p=pcolor(xloc./1000,yloc./1000,cdata);shading flat;axis equal tight;hold on;alpha(p,0.4)
caxis([0,2]); vt = [0:1/5.5:1];c=colorbar('Position',[0.65,0.1,0.01,0.3]);
set(c,'ylim',[0,1],'ytick',vt,'yticklabel',label_litho);xlim([min(xloc) max(xloc)]./1000);ylim([min(yloc) max(yloc)]./1000);
set(get(c,'ylabel'),'string','Geology');
view([90 90]);
set(gca,'Xdir','reverse','Fontsize',12);
xlabel('Northing (km)');ylabel('Easting (km)');
% Foregroung geology (centerline)
mindata=0; maxdata=5;
cdata=0.97*(data-mindata)/(maxdata-mindata)+1;
scatter(xc(ind)./1000,yc(ind)./1000,10,cdata,'fill'); vt = [0,0.25,0.5,0.75,1];c2=colorbar('Position',[0.65,0.5,0.01,0.3]);
set(c2,'ylim',[0,2],'ytick',vt+1,'yticklabel',mindata+vt*(maxdata-mindata),'Visible','off');
set(c,'ylim',[0,0.97]);


% --------- Scatter plots with simple geology (width vs relief) -----------
% Fjord only
% geol5: 2=BASEMENT - 3=METASEDIMENT - 4=Volcanics - 5=Sediment
figure,set(gcf,'units','normalized','outerposition',[0 0 1 1]);
ind=find(geol5==2&Dv~=0&zc_min<0); subplot(1,4,1);scatter(wc_max(ind)./1000,Dv(ind)./1e3,10,geol5(ind),'fill');colormap(geolmap5);caxis([0 5]);title('Basement'); xlim([0 50]); ylim([0 2.1]); axis square;xlabel('Valley width (km)');ylabel('Valley depth (km)');
ind=find(geol5==4&Dv~=0&zc_min<0); subplot(1,4,2);scatter(wc_max(ind)./1000,Dv(ind)./1e3,10,geol5(ind),'fill');colormap(geolmap5);caxis([0 5]);title('Metasediment'); xlim([0 50]); ylim([0 2.1]); axis square;xlabel('Valley width (km)');ylabel('Valley depth (km)');
ind=find((geol5==5)&Dv~=0&zc_min<0); subplot(1,4,4);scatter(wc_max(ind)./1000,Dv(ind)./1e3,10,geol5(ind),'fill');colormap(geolmap5);caxis([0 5]);title('Sedimentaries');xlim([0 50]); ylim([0 2.1]); axis square;xlabel('Valley width (km)');ylabel('Valley depth (km)');
ind=find(geol5==3&Dv~=0&zc_min<0); subplot(1,4,3);scatter(wc_max(ind)./1000,Dv(ind)./1e3,10,geol5(ind),'fill');colormap(geolmap5);caxis([0 5]);title('Volcanics'); xlim([0 50]); ylim([0 2.1]); axis square;xlabel('Valley width (km)');ylabel('Valley depth (km)');
xlim([0 50]); ylim([0 2.1])
xlabel('Valley width (km)');ylabel('Valley depth (km)');

% --------- Scatter plots with simple geology (width vs relief) -----------
% Fjords and glacial valleys above sea level
% geol5: 2=BASEMENT - 3=METASEDIMENT - 4=Volcanics - 5=Sediment
figure,set(gcf,'units','normalized','outerposition',[0 0 1 1]);
ind=find(geol5==2&Dv~=0); subplot(1,4,1);scatter(wc_max(ind)./1000,Dv(ind)./1e3,10,geol5(ind),'fill');colormap(geolmap5);caxis([0 5]);title('Basement'); xlim([0 50]); ylim([0 2.1]); axis square;xlabel('Valley width (km)');ylabel('Valley depth (km)');
ind=find(geol5==4&Dv~=0); subplot(1,4,2);scatter(wc_max(ind)./1000,Dv(ind)./1e3,10,geol5(ind),'fill');colormap(geolmap5);caxis([0 5]);title('Metasediment'); xlim([0 50]); ylim([0 2.1]); axis square;xlabel('Valley width (km)');ylabel('Valley depth (km)');
ind=find(geol5==5&Dv~=0); subplot(1,4,4);scatter(wc_max(ind)./1000,Dv(ind)./1e3,10,geol5(ind),'fill');colormap(geolmap5);caxis([0 5]);title('Sedimentaries');xlim([0 50]); ylim([0 2.1]); axis square;xlabel('Valley width (km)');ylabel('Valley depth (km)');
ind=find(geol5==3&Dv~=0); subplot(1,4,3);scatter(wc_max(ind)./1000,Dv(ind)./1e3,10,geol5(ind),'fill');colormap(geolmap5);caxis([0 5]);title('Volcanics'); xlim([0 50]); ylim([0 2.1]); axis square;xlabel('Valley width (km)');ylabel('Valley depth (km)');
xlim([0 50]); ylim([0 2.1])
xlabel('Valley width (km)');ylabel('Valley depth (km)');



%------------------ Power-law relationships -------------------------------
ind_As = find(As_log>4.8);
As_max2 = As_max;As_max2(ind_As) = [];
Dv2 = Dv; Dv2(ind_As) = []; geol5b = geol5; geol5b(ind_As) = [];zc_min2 = zc_min; zc_min2(ind_As) = [];
% Limit range for drainage area
ind = find(Dv2>0 & zc_min2<0); %Only Fjords
%ind = find(Dv2>0); %Fjord and valleys
labelax = {'Drainage area (km²)','Valley Depth (km)'};
% Compute and plot power law relationships
plot_results(geol5b(ind),Dv2(ind)./1e3,As_max2(ind),geolmap5,[3,3],labelax,'Area_Depth.xls'); 

%%%%%%%%% Scatter plot with simple geology (drainage area vs Width)%%%%%%%%%%%%%%%%%%%%%%%%%%                
labelax = {'Drainage area (km²)','Valley Width (km)'};
wc_max2 = wc_max; wc_max2(ind_As) = [];
plot_results(geol5b(ind),wc_max2(ind)./1e3,As_max2(ind),geolmap5,[1,3],labelax,'Area_Width.xls'); 
