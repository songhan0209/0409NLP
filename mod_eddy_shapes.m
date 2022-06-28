function mod_eddy_shapes(nc_u,nc_v,nc_dim,rad,a,path_out)
% mod_eddy_shapes.m
%  mod_eddy_shapes(nc_u,nc_v,nc_dim,rad,a,path_out) computes the shapes of
%  the eddies identified by the centers detected with mod_eddy_centers.m
%  and saves them in [path_out,'eddy_shapes'];
%
%  - nc_u and nc_v are the time series of the 2D u and v velocity field;
%  - nc_dim defines the domain dimensions;
%  - rad defines the area where the streamfunction (PSI) is first computed
%  - a is one of the parameters from the eddy detection constraints
%  - path_out is the path where the algorithm output is saved
%
%  (For a description of the input parameters see param_eddy_tracking.m)
%
%  Eddy shapes are saved in the structure array [path_out,'eddy_shapes']:
%
%  - shapes(t).lonlat(n) : lonlat is a [2xm double] cell array, containing
%    the longitude and latitude position of the m vertices that define the
%    boundaries of the n-th eddy of the t-th day;
%
%  Another file which contains information on the process to  compute eddy
%  shapes is saved as [path_out,'warnings_shapes']:
%
%  - warn_shapes(t).land(n) : 1 if land was found, and the area where the
%                             eddy shape was computed was resized;
%  - warn_shapes(t).no_curve(n): 1 if no closed contour of PSI was found
%                             around the eddy center. Eddy is assumed
%                             circular with radius "a-1";
%  - warn_shapes(t).large_curve(n): 1 if no closed contour around the
%                             center with increasing velocity across. Shape
%                             is defined by just the larges closed contour;
%  - warn_shapes(t).fac(n): number of times the area where the shape is
%                             computed was enlarged;
%  (n is the total number of eddy for a given day; t the total number of
%  days)
%
%  The function returns also a log file [path_out,'log_eddy_shapes.txt']
%  which contains additional information on the process to compute eddy
%  shapes.
%  (NOTE: if the file already exist, the new log file will be append to it)
%
%  'eddy_dim' is used to compute eddy shapes. Check the documentation in
%  eddy_dim.m for further details.
%
%-------------------------
%   Ver. 2.0 Jan.2012
%   Ver. 1.3 Apr.2011
%   Ver. 1.2 May.2010
%   Ver. 1.1 Dec.2009
%   Authors: Francesco Nencioli, francesco.nencioli@univ-amu.fr
%            Charles Dong, cdong@atmos.ucla.edu
%-------------------------
%
% Copyright (C) 2009-2012 Francesco Nencioli and Charles Dong
%
% This file is part of the Vector-Geometry Eddy Detection Algorithm.
%
% The Vector-Geometry Eddy Detection Algorithm is free software:
% you can redistribute it and/or modify it under the terms of the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
%
% The Vector-Geometry Eddy Detection Algorithm is distributed in the
% hope that it will be useful, but WITHOUT ANY WARRANTY; without even
% the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE. See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the Vector-Geometry Eddy Detection Algorithm.
% If not, see <http://www.gnu.org/licenses/>.
%
%=========================

% begin the log file
diary([path_out,'log_eddy_shapes.txt']);
%----------------------------------------------
% load eddy centers
load([path_out,'eddy_centers']);
% load time variable
time=nc_varget(nc_u,'day');
% load coordinates
lon=nc_varget(nc_dim,'lon');
lat=nc_varget(nc_dim,'lat');

LON=lon; % no expantion lon
LAT=lat; % no expantion lat
[im jm]=size(LAT);
% load mask
mask=nc_varget(nc_dim,'mask');
MASK=mask;
%----------------------------------------------
% Added option for periodic East-West boundaries
% (first and last column of the domain are identical)
global periodic
if periodic==1
    % check if minimum and maximum longitude are constant with latitude
    if ~isempty(find(diff(min(lon,[],2))~=0)) | ...
            ~isempty(find(diff(max(lon,[],2))~=0))
        error('Problem with your domain: Lonmax and Lonmin vary with Lat!')
    end
    mlon=min(lon(1,:));
    Mlon=max(lon(1,:));
    % shifted longitudinal positions
    Dlon=lon+(Mlon-mlon);
    dlon=lon-(Mlon-mlon);
    % lat and lon are expanded by adding the shifted positions to the
    % beginning and to the end of the domain;
    % (they are now 3x their original size: in case you wanna reduce this
    % make sure to change also u and v expansion in eddy_dim.m)
    lon=[dlon(:,1:end-1),lon,Dlon(:,2:end)];
    lat=[lat(:,1:end-1),lat,lat(:,2:end)];
end
% coordinates extremes
Lonmin=min(min(lon));
Lonmax=max(max(lon));
Latmin=min(min(lat));
Latmax=max(max(lat));

% Compute eddy shape --------------------------
% preallocate shape and warning array
shapes = repmat(struct('lonlat',{},'eddy_vor',[],'eddy_mvor',[],'eddy_eke',[],...
    'area',[],'radius',[],'EI',[],'Q',[]),1,length(centers));
warn_shapes = repmat(struct('land',{},'no_curve',{}, ...
    'large_curve',{},'fac',{}),1,length(centers));
% loop through all days of the time series
for i=1:length(centers)
    disp(['Day ',num2str(time(i)),' %-------------'])
    iday=find(time==centers(i).day);
    % load 2D velocity field -----------------------------------
    U=nc_varget(nc_u,'ssu',[iday-1 0 0],[1 Inf Inf]);
    V=nc_varget(nc_v,'ssv',[iday-1 0 0],[1 Inf Inf]);
    %-----------------------------------------------------------
    % mask velocity data
    U(MASK==0)=NaN;
    V(MASK==0)=NaN;
    % loop through all centers detected for a given day
    for ii=1:length(centers(i).type)
        disp([' === Eddy ',num2str(ii),' ==='])
        % factor to increase the area where PSI is computed, if eddy
        % dimesions are too big!!!
        fac=1; % initially set to no increase
        % lonlat is the computed shape;
        % the others are flags output by eddy_dim, and saved in
        %  'warnings_shapes';
        % box is a flag that indicates if eddy shape is close to the area
        %  boundaries (1.5 grid points; limit set in eddy_dim);
        [lonlat warn land box large]=eddy_dim(nc_u,nc_v, ...
            centers(i).day,centers(i).i(ii),centers(i).j(ii),...
            lon,lat,Lonmin,Lonmax,Latmin,Latmax, ...
            mask,fac,rad,a);
        % temporary save eddy_shape
        tmp_lonlat=lonlat;
        % log file messages -------------------------
        if large==1 || warn==1 || land==1 || box==1
            if land==1
                disp('   Found land!!!')
            end
            if warn==1
                disp(['   Warning: No streamlines ', ...
                    'closed around the center!!!!!!'])
            end
            if large==1
                disp('   Largest closed curve')
            end
            if large==1 || warn==1 || land==1
                disp(' ')
            end
        end
        %----------------------------------------------------------
        %----------------------------------------------------------
        % area where psi is computed is enalrged while the eddy shape
        % is close to the area boundaries (1.5 grid points; limit set in
        % eddy_dim), and no land is found within the area
        while land==0 && box==1 && large==0 && warn==0
            fac=fac+1; % this determine larger area
            disp(['   Big eddy: going to fac = ',num2str(fac)])
            % compute eddy shape in the larger area
            [lonlat warn land box large]=eddy_dim(nc_u,nc_v, ...
                centers(i).day,centers(i).i(ii),centers(i).j(ii),...
                lon,lat,Lonmin,Lonmax,Latmin,Latmax, ...
                mask,fac,rad,a);
            
            % flags -----------------------------------
            if land==1
                disp('     Found land!!!')
            end
            % if no closed curve in the larger area then final eddy shape
            % is the one computed in the smaller area
            if warn==1 || large==1
                disp(['     No closed or largest curve at fac ', ...
                    num2str(fac),' back to curve at fac ', ...
                    num2str(fac-1),'!!!'])
                lonlat=tmp_lonlat;
                % if eddy shape still close to the area boundaries then
                % temporary save lonlat
            elseif land==0
                tmp_lonlat=lonlat;
            end
            if warn==1 || large==1 || box==0 || land==1
                disp(' ')
            end
            %-------------------------------------------
        end
        %----------------------------------------------------------
        % This is to fix the flags saved in case going to larger area
        % resulted in no closed conotur of PSI around the center
        if (large==1 || warn==1) && fac>1
            fac=fac-1;
            large=0;
            warn=0;
        end
        % in case of periodic domain the eddy dimensions are shifted back
        % to fit within the initial domain boundaries
        if periodic==1
            lonlat(1,lonlat(1,:)>Mlon)=lonlat(1,lonlat(1,:)>Mlon)-(Mlon-mlon);
            lonlat(1,lonlat(1,:)<mlon)=lonlat(1,lonlat(1,:)<mlon)+(Mlon-mlon);
        end
        
        %%%%%%%%%%%%%% written by Hu %%%%%%%%%%%%%%%%%%%
        
        % boundary indx of rectangle containing the eddy
        lon01=max(fix(min(lonlat(3,:))-1),1);
        lon02=min(fix(max(lonlat(3,:))+1),jm);
        lat01=max(fix(min(lonlat(4,:))-1),1);
        lat02=min(fix(max(lonlat(4,:))+1),im);
        indx=zeros(lat02-lat01+1,lon02-lon01+1);
        q=1;
        for il=lon01:lon02
            p=1;
            for jl=lat01:lat02
                if inside(il,jl,lonlat(3,:),lonlat(4,:))==1 ...
                        % decide whether the point is inside the eddy.  refer
                        % to inside.m for detail
                    indx(p,q)=1; % index of the point in the eddy
                end
                p=p+1;
            end
            q=q+1;
        end
        
        %--- compute area & radius of eddy---%
        area=cal_area(lonlat(1,:),lonlat(2,:));
        area=area*(111*111*cos(centers(i).lat(ii)*pi/180));  % km^2
        radius=sqrt(area/pi);   % km
        
        %--- compute vorticity ---%
        %compute vorticity of rectangle
        vor_rect=vor(LAT(lat01:lat02,1),LON(1,lon01:lon02),...
            U(lat01:lat02,lon01:lon02),V(lat01:lat02,...
            lon01:lon02));  %need vor.m! copy it to this path
        eddy_vor=nanmax(nanmax(abs(vor_rect(indx==1))))...
            *sign(centers(i).type(ii));  %mean vorticity of the eddy
        vortemp=vor_rect(indx==1);
        eddy_mvor=nanmean(vortemp(:));
        
        %--- compute eke ---%
        u0=U(lat01:lat02,lon01:lon02);
        v0=V(lat01:lat02,lon01:lon02);
        eddy_eke=nanmean((u0(indx==1).^2+v0(indx==1).^2)/2); %mean eke of the eddy
        clear u0 v0;
        
        %--- compute EI ---%
        eddy_eke_cm=eddy_eke*1e4;
        EI=eddy_eke_cm/area;  % cm^2/s^2km^2
        
        %--- compute Q ---%
        ux=zeros(size(indx));
        uy=zeros(size(indx));
        vx=zeros(size(indx));
        vy=zeros(size(indx));
        for jl=lat01:lat02
            dx=1.11e5*cos(LAT(jl,1)*pi/180)*.25;
            ux(jl-lat01+1,:)=gradient(U(jl,lon01:lon02)/100,dx);
            vx(jl-lat01+1,:)=gradient(V(jl,lon01:lon02)/100,dx);
        end
        yy=(LAT(lat01:lat02,1)-LAT(lat01,1))*1.11e5;
        [~,uy]=gradient(U(lat01:lat02,lon01:lon02),1,yy);
        [~,vy]=gradient(V(lat01:lat02,lon01:lon02),1,yy);
        Q=-ux.^2-vx.*uy;
        Q=nanmean(Q(indx==1));  % unit: s^-2
        
        %%%%%%%%%%%%%% written by Hu %%%%%%%%%%%%%%%%%%%
        
        % save dim in a struct array
        shapes(i).lonlat(ii)={lonlat};
        shapes(i).eddy_vor(ii)=eddy_vor;
        shapes(i).eddy_mvor(ii)=eddy_mvor;
        shapes(i).eddy_eke(ii)=eddy_eke;
        shapes(i).area(ii)=area;
        shapes(i).radius(ii)=radius;
        shapes(i).EI(ii)=EI;
        shapes(i).Q(ii)=Q;
        % warnings from shape computation
        warn_shapes(i).land(ii)=land;
        warn_shapes(i).no_curve(ii)=warn;
        warn_shapes(i).large_curve(ii)=large;
        warn_shapes(i).fac(ii)=fac;
    end
end
% output eddy_shapes and shapes warnings
save([path_out,'eddy_shapes'],'shapes')
save([path_out,'warnings_shapes'],'warn_shapes')
% close log file
diary off
