function mod_eddy_centers(nc_u,nc_v,nc_dim,a,b,path_out)
% mod_eddy_centers.m
%  mode_eddy_centers(nc_u,nc_v,nc_dim,a,b,path_out) detect the eddy centers  
%  present in the domain, defined by nc_dim, for each time step of the time 
%  series of the 2-D velocity field defined by nc_u and nc_v, using the 
%  parameters a and b.
%
%  For a description of the input parameters see param_eddy_tracking.m.
%
%  Eddy centers are saved in the structure array [path_out,'eddy_centers']:
%
%  - centers(t).day : day when the eddy was detected
%  - centers(t).type(n) : eddy type (1 => cyclonic; -1 => anticyclonic)
%  - centers(t).lat(n) : eddy center latitude
%  - centers(t).lon(n) : eddy center longitude
%  - centers(t).i(n) : eddy center column index
%  - centers(t).j(n) : eddy center row index
%  
%  (t is the time index; n is the number of eddies detected at t)
%
%  'uv_search' is used to determine the points in the domain that satisfy 
%  all 4 constraints. Check the documentation in uv_search.m for further 
%  details.
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

%------- load time, coordinates and mask -------
Lat=nc_varget(nc_dim,'lat');
Lon=nc_varget(nc_dim,'lon');
mask=nc_varget(nc_dim,'mask');

time=nc_varget(nc_u,'day');
%---------------------------------------------
% preallocate centers array
centers = repmat(struct('day',{},'type',{},'lat',{},'lon',{}, ...
    'j',{},'i',{}),1, 10);

% cycle through time steps
for i=1:length(time)
    disp(['Searching day ',num2str(time(i)),' ... '])
    % eddy centers for a given day
    
    [eddy_fld_uv eddy_fld_c eddy_fld]=uv_search(nc_u,nc_v,Lon,Lat,mask,i,a,b);
    % save eddy positions in struct array
    k=1;
    if ~isempty(eddy_fld)
        for ii=1:length(eddy_fld(:,1))
            if eddy_fld(ii,1)>5 && eddy_fld(ii,1)<20 && eddy_fld(ii,2)>108 && eddy_fld(ii,2)<120
                centers(i).type(k)=eddy_fld(ii,3);
                centers(i).lat(k)=eddy_fld(ii,1);
                centers(i).lon(k)=eddy_fld(ii,2);
                [centers(i).j(k) centers(i).i(k)]=find(Lon==eddy_fld(ii,2) ...
                    & Lat==eddy_fld(ii,1));
                centers(i).day=time(i);
                k=k+1;
            end
        end
    end
        uv_centers(i).day=time(i);
    if ~isempty(eddy_fld_uv)
        for ii=1:length(eddy_fld_uv(:,1))
            uv_centers(i).lat(ii)=eddy_fld_uv(ii,1);
            uv_centers(i).lon(ii)=eddy_fld_uv(ii,2);
        end
    end
    
        c_centers(i).day=time(i);
    if ~isempty(eddy_fld_c)
        for ii=1:length(eddy_fld_c(:,1))
            c_centers(i).lat(ii)=eddy_fld_c(ii,1);
            c_centers(i).lon(ii)=eddy_fld_c(ii,2);
        end
    end
    
    
end

save([path_out,'eddy_centers'],'centers')
save([path_out,'c_centers'],'c_centers')
save([path_out,'uv_centers'],'uv_centers')
