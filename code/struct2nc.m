function struct2nc(x,ncfile,ncfiletype,deflate_lev)
% STRUCT2NC writes all float,double and character vars to netcdf
% Usage: struct2nc(x,ncfile,[ncfiletype],[deflate_lev])
% x = structure
% ncfile = name of netcdf output file (e.g. 'test.nc')
% ncfiletype = netcdf file type (e.g. 'classic','netcdf4_classic')
% deflate_lev = deflate level (0-9, 0 is none)
%
% This function writes all 'double','single' and 'char' variables
% to NetCDF using the native Matlab NetCDF interface.  It skips all
% other classes in the struct (e.g. structs, cell arrays, etc).  It
% also only handles scalar, 1D, 2D, and 3D arrays currently, although
% this could easily be extended.

if nargin==2,
    ncfiletype='classic';
    deflate_lev=0;
elseif nargin==3;
    switch ncfiletype
        case {'netcdf4','netcdf4_classic'}
            deflate_lev=6;
        otherwise
            deflate_lev=0;
    end
end
s      = fieldnames(x);
layers = length(x);
k=0;
% create variables first, but don't write data
for l=1:layers
    for i=1:length(s)
        vname=char(s(i));
        if l>1
           uname=[vname,num2str(l-1)];
        else
           uname=vname;
        end
            var   =x(l).(vname);
            vtype = class(var);
            vshape= size(var);
            ndims = length(vshape);
            vlen  = length(var(:));
    switch vtype;
        case {'double','single','int32','uint8'},
            if vlen==1,
                nccreate(ncfile,uname,...
                    'Datatype',vtype,'format',ncfiletype);
                k=k+1;
                vnames{k}=vname;
                unames{k}=uname;
            else
                if min(vshape)==1,
                    nccreate(ncfile,uname,...
                        'Datatype',vtype,...
                        'DeflateLevel',deflate_lev,...
                        'Dimensions',{[uname '1'] vlen},...
                        'format',ncfiletype);
                    k=k+1;
                    vnames{k}=vname;
                    unames{k}=uname;                    
                elseif ndims==2,
                    nccreate(ncfile,uname,...
                        'Datatype',vtype,...
                        'DeflateLevel',deflate_lev,...
                        'Dimensions',{[uname '1'] vshape(1) [uname '2'] vshape(2)},...
                        'format',ncfiletype);
                    k=k+1;
                    vnames{k}=vname;
                    unames{k}=uname;                    
                elseif ndims==3,
                    nccreate(ncfile,uname,...
                        'Datatype',vtype,...
                        'DeflateLevel',deflate_lev,...
                        'Dimensions',...
                        {[uname '1'] vshape(1) [uname '2'] vshape(2) [uname '3'] vshape(3)},...
                        'format',ncfiletype);
                    k=k+1;
                    vnames{k}=vname;
                    unames{k}=uname;                    
                else,
                    disp('Skipping variable with more than 3 dimensions');
                end
            end
       case {'char'},
         if min(vshape)==1,
            nccreate(ncfile,uname,...
                'Datatype',vtype,...
                'Dimensions',{[uname '1'] vlen},.....
                'format',ncfiletype);
            k=k+1;
            vnames{k}=vname;
            unames{k}=uname;
         elseif ndims==2,
             nccreate(ncfile,uname,...
                      'Datatype',vtype,...
                      'Dimensions',{[uname '1'] vshape(1) [uname '2'] vshape(2)},...
                      'format',ncfiletype);
             k=k+1;
             vnames{k}=vname;
             unames{k}=uname;                    
         end
        otherwise,
            disp(['skipping ' vname])
    end
end
end
%write all the data at the end
for l=1:layers
for i=1:length(vnames)
    ncwrite(ncfile,unames{i},x(l).(vnames{i}));
end
end