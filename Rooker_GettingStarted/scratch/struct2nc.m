function struct2nc(data, name, mode)
%STRUCT2NC Save a MATLAB struct to a NetCDF file.
%
%   struct2nc(data, name, mode)
%   - data: struct with numeric fields
%   - name: filename for the .nc file (e.g., 'output.nc')
%   - mode: NetCDF creation mode (e.g., 'CLOBBER' or bitwise bitor of constants)

ncid = netcdf.create(name, mode);
fields = fieldnames(data);
varInfo = struct;

for i = 1:numel(fields)
    fname = fields{i};
    val = data.(fname);

    if isnumeric(val) || islogical(val)
        if islogical(val)
            val = uint8(val); % Convert logical to numeric
        end

        dims = size(val);
        dimIDs = [];

        for d = 1:length(dims)
            dimName = sprintf('%s_dim%d', fname, d);
            dimIDs(d) = netcdf.defDim(ncid, dimName, dims(d));
        end

        varid = netcdf.defVar(ncid, fname, 'double', dimIDs);
        varInfo(i).id = varid;
        varInfo(i).val = val;

    elseif ischar(val) || isstring(val)
        % Save strings as global attributes
        netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), fname, char(val));

    elseif isscalar(val)
        % Save scalar numerics as global attributes
        if isnumeric(val)
            netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), fname, val);
        end
    else
        warning('Skipping field "%s" (unsupported type or nested struct)', fname);
    end
end

netcdf.endDef(ncid);

for i = 1:numel(fields)
    if isfield(varInfo, 'id')
        try
            netcdf.putVar(ncid, varInfo(i).id, varInfo(i).val);
        catch
            warning('Could not write field "%s"', fields{i});
        end
    end
end

netcdf.close(ncid);
fprintf('✅ NetCDF file "%s" written.\n', name);
