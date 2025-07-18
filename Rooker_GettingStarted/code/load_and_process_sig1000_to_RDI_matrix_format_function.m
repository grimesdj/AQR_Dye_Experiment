function A = load_and_process_sig1000_to_RDI_matrix_format_function(Config, rootDIR, fRoot, L0dir, filePrefix, hab, echo_mode, deployTime, recoverTime, HeadingOffset)

files = dir([rootDIR,filesep,fRoot,'*.mat']);
Nf    = length(files);
%
%
% average the start/end pressures
ATM_Pressure = mean(Config.ATM_Pressure);
%
% get bin sizing etc for currents
blnk = Config.Burst_BlankingDistance;
binw = Config.Burst_CellSize;
Nc   = double(Config.Burst_NCells);
mab  = hab+blnk+binw.*[1:Nc]';
%
%  get bin sizing etc for echo
if echo_mode
blnkE= Config.EchoSounder_BlankingDistance;
binwE= Config.EchoSounder_CellSize;
NcE  = double(Config.EchoSounder_NCells);
mabE = hab+blnkE+binwE.*[1:NcE]';
end
%
%save the config info
outFile = [L0dir, filePrefix, 'config.mat'];
if ~exist(L0dir,'dir')
    eval(['!mkdir ',L0dir])
end
save(outFile,'Config')
% outFile = [L0dir, filePrefix, 'config.nc'];
% if exist(outFile,'file')
%     eval(['!rm ', outFile])
% end
%struct2nc(Config,outFile,'NETCDF4')


%
% limit archive file to 24hr length
maxDuration = 24*3600;
%
% get sample rate
fs = double(Config.Burst_SamplingRate);
%
% maximum number of samples per file
Nmax = fs*maxDuration;
% number of measurements in current file
N    = 0;
%
% flag to load next file
loadFlag = 1;
ii = 0;
outNf = 0;
fprintf('\n \n')



while ii<=Nf
    %
    if loadFlag
        ii  = ii+1;
        loadFlag=0;
        fileName = [fRoot,num2str(ii),'.mat'];
        fin = [files(ii).folder,'/',fileName];
        fprintf(['pre-processing:       %s \n'], fileName)
        load(fin)
        % use 5th beam time; the other beams are offset by dt = 1/(2*fs)
        t    = Data.IBurst_Time;
        nt   = length(t);
        % check if instrument is deployed
        is   = find(t>=deployTime,1,'first');
        if isempty(is)
            disp(['--> instrument is not yet deployed, skipping this file.'])
            loadFlag=1;
            continue
        end
        % check if instrument was pulled
        iStop= find(t>=recoverTime,1,'first');
        if ~isempty(iStop)
            disp(['trimming data at recovery time: nt = ',num2str(nt),'--> ',num2str(iStop-1)]);
            nt = iStop-1;
        end
    end
    %
    if N==0
        % initialize output structure
        out  = struct('Time',[],'Heading',[],'Pitch',[],'Roll',[],'Temperature',[],'Pressure',[],'Battery',[],...
                      'Velocity_Beam',[],...,'VelBeam2',[],'VelBeam3',[],'VelBeam4',[],'VelBeam5',[],...
                      'Amplitude_Beam',[],...'AmpBeam2',[],'AmpBeam3',[],'AmpBeam4',[],'AmpBeam5',[],...
                      'Correlation_Beam',[],...'CorBeam2',[],'CorBeam3',[],'CorBeam4',[],'CorBeam5',[],...,
                      'Velocity_East',[],'Velocity_North',[],'Velocity_Up',[],'Velocity_Error',[],...
                      'bin_mab',mab);
        if echo_mode
            out.bin_mab_Echo = mabE;
            out.Echo1=[];
            out.Echo2=[];
        end
    end
    %
    if N+nt-(is-1)<=Nmax
        ie   = nt;
        loadFlag = 1;
    else
        ie   = Nmax-N-(is-1);
        loadFlag = 0;
    end
    %
    out.Time       (1,N+1:N+(ie-(is-1)))     = Data.IBurst_Time(is:ie,1);
    out.Heading    (1,N+1:N+(ie-(is-1)))     = Data.Burst_Heading(is:ie,1);
    out.Pitch      (1,N+1:N+(ie-(is-1)))     = Data.Burst_Pitch(is:ie,1);
    out.Roll       (1,N+1:N+(ie-(is-1)))     = Data.Burst_Roll(is:ie,1);
    out.Temperature(1:2,N+1:N+(ie-(is-1)))   = [Data.IBurst_Temperature(is:ie,1)'; Data.Burst_Temperature(is:ie,1)'];
    out.Pressure   (1:2,N+1:N+(ie-(is-1)))   = [Data.IBurst_Pressure(is:ie,1)';    Data.Burst_Pressure(is:ie,1)'] - ATM_Pressure;
    out.Battery    (1,N+1:N+(ie-(is-1)))     = Data.IBurst_Battery(is:ie,1);
    out.Velocity_Beam    (1:Nc,N+1:N+(ie-(is-1)),1:5) = permute(cat(2,Data.Burst_Velocity_Beam(is:ie,:,1:Nc) ,Data.IBurst_Velocity_Beam(is:ie,1,1:Nc)), [3 1 2]) ;
    out.Amplitude_Beam   (1:Nc,N+1:N+(ie-(is-1)),1:5) = permute(cat(2,Data.Burst_Amplitude_Beam(is:ie,:,1:Nc),Data.IBurst_Amplitude_Beam(is:ie,:,1:Nc)), [3 1 2]);
    out.Correlation_Beam (1:Nc,N+1:N+(ie-(is-1)),1:5) = permute(cat(2,Data.Burst_Correlation_Beam(is:ie,:,1:Nc),Data.IBurst_Correlation_Beam(is:ie,:,1:Nc)), [3 1 2]);
% $$$     out.VelBeam1   (N+1:N+(ie-(is-1)),1:4,1:Nc)  = Data.Burst_Velocity_Beam(is:ie,1,1:Nc);
% $$$     out.VelBeam2   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Velocity_Beam(is:ie,2,1:Nc);
% $$$     out.VelBeam3   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Velocity_Beam(is:ie,3,1:Nc);
% $$$     out.VelBeam4   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Velocity_Beam(is:ie,4,1:Nc);    
% $$$     out.VelBeam5   (N+1:N+(ie-(is-1)),1:Nc)  = Data.IBurst_Velocity_Beam(is:ie,1,1:Nc);
% $$$     out.AmpBeam1   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Amplitude_Beam(is:ie,1,1:Nc);
% $$$     out.AmpBeam2   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Amplitude_Beam(is:ie,2,1:Nc);
% $$$     out.AmpBeam3   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Amplitude_Beam(is:ie,3,1:Nc);
% $$$     out.AmpBeam4   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Amplitude_Beam(is:ie,4,1:Nc);    
% $$$     out.AmpBeam5   (N+1:N+(ie-(is-1)),1:Nc)  = Data.IBurst_Amplitude_Beam(is:ie,1,1:Nc);
% $$$     out.CorBeam1   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Correlation_Beam(is:ie,1,1:Nc);
% $$$     out.CorBeam2   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Correlation_Beam(is:ie,2,1:Nc);
% $$$     out.CorBeam3   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Correlation_Beam(is:ie,3,1:Nc);
% $$$     out.CorBeam4   (N+1:N+(ie-(is-1)),1:Nc)  = Data.Burst_Correlation_Beam(is:ie,4,1:Nc);    
% $$$     out.CorBeam5   (N+1:N+(ie-(is-1)),1:Nc)  = Data.IBurst_Correlation_Beam(is:ie,1,1:Nc);
    out.Velocity_East    (1:Nc,N+1:N+(ie-(is-1)))  = permute(Data.Burst_Velocity_ENU(is:ie,1,1:Nc),[3,1,2]);
    out.Velocity_North   (1:Nc,N+1:N+(ie-(is-1)))  = permute(Data.Burst_Velocity_ENU(is:ie,2,1:Nc),[3,1,2]);
    out.Velocity_Up      (1:Nc,N+1:N+(ie-(is-1)))  = permute(Data.Burst_Velocity_ENU(is:ie,3,1:Nc),[3,1,2]);
    out.Velocity_Error   (1:Nc,N+1:N+(ie-(is-1)))  = permute(Data.Burst_Velocity_ENU(is:ie,4,1:Nc),[3,1,2]);
    out.AST              (1,N+1:N+(ie-(is-1)))     = Data.Burst_AltimeterDistanceAST(is:ie,1);
    out.AST_Offset       (1,N+1:N+(ie-(is-1)))     = Data.Burst_AltimeterTimeOffsetAST(is:ie,1);

    if echo_mode
    out.Echo1            (1:NcE,N+1:N+(ie-(is-1))) = Data.Echo1Bin1_1000kHz_Echo(is:ie,1:NcE)';
    out.Echo2            (1:NcE,N+1:N+(ie-(is-1))) = Data.Echo2Bin1_1000kHz_Echo(is:ie,1:NcE)';
    end
    %
    N = N+ie-(is-1);
    is = ie+1;
    %
    if N==Nmax || (ii==Nf & loadFlag)
        % process data
        minC  = min( out.Correlation_Beam, [], 3);%min( cat(3,out.CorBeam1, out.CorBeam2, out.CorBeam3, out.CorBeam4),[], 3);
        minA  = min( out.Amplitude_Beam  , [], 3);  %cat(3,out.CorBeam1, out.CorBeam2, out.CorBeam3, out.CorBeam4),[], 3);
        maxRNG= out.Pressure(2,:)*cosd(25)-binw;
        qcFlag= (minC>50 & minA>30 & maxRNG>mab);
        %
        out.Correlation_Minimum = minC;
        out.Amplitude_Minimum = minA;
        out.maxRNG = maxRNG;
        out.qcFlag = qcFlag;
        %
        disp('not applying heading offset correction, velocities are relative to magnetic north/south')
% $$$         % now rotate to xyz and ENU
% $$$         [out,Config2,beam2xyz] = beam2earth_sig1000_DG_FRFcoords(out,Config,'');        
        out.HeadingOffset = HeadingOffset;
% $$$         out.beam2xyz = beam2xyz;
        A = out;
        % save output
         outNf = outNf+1;
         outFileName = sprintf([filePrefix,'%03d.mat'],outNf);
         fout  = [L0dir,outFileName];
         fprintf(['saving output file:   %s \n'],outFileName)
         save(fout,'-struct','out')
%         outFileName = sprintf([filePrefix,'%03d.nc'],outNf);
 %        fout = [L0dir,outFileName];
  %       if exist(fout,'file')
   %         eval(['!rm ', outFileName])
    %     end
%         struct2nc(out,outFileName,'NETCDF4')
%         %
         N = 0;
%         %
         if (ii==Nf & loadFlag )
         fprintf('done! \n')
         break
         end
    end
    %
    %
end
end