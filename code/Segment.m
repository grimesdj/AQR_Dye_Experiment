function [startINDs, endINDs] = Segment(time,cuttoff);
% 
% Usage: [startINDs, endINDs] = Segment(time,cuttoff);
%
% This function will log the start and end indices for nearly
% monotonic sections amidst large steps. time is the vector of time
% stamps to be split up and cuttof is the max number of timesteps to
% allow before generating new segment. 

%Clip the timeseries into monotonic segments
[r,c] = size(time);

dt = diff(time);% in seconds

nominal_sample_rate = nanmedian(dt);

seg_ind = find(dt >= ones(size(dt))*cuttoff*nominal_sample_rate);

startINDs = [0; seg_ind]+1;
endINDs = [seg_ind; r];