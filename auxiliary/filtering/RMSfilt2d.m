function [ RMStrace ] = RMSfilt2d( xy_trace, wSize )
%RMSfilt2d: Returns rms trace for 2dim raw data (eg fitted positions)

% input and output size
if size(xy_trace,2)>2 && size(xy_trace,1)==2
    xy_trace = xy_trace';
    resh = true;
else
    resh = false;
end

% set up
delta=floor(wSize/2);
X=(xy_trace(:,1));
Y=(xy_trace(:,2));
RMStrace=zeros(length(xy_trace),1);

% main part
for i=delta+1:length(xy_trace)-delta
    x_mean = mean(X(i-delta:i+delta));
    y_mean = mean(Y(i-delta:i+delta));
    RMStrace(i)=sqrt(mean((X(i-delta:i+delta)-x_mean).^2+(Y(i-delta:i+delta)-y_mean).^2));
end
% first delta frames
RMStrace(1:delta) = RMStrace(delta+1).*ones(delta,1);
% last delta frames
RMStrace(end-delta+1:end) = RMStrace(end-delta).*ones(delta,1);

% adjust output size
if resh
    RMStrace = RMStrace';
end
end

