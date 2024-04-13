function ramp_filter = ramp_filter(ramp_length,signal_lenght)
% Ramp-filter based on Hanning Window Filter
hann_filter = hann(2*ramp_length);
ramp_filter = ones([signal_lenght,1]);
first_ramp = hann_filter(1:ramp_length);
last_ramp = hann_filter(ramp_length+1:end);
ramp_filter(1:ramp_length) = first_ramp;
ramp_filter(signal_lenght-ramp_length+1:end) = last_ramp;

end