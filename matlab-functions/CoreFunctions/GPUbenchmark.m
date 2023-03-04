function GPUbenchmark(deviceID)
% GPUbenchmark(1) runs on the benmark on GPU 1
% GPUbenchmark(1) runs on the benmark on GPU 2

    if isempty(deviceID)
        deviceID = 1;
    end
    disp(gpuDevice(deviceID));
    
    a = complex(randn(4*4096,4096),randn(4*4096,4096));   % Data input  % was 4096 100
    b = randn(24,1);                                % Filter input
    c = fastConvolution_v2(a,b);                    % Calculate output
    ctime = timeit(@()fastConvolution_v2(a,b));     % Measure CPU time
    disp(['Execution time on CPU = ',num2str(ctime)]);

    ga = gpuArray(a);                               % Move data to GPU
    gb = gpuArray(b);                               % Move filter to GPU
    gc = fastConvolution_v2(ga, gb);                % Calculate on GPU
    gtime = gputimeit(@()fastConvolution_v2(ga,gb));% Measure GPU time
    gerr = max(max(abs(gather(gc)-c)));             % Calculate error
    disp(['Execution time on GPU = ',num2str(gtime)]);
    disp(['Maximum absolute error = ',num2str(gerr)]);
    
end

function y = fastConvolution_v2(data,filter)
    m = size(data,1);
    % Zero-pad filter to the length of data, and transform
    filter_f = fft(filter,m);

    % Transform each column of the input
    af = fft(data);

    % Multiply each column by filter and compute inverse transform
    y = ifft(bsxfun(@times,af,filter_f));
end