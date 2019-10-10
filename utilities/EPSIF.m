function output = EPSIF(input,threshold,w_size,iter_max,mode,sigma_s,sigma_r,guide)

%% init
if ~exist('threshold','var')
    threshold = 0.45;
    mode = 'simple';
    iter_max = 3;
    w_size = 9;
elseif ~exist('w_size','var')
    threshold = double(threshold);
    mode = 'simple';
    iter_max = 3;
    w_size = 9;
elseif ~exist('iter_max','var')
    threshold = double(threshold);
    mode = 'simple';
    iter_max = 3;
    w_size = ceil(w_size/2)*2-1;
elseif ~exist('mode','var')
    threshold = double(threshold);
    mode = 'simple';
    iter_max = ceil(iter_max);
    w_size = ceil(w_size/2)*2-1;
end

if ~isa(input, 'double')
    input = im2double(input);
end

%% process

if strcmp( mode , 'simple')   % fast mode
    
    if w_size == 15
        output = ep2d_simple_15(input,threshold,iter_max);
    elseif w_size == 9
        output = ep2d_simple_9(input,threshold,iter_max);
    elseif w_size == 31
        output = ep2d_simple_31(input,threshold,iter_max);
    else
        output = ep2d_simple(input,w_size,threshold,iter_max);
    end
    
    
elseif   strcmp( mode , 'full') % full mode
    fprintf('Processed in full mode.\n')
    output = ep2d_full(input,w_size,sigma_s,sigma_r,threshold,iter_max);
    
elseif strcmp( mode , 'guided') %guided mode
    if ~exist('guide','var')
        guide = input;
    end
    
    if abs(size(input,1)-size(guide,1)) + abs(size(input,2)-size(guide,2)) == 0
        output = ep2d_guided(input,w_size,threshold,iter_max,guide);
    else
        output = zeros(size(input));
    end
    
else
    
    output = zeros(size(input));
    
end

end