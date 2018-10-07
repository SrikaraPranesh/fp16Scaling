function z = horzcat(varargin)
    z = fp16(varargin{1});
    for k = 2:nargin
        x = fp16(varargin{k});
        z.u = [z.u x.u];
    end
end