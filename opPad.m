classdef opPad2D < opSpot
    % Operator for zero padding a signal
    
    
    properties
        
       direction, pad_value, dim_in, dim_out
    end
    
    
    methods
         function op = opPad2D(dim_in, dim_out, direction, pad_type)    
            op = op@opSpot('opPad', prod(dim_out), prod(dim_in));
            op.cflag = 1;
            op.linear = 1;
            op.children = [];

            % Store the natural dimensions         
            op.direction = direction;
            op.pad_type = pad_type;
            op.dim_in = dim_in;
            op.dim_out = dim_out;
         end %constructor
        
    end
    
     methods(Access = protected)
        function y = multiply(op,x,mode)
       
            if mode==1
                x = reshape(x, dim_in);
                pad = dim_out - dim_in;
                
                y = padarray(x, pad, op.pad_value, op.direction);
            else % adjoint
                x = reshape(x, dim_out);
                pad = dim_out-dim_in;
                
                
            end
                
                
