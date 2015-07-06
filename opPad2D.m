classdef opPad2D < opSpot
    % Operator for zero padding a signal
    
    
    properties
        
      pad_value, dim_in, dim_out
    end
    
    
    methods
         function op = opPad2D(dim_in, dim_out, pad_value)    
            op = op@opSpot('opPad', prod(dim_out), prod(dim_in));
            op.cflag = 0;
            op.linear = 1;
            op.children = [];

            % Store the natural dimensions         
            op.pad_value = pad_value;
            op.dim_in = dim_in;
            op.dim_out = dim_out;
         end %constructor
        
    end
    
     methods(Access = protected)
        function y = multiply(op,x,mode)
       
            if mode==1
                x = reshape(x, op.dim_in);
                pad = op.dim_out - op.dim_in;
                
                y = padarray(x, pad, op.pad_value, 'post');
            else % adjoint
                x = reshape(x, op.dim_out);
                pad = op.dim_out - op.dim_in;
                
                y = x(1:end-pad(1), 1:end-pad(2));
            end
            
            y = y(:);
                
        end
     end
end      
                
