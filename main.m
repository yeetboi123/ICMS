function [y,transform_data,varargout] = main(x1,x2,varargin)
transform_data = @nestedfun1;

    disp('Hello world!')
    DATA_DIR = uigetdir(pwd);
    data=load(fullfile(DATA_DIR,'data.mat')).data;
    
    y = x1+x2;
    
    preprocessed_data = preprocessData('lpf');
    transformed_data = nestedfun1('preprocess',3,2);

    switch nargout 
        case 3
            varargout{1} = data;
        case 4
            varargout = {data transformed_data};
        otherwise
            disp('Error!')
    end
    
    switch nargin 
        case 2
            disp('2 inputs')
        case 3
            fprintf('3rd input: %s\n', varargin{1})
        otherwise
            disp('All optional arguments:')
            fprintf(1, '%s\n', varargin{:})
            fprintf(1, '\n')
    end
    
    function preprocessed_data = preprocessData(method)
        switch method
            case 'lpf'
                preprocessed_data = data * 0.5;
            case 'hpf'
                preprocessed_data = data * 5;
            otherwise
                disp('Error!')
        end        
    end
    
    function transformed_data = nestedfun1(preprocess_flag,idx,exp_value)
        if strcmp('preprocess',preprocess_flag) == 1
            data = preprocessed_data;
        end
        data(idx) = data(idx) + y^exp_value;
        transformed_data = data;
    end
   
    
    
end
