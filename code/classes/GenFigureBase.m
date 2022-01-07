%% Description  -- GenFigureBase Class
%     base class of Figure handle 
%%
classdef GenFigureBase
% properties: 
%     face [double array]   -- connectivity of mesh
%	  vertex [double array] -- vertex of mesh
%     EValue[double array]  -- extra value of mesh
%	  renderColor[struct]   -- struct of render color include three fields:
%										color -- render color can
%										dircte use
% 
%										value -- value of prf 
%										handle
% 
%										valueType -- type of prf
%										handle
% 
    properties
        face;
        vertex;
        EValue;
        renderColor;
    end
    
    methods
%%
% Description  -- construt function:  
%		function obj = GenFigureBase(param)
%		initial parameters
%
% Parameter(s): 
%     param [struct]  -- including some required parameters to create a figure. 
%
% return: 
%      obj[GenFigureBase]  -- a GenFigureBase after initialize
%%
        function obj = GenFigureBase(param)
            obj.face = param.face;
            obj.vertex = param.vertex;
			if(isfield(param, 'EValue'))
				obj.EValue = param.EValue;
			end
			if(isfield(param, 'renderColor'))
				obj.renderColor = param.renderColor;
			end
        end
	end
	
    methods
%%		
% Description  -- function result = setConfiguration(obj)
%		setup configuration of figure process.
%
% Parameter(s): 
%		obj[GenFigureBase]  --  GenFigureBase object 
%
% return: 
%		result[GenFigureBase]  -- GenFigureBase object
%
% NOTE:
%		need to override this function in sub class if use configuration
%		specificly.
%%
        function result = setConfiguration(obj)
			figure;
			axis off;
            result = obj;
		end
%%		
% Description  -- function result = setRenderColor(obj)
%		setup render color of figure process.
%
% Parameter(s): 
%		obj[GenFigureBase]  --  GenFigureBase object 
%
% return: 
%		result[GenFigureBase]  -- GenFigureBase object
% 
% NOTE:
%		need to override this function in sub class if use render color specificly.
%		
%%
        function result = setRenderColor(obj)
			if (isfield(obj.renderColor, 'maxValue'))
				if (obj.renderColor.maxValue > 0)
					rColor = prf_value_2_color(obj.renderColor.valueType, obj.renderColor.value,obj.renderColor.maxValue);
				else
					rColor = prf_value_2_color(obj.renderColor.valueType, obj.renderColor.value);
				end
			else
				rColor = prf_value_2_color(obj.renderColor.valueType, obj.renderColor.value);
			end
			obj.renderColor.color = rColor;
            result = obj;
		end
%%		
% Description  -- function result = setLightConfiguration(obj)
%		setup light information of figure process.
%
% Parameter(s): 
%		obj[GenFigureBase]  --  GenFigureBase object 
%
% return: 
%		result[GenFigureBase]  -- GenFigureBase object
% 
% NOTE:
%		need to override this function in sub class if use light information specificly.
%		
%%
        function result = setLightConfiguration(obj)
            result = obj;
		end
%%		
% Description  -- result = draw(obj, varargin)
%		draw processing with Configuration and light Configuration.
%
% Parameter(s): 
%		obj[GenFigureBase]  --  GenFigureBase object 
%		varargin -- if there are more than 1 varargin, execute the function
%		in sub classes.
% return: 
%		result[GenFigureBase]  -- GenFigureBase object
% 
% NOTE:
%		need to override this function in sub class if use custom draw approach specificly.
%		
%%
        function result = draw(obj, varargin)
            obj = obj.setRenderColor().setConfiguration().setLightConfiguration();
            if (nargin == 1)
                plot_surf(obj.face,obj.vertex,obj.renderColor.color); 
			end
			axis off;
            result = obj;
		end
%%
% Description  -- function result = save(obj, name)
%		save the current figure to figure path
%
% Parameter(s): 
%		obj[GenFigureBase]  --  GenFigureBase object 
%		name[String] -- file name which you want to save.
% return: 
%		result[GenFigureBase]  -- GenFigureBase object
% 
%%
        function result = save(obj, name)
            saveas(gca,['../Figure/subfigs/', name, '.svg']);
            result = obj;
		end
		function result = setTitle(obj, text)
			title(text);
			result = obj;
		end
    end
end