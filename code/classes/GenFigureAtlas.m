classdef GenFigureAtlas < GenFigureBase
    methods
		function obj = GenFigureAtlas(param)
            obj = obj@GenFigureBase(param);
		end
		
		function result = setRenderColor(obj)
			if (~isfield(obj.renderColor, 'color'))
				obj = setRenderColor@GenFigureBase(obj);
				if (isfield(obj.renderColor, 'ecc'))
					obj.renderColor.color(obj.EValue.ecc==0,:) = obj.renderColor.ecc;
				end
			end
			result = obj;
		end
		
		function result = setConfigurationExtra(obj)
			camva(8.590709);
			set(gca,'CameraPosition',[965.264891 -1364.624697 410.385565]);
			material dull;
			result = obj;
		end
		function result = setLightConfigurationExtra(obj)
			light('Position',[965.264891 -1364.624697 410.385565],'Style','local');
			result = obj;
		end
		function result = subPlot(obj, colorType, colorValue)
			[color, obj] = subPlotColor(obj, colorType, colorValue);
			plot_surf(obj.face,obj.vertex,color); 
			axis off;
			result = obj;
		end
		function [color, result] = subPlotColor(obj, colorType, colorValue)
			color = prf_value_2_color(colorType, colorValue);
			color(obj.EValue.ecc==0,:) = 0.8; 
			result = obj;
		end
	end
end


