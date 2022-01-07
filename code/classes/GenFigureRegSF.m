classdef GenFigureRegSF < GenFigureBase
	methods
		function result = setRenderColor(obj)
			if (~isfield(obj.renderColor, 'color'))
				obj = setRenderColor@GenFigureBase(obj);
				if (isfield(obj.renderColor, 'ecc'))
					id =  isnan(obj.renderColor.ecc);
					obj.renderColor.color(id,:) = obj.renderColor.color(id,:)*0+[200 160 140]/255;
				end
			end
			result = obj;
		end
		function result = setConfigurationExtra(obj)
			camva(8.173459);
			set(gca,'CameraPosition',[822.124821 -1359.197318 -425.643846])
			set(gca,'Position',[0.130000 0.110000 0.775000 0.815000])
			set(gca,'Xlim',[-49.936703 49.933779])
			set(gca,'Ylim',[-131.968077 131.985722])
			set(gca,'Zlim',[-84.466962 84.453975])
			set(gcf,'position',[479.000000 414.000000 560.000000 420.000000])
			hLegend = findobj(gcf, 'Type', 'Legend');
			axis off;
			material([0.5 0.5 0.5]) 
			result = obj;
		end
		function result = setLightConfigurationExtra(obj)
			light('Position',[500 0 0],'Style','local')
			light('Position',[0 500 0],'Style','local')
			light('Position',[0  0 500],'Style','local')
			light('Position',[-500 0 0],'Style','local')
			light('Position',[0 -500 0],'Style','local')
			light('Position',[0  0 -500],'Style','local')
			result = obj;
		end
		function result = subPlot(obj, F, V, colorType, colorValue, maxValue)
			hold on;
			if exist('maxValue', 'var')
				[color, object] = subPlotColor(obj, colorType, colorValue, maxValue);
			else
				[color, object] = subPlotColor(obj, colorType, colorValue, -1);
			end
			plot_surf(F, V, color);
			axis off;
			result = object;
		end
		function [color, result] = subPlotColor(obj, colorType, colorValue, maxValue)
			if (maxValue ~= -1)
				color = prf_value_2_color(colorType, colorValue, maxValue);
			else
				color = prf_value_2_color(colorType, colorValue);
			end
			result = obj;
		end
		function result = subPlotBoundary(obj, tem_hemi)
			hold on;
			areaids = unique(tem_hemi.varea);
			for i=1:length(areaids)
				id2del = areaids(i)  ~= tem_hemi.varea;
				[fregion,uv,] = gf_remove_mesh_vertices(tem_hemi.F,tem_hemi.uv, id2del);
				b = compute_bd(fregion);
				plot(uv(b,1),uv(b,2),'-w','linewidth',0.5)
				result = obj;
			end
			axis off;
		end
	end
end

