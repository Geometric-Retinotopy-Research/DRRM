classdef GenFigureRegHCP < GenFigureBase
	methods
		function result = setRenderColor(obj)
			obj = setRenderColor@GenFigureBase(obj);
			if (isfield(obj.renderColor, 'ecc'))
				id =  isnan(obj.renderColor.ecc);
				obj.renderColor.color(id,:) = obj.renderColor.color(id,:)*0+[200 160 140]/255;
			end
			result = obj;
		end
		
		% % extract area boundaries 
		function result = subPlot(obj, template)
			hold on;
			areaids = unique(template.varea);
			for i=1:length(areaids)
			   id2del = areaids(i)  ~= template.varea;
			   [fregion,uv,] = gf_remove_mesh_vertices(template.F, template.uv, id2del);
			   b = compute_bd(fregion);   
			   plot(uv(b,1),uv(b,2),'-w','linewidth',0.5);
			end
			result = obj;
		end
	end
end
