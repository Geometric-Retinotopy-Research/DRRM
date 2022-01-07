classdef GenFigureComBoundary < GenFigureBase
	methods
		function [color, result] = setSubRenderColor(obj, static) 
			idvalid = ~isnan(static.feature(:,1)) & ~isnan(static.feature(:,2));
			template_vis = unit_cvt(static.feature,'p2c');
			visxfun =  scatteredInterpolant(static.uv(idvalid,:), template_vis(idvalid,1),'natural'); % rampa
			visyfun =  scatteredInterpolant(static.uv(idvalid,:), template_vis(idvalid,2),'natural');
			eccang = unit_cvt([visxfun(obj.vertex) visyfun(obj.vertex)], 'c2p');
			angchop = eccang(:,2); angchop(eccang(:,1)>8) =0;
			color = prf_value_2_color('lh',angchop);   
			result = obj;
		end
		
		function result = subPlot(obj, static, template_hemi)
			setSubRenderColor(obj, static);
			hold on;
			% extract area boundaries 
			areaids = unique(template_hemi.varea);
			for i=1:length(areaids)
			   id2del = areaids(i)  ~= template_hemi.varea;
			   [fregion,uv,] = gf_remove_mesh_vertices(template_hemi.F, template_hemi.uv, id2del);
			   b = compute_bd(fregion);   
			   plot(uv(b,1),uv(b,2),'-k','linewidth',0.5)   
			end
			axis off;
			result = obj;
		end
	end
end