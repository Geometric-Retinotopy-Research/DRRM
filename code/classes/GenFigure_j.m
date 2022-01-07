classdef GenFigure_j< GenFigureBase
    properties
        father;
    end
    
    methods
        function obj = GenFigure_j(param)
            obj = obj@GenFigureBase(param);
            obj.father = param.father;
        end
    end
    
    methods
        function result = setRenderColor(obj)
			obj = setRenderColor@GenFigureBase(obj);
            id = find(obj.EValue.Vertex_prf(:,5)<1);
            obj.renderColor.color(id,:)= repmat([189 140 124 ]/255,length(id),1); 
            result = obj;
        end
        function result = draw(obj)
            obj = draw@GenFigureBase(obj,[]);
            plot_surf(obj.face,obj.vertex,'FaceVertexCData',obj.renderColor.color(obj.father,:));
            result = obj;
        end
        
        function result = sub_plot(obj,areaids)
			hold on;
			for i=1:length(areaids)
				id2del = areaids(i) ~= obj.EValue.Vertex_atlashcp(obj.father);
				[fregion,uvb,father_r] = gf_remove_mesh_vertices(obj.face,obj.vertex, id2del);   
				b = compute_bd(fregion);   
				plot(uvb(b,1),uvb(b,2),'-k','linewidth',1)
				mu = compute_bc(fregion, uvb, obj.EValue.Vertex_prf(obj.father(father_r),1:2));
				if i==1  
					face = fregion(abs(mu)>1,:);
				else
					face = fregion(abs(mu)<1,:);
				end
				plot_surf(face,uvb,'FaceColor','k','EdgeColor','k');
			end
			axis off;
            result = obj;
        end
        
    end
end

