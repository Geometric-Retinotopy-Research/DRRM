classdef GenFigure_g < GenFigureBase
   
    methods
        function obj = GenFigure_g(param)
            obj@GenFigureBase(param);
        end
    end
    
    methods

		function result = setRenderColor(obj)
			obj = setRenderColor@GenFigureBase(obj);
            id = find(obj.EValue.Vertex_prf(:,5)<1);
            obj.renderColor.color(id,:)= repmat([189 140 124 ]/255,length(id),1); 
            result = obj;
		end

		function result = setConfiguration(obj)
			obj = setConfiguration@GenFigureBase(obj);
            camva(9.156100);
            set(gca,'CameraPosition',[685.788014 -940.704954 46.745327]);
            set(gca,'Position',[0.130000 0.110000 0.775000 0.815000]);
            set(gca,'Xlim',[-70.716545 3.646090]);
            set(gca,'Ylim',[-107.070717 73.462669]);
            set(gca,'Zlim',[-49.762669 79.694069]);
            set(gcf,'position',[799.000000 483.000000 560.000000 420.000000]);
            result = obj;
        end
        
        function result = setLightConfiguration(obj)
			hold on;
            light('Position',[500 500 0],'Style','local');
            light('Position',[-500 500 0],'Style','local');
            light('Position',[0 500 0],'Style','local');
            light('Position',[0 -500 0],'Style','local');
			alpha(0.8)
            result = obj;
        end
        function result = draw(obj)
            obj = draw@GenFigureBase(obj,[]);
            plot_surf(obj.face,obj.vertex,'FaceVertexCData',obj.renderColor.color, 'Edgecolor','none','FaceLighting','gouraud',...
                    'AmbientStrength',0.5);
			axis off;
            result = obj;
        end
    end
end

