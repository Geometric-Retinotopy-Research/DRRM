classdef GenFigureRegSyn < GenFigureBase
    properties (Access = private)
        noise = 0;
    end
    methods
        function obj = GenFigureRegSyn(param)
            obj = obj@GenFigureBase(param);
        end
        
        function result = setConfiguration(obj)
			obj = setConfiguration@GenFigureBase(obj);
            xlim([-6,-2]);
            ylim([-2.5 2.5]);
            colormap('hsv');
            result = obj;
        end
%%		
% Description  -- function result = overlayPlot(obj)
%		overlay to plot line on the figure
%
% Parameter(s): 
%		obj[GenFigureBase]  --  GenFigureBase object 
% 
% return: 
%		result[GenFigureBase]  -- GenFigureBase object
%		
%%
		function result = overlayPlot(obj)
			hold on; axis off;
            if (obj.noise == 0)
                overlay_levelset(gca,obj.face,obj.vertex,obj.EValue.Vertex_vis);
            else 
                overlay_levelset(gca,obj.face,obj.vertex, obj.noise);
                obj.noise = 0;
            end
            result = obj;
        end
%%		
% Description  -- function result = subPlot(obj, xy)
%		subplot xy to the figure 
%
% Parameter(s): 
%		obj[GenFigureBase]  --  GenFigureBase object 
%		xy[double array] -- parameters of subplot. 
% 
% return: 
%		result[GenFigureBase]  -- GenFigureBase object
%		
%%       
		function result = subPlot(obj, xy)
            plot(xy(:,1), xy(:,2), 'wo','MarkerSize',6,'Linewidth',2);
            obj.addDecToPlot(xy);
            result = obj; 
        end
%%		
% Description  -- function result = addDecToPlot(obj, xy)
%		add description(text number) on the point
%
% Parameter(s): 
%		obj[GenFigureBase]  --  GenFigureBase object
%		xy[double array] -- parameters of subplot. 
% 
% return: 
%		result[GenFigureBase]  -- GenFigureBase object
%		
%%  
		function result = addDecToPlot(obj, xy)
            text(xy(:,1)+0.05, xy(:,2), (num2str((1:size(xy,1))')),'color','k','FontSize',16);
            result = obj;
        end
        
        function result = addSmallNoisy(obj)
            rng(0)
            Enoise = obj.EValue.Vertex_vis;
            Enoise(:,1:2) = Enoise(:,1:2) + 0.173*rand(size(Enoise(:,1:2)));
            psnr(obj.EValue.Vertex_vis(:,1:2), Enoise(:,1:2))
            obj.noise = Enoise;
            result = obj; 
        end
        
        function result = addStrongNoisy(obj)
            rng(0)
            Enoise = obj.EValue.Vertex_vis;
            Enoise(:,1:2) = Enoise(:,1:2) + 0.5473*rand(size(Enoise(:,1:2)));
            psnr(obj.EValue.Vertex_vis(:,1:2), Enoise(:,1:2))
            obj.noise = Enoise;
            result = obj; 
        end
%%		
% Description  -- function result = addQuiver(obj, roii)
%		add description(text number) on the point
%
% Parameter(s): 
%		obj[GenFigureBase]  --  GenFigureBase object 
% 		roii[int array] -- index of quiver
% 
% return: 
%		result[GenFigureBase]  -- GenFigureBase object
%		
%%  
        function result = addQuiver(obj, roii)
			hold on; axis off;
            quiver(obj.vertex(roii,1), obj.vertex(roii,2), obj.EValue.Vertex_targe(roii,1) - obj.vertex(roii,1), obj.EValue.Vertex_targe(roii,2)-obj.vertex(roii,2),1.0, 'Linewidth',1.2,'color',[0 0 0], 'MaxHeadSize',5);
            result = obj; 
		end
		
		function result = draw(obj)
			obj = draw@GenFigureBase(obj);
			alpha(0.8);
			view(0,-90);
			result = obj;
		end

    end
end