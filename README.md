# DRRM



This project is the core program of Diffeomorphic Registration for Retinotopic Maps of Multiple Visual Regions. 

 * [**Data**](https://osf.io/s25pe/) is availale here. 

---

## Introduction:

### The construction of project:

```c
|- code/       // including code and anything files support the project could run
|- data/       // some data which project need to plot figures
|- Figure/     // figures will save in this floder after plot results
|- README.md   // the document introduce the whole project (current file)
|- LICENSE.txt // LICENSE file
``` 

***details of each folders:***

* **code/:**

  Sub-folders:

  <span style="color: orange">**- classes/:** </span> including classes use to plot figures.

  <span style="color: orange">**- libiaries/:** </span> third-part libiaries provide some extra function efficiently. Such as "BensonRegister", "pRF-decoder", "toolbox_fast_marching" and so on.

  <span style="color: orange">**- resources/:** </span> resources support project, including some sub folders: text files, images generating during runing processing, templates will use during ploting processing, and other resources.

  <span style="color: orange">**- utilities/:** </span>some small utilities provide specific functions to handle problems more effcient and quickly.

  <u>**Other .m file could run directly and show figures**</u>

  

* **data/:** **<u> (required downloading data and put in these folders）</u>**

  Sub-folders:

  * <span style="color: orange">**- HCP/:** </span> The Human Connectome Project (*HCP*).

    * <span style="color: yellow">**- mesh/:** </span> mesh data.

  * <span style="color: orange">**- HCP32Kmat/:** </span> HCP32K data.

  * <span style="color: orange">**- Studyforrest/:** </span> study forrest data.
 
  
  ---  


# How to use our code for your own work?
The retinotopic code consist of geometry processing, fMRI fitting, mathematics concepts. Therefore, it is relatively complicate to adopt without an explaination.  Here we provide some information for starters to begin. 

## Step 1: Get familar with the geometry processing package. 
You can read the tutorials in the geometry processing package. The fundamental data structure is the mesh, which is consisted of a list a vertices and a list of faces. 
For instance,  
`[F, V, E]=read_mfile(mfilename)` 
This read a mesh with face list to F, vertices list to v, and maybe(very ofen) some extra features associated to vertices or faces(typically on the vertex). 

** Whenever we process on the mesh, the input is F,V,E (may stored in a struct) and the output should be Fn,Vn,En. If we change F, we change the connectiivy of the mesh. In our paper, only cutting the cortical surface need change the F. It is very ofen not change the V, but record a new veriable, say uv to do some mapping. 
e.g. `uv = disk_conformal_map(F,V)`
This does not change F, and V. It generate a new vertice list. Then we can use (F,uv) to form a new mesh. In this case, it is the conformal map of input mesh (F,V)
** 

## Step 2: Use our core funtionality: disk_registration
After preparing two meshes on the planar disk. 
Namely, static = (F1,u1,E1) and moving = (F2,u2,E2). 
The disk_registration try to align them diffeomorphiclly. The diffeomorphic condition is to say, no flip of triangle orientation during the alginment. 
We assume the feature is same as the number of vertices. Nxf: N is the number of points, f is the features. Here, we use the eccentricity and polar angle. So f=2. 

  ```matlab
  %input parameters of: static mesh structure, moving mesh structure, specify the landmark including boundary index, target of landmark.
  
  % generate a static struture
  static.face = face;
  static.vertex = vertex;
  static.feature = feature;
  ... ...
  
  % generate a moving structure
  moving.face = face;
  moving.vertex = vertex;
  moving.feature = feature;
  ... ...
  
  % make registeraion for below, the rest of parameters set [] if you don't have these value.
  [result,bc]  = disk_registration(static, moving, [], []);
  
  % results: registration result
  % bc: corresponding Beltrami coefficients
  
  ```
After the execuation of disk_registration. The output is a map of moved points. bc is the associated with the map.


* [**function [map,map_mu] = registration(moving, static, landmark, target, maxiter, maxfixstep)**](https://github.com/colinxiny/BrainStructuralAndFunctional/tree/feature/code_reconstrution/tutorial/disk_registrationfunction)

## Step 3: Prepration & Evaulation performance, 

Most of time,  you don't have mesh on the disk to registere. On the contract, you have a 3D mesh of genus zero. In this case, you do the following:
** Cut the ROI : We provide  `cut_on_sphere_32k`
** Map to disk: We provide  `disk_conformal_map`
Then you the mesh on the disk. 

Evaulation the regsitration results, can be tricky. This requires more knowledge of pRF decoding. 
We recommand the readers to read `compute_fMRI_metrics.m`.

## What to do if I still cannot process my data?
Maybe you can submit a issue in github or contact us. 


---

  ## Requirements:

  * Matlab: 2020b and later version is recommanded. 


---

  ## Architecture：

  * Classes:
  
    * `GenConsts`:  definition of some const in project, like file path.
    
    * `GenDataManager`: provide data from outside resouces to runing processing which need to use, like funtions: **getTemplateTxt()**
    
    * `GenFigureBase`: Base class of figure handle. Provide some basic ability to plot a figure: 
    
      Properties:
    
      * **face** 
      * **vertex**
      * **EValue**
      * **renderColor**
    
      Functions:
    
      - **setConfiguration():** configurations of figures, like camera position, Xlim, and so on.
      - **setLightConfiguration():** set light position of figures.
      - **setRenderColor():** set render color with value and type, then return a pRF color.
      - **draw():** finish plot surf handle with face, vertex and color, before plot, use above configurations automatically.
    
    * These are sub-class from GenFigureBase:
    
      `GenFigure_g` : load hemisphere and find the retina map on the hemisphere.
    
      `GenFigure_j`: conformal map of hemisphere.
    
      `GenFigureRegSyn`: plot register data with synthetic way.
    
      `GenFigureAtlas`: generating template with benson's template.
    
      `GenFigureRegHCP`: plot with HCP template saved before.
    
      `GenFigureRegSF` : plot with study forrest subject.
    
      These classes inherit base and override some function in base classes in order to implement different function seperately.
    
      

---

  ## Usage of GenFigure：

  * create a class inherit the GenFigureBase :

  ```matlab
  classdef GenFigureCustom < GenFigureBase
  ...
  ...
  ...
  end
  ```

  * Besides fixed properties in the base class, claiming some properties which you need:

  ```matlab
      properties
          propertryNameA;
          propertryNameB;
          propertryNameC;
      end
  ```

  * Then, override construction if you want to set value or make some initial handle in before starting use this class. **(Tip: use figure base class constrution funtion with the first line)**

  ```matlab
  methods
      function obj = GenFigureA(param)
      % use base class constrution function before custom opaeration.
      obj = obj@GenFigureBase(param);  
      obj.propertryNameA = param.propertryNameA;
      ...
      ...
      % initial handle in there
      end
  end
  ```

  

  * Implement functions if you need using custom handle, like:
  
    > **setRenderColor()**

  ```matlab
  function result = setRenderColor(obj)
      % use base class setRenderColor(obj) function before custom. Through the function of
      % set renderColor in the base class. We get the color and store in the struct 
      % 'renderColor' with the varible color.
      obj = setRenderColor@GenFigureBase(obj);
      % modifing the color when get the color from base class.
      id = find(obj.EValue.Vertex_prf(:,5)<1);
      obj.renderColor.color(id,:)= repmat([189 140 124 ]/255,length(id),1); 
      result = obj;
  end
  
  %-----------------------------setRenderColor@GenFigureBase(obj)---------------------------
  
  function result = setRenderColor(obj)
      % max value is supported in the base class in setRenderColo().
      if (isfield(obj.renderColor, 'maxValue'))
          if (obj.renderColor.maxValue > 0)
          % using type and value to calculate the prf color with the function
          % 'prf_value_2_color'.
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
  ```

  

  > **setConfiguration()**

  ```matlab
  function result = setConfiguration(obj)
      % optional
      obj = setConfiguration@GenFigureBase(obj);
      % figure;
      camva(9.156100);
      set(gca,'CameraPosition',[685.788014 -940.704954 46.745327]);
      set(gca,'Position',[0.130000 0.110000 0.775000 0.815000]);
      set(gca,'Xlim',[-70.716545 3.646090]);
      set(gca,'Ylim',[-107.070717 73.462669]);
      set(gca,'Zlim',[-49.762669 79.694069]);
      set(gcf,'position',[799.000000 483.000000 560.000000 420.000000]);
      result = obj;
  end
  ```

  > **setLightConfiguration()**

  ```matlab
  % this function don't need use base funtion setLightConfiguration(), But please keep the 
  % name of function correctly "setLightConfiguration"
  function result = setLightConfiguration(obj)
      hold on;
      light('Position',[500 500 0],'Style','local');
      light('Position',[-500 500 0],'Style','local');
      light('Position',[0 500 0],'Style','local');
      light('Position',[0 -500 0],'Style','local');
      alpha(0.8)
      result = obj;
  end
  ```

  > **draw()**

  ```matlab
  function result = draw(obj)
      % the second parameter pass to a non-value parameter('[]') to skip the draw function in 
      % the base class and force use the custom function.
      obj = draw@GenFigureBase(obj,[]);
      plot_surf(obj.face,obj.vertex,'FaceVertexCData',obj.renderColor.color, 		'Edgecolor','none','FaceLighting','gouraud',...
    'AmbientStrength',0.5);
      axis off;
      result = obj;
  end
  
  %-----------------------------setRenderColor@GenFigureBase(obj)---------------------------
  function result = draw(obj, varargin)
      % set the cofiguration with chain call. (Automatically find the sub-class and the base
      % class)
      obj = obj.setRenderColor().setConfiguration().setLightConfiguration();
      if (nargin == 1)
          plot_surf(obj.face,obj.vertex,obj.renderColor.color); 
      end
      axis off;
      result = obj;
  end
  ```

  *In addtion, you can set any functions if you need in the sub-class.*

  

  ### How to use your own figure class?

  * basic usage: 

  ```matlab
  % basic parameters
  param.face = valueFace;
  param.vertex = valueVertex;
  param.EValue = valueEValue;
  
  % render color settingz
  param.renderColor.valueType = colorType;
  param.renderColor.value = colorValue;
  
  % generate a Figure class
  figure = GenFigureCustom(param);
  figure.draw();
  ```

  * with maxValue in render color:

  ```matlab
  param.renderColor.valueType = colorType;
  param.renderColor.value = colorValue;
  param.renderColor.maxValue = colorMaxValue;
  ```

  * with extra configuration function or sub plot in the current figure:

  ```matlab
  % cutom cofiguration
  figure.draw().extraConfigurationA(). ... ... .extraConfigurationZ();
  
  % sub plot
  figure.draw().subPlotA(). ... ... .subPlotZ();
  
  ```

  * base figure class support some basic ability like:

  ```matlab
  % save figure to the specific file path (the path defined in the GenConsts.m)
  % the type of file maybe different from different platform(Mac OS, Windows)
  % you can change the suffix of the file name like 'figure.csv'
  figure.draw().save('figure.png');
  
  % set title
  figure.draw().setTitle('figureName');
  ```

---

  ## License

  Released under the [ASU GitHub Project License](https://github.com/Retinotopy-mapping-Research/DRRM/blob/master/LICENSE.txt).
