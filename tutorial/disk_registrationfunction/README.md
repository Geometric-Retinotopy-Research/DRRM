* disk_registration

  <span style="color: orange">**function [map,map_mu] = registration(moving, static, landmark, target, maxiter, maxfixstep)**</span>

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

  