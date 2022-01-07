- perform_fast_marching_mesh

  <span style="color: orange">**function [D,S,Q] = perform_fast_marching_mesh(vertex, faces, start_points, options)**</span>

```matlab
% input: vertex, face, start_points, options of mesh processing
% output: D the distance function to the set of starting points.
%   			S is the final state of the points : -1 for dead (ie the distance
%   	    has been computed), 0 for open (ie the distance is only a temporary
%     	  value), 1 for far (ie point not already computed). Distance function
%       	for far points is Inf.
% 			  Q is the index of the closest point. Q is set to 0 for far points.
[D,S,Q] = perform_fast_marching_mesh(Vertex, Face, start);
```

* gf_remove_mesh_vertices

<span style="color: orange">**function [Fout, Vout, new2OldMap] = gf_remove_mesh_vertices(F, V, id2delete)**</span>

````matlab
% input: Face, Vetex and the removal of mesh position
% output: Fout: face after removal operation
% 				Vout: vetex after removal operation
%         new2OldMap: map from old indices to new indices
[Fout,Vout,father] = gf_remove_mesh_vertices(Ffull,Vfull, ind2del); 
````

