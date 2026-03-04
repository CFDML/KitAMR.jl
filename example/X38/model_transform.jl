using FileIO, MeshIO, GeometryBasics
mesh = load("./example/X38/X38-surface-fine.stl")
points = [x[i] for i in 1:3, x in mesh.position]
mins = minimum(points,dims=2)
maxes = maximum(points,dims=2)
L = maximum(maxes-mins)
α = -π/9# AOA 20 degrees
Rot = [cos(α) -sin(α) 0;sin(α) cos(α) 0;0 0 1.]
normalized_points = Rot*points./L
normalized_position = [Point(x[1],x[2],x[3]) for x in eachcol(normalized_points)]
new_mesh = Mesh(normalized_position,mesh.faces)
save("./example/X38/X38_normalized.stl",new_mesh)
