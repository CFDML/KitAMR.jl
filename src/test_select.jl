using LinearAlgebra
points = [[Float64(i),Float64(j)] for i in 0:4 for j in 0:4]
ref_dir = [0.,1.0]
weights = [dot(x,ref_dir)/norm(x)^2 for x in points]
id = sortperm(weights;rev=true)
sorted = @views points[id]
function make_basis(x)
    return [x[1]^2,x[2]^2,x[1],x[2],1.0]
end
function select_linear_independent_rows(points)
    k = 0
    index = @SVector zeros(Int,5)
    selected = Nothing
    for i in eachindex(points)
        if k == 0
            selected = make_basis(points[i])
            k = 1
            index[k] = i
        else
            temp = hcat(selected,make_basis(points[i]))
            r = rank(temp)
            if r > k
                selected = temp
                k = r
                index[k] = i
            end
        end
        if k == 5
            break
        end
    end
    
    return index,permutedims(selected)
end
select_linear_independent_rows(sorted)