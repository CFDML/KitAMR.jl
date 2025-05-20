function shock_wave_region(midpoint,ds,global_data,level)
    if midpoint[1]>-1.0&&midpoint[1]<3.0&&midpoint[2]<1.0&&midpoint[2]>-1.0&&level<4
        return true
    end
    return false
end
# function shock_wave_region(midpoint,ds,global_data)
#     if midpoint[1]>-1.5&&midpoint[1]<1.0&&midpoint[2]<1.0&&midpoint[2]>-1.0&&√sum(midpoint.^2)>0.5
#         return true
#     end
#     return false
# end