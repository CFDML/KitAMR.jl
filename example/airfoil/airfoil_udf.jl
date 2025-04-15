function shock_wave_region(midpoint,ds,global_data,level)
    if midpoint[1]>-0.4&&midpoint[1]<1.2&&midpoint[2]<0.4&&midpoint[2]>-0.4&&level<3
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