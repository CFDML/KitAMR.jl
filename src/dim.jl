
if DIM === 2 && NDF === 2
    const discrete_maxwell = discrete_maxwell_2DF_2D
    const shakhov_part = shakhov_part_2DF_2D
    const calc_fwb = calc_fwb_2DF_2D
    const calc_w0 = calc_w0_2DF_2D
    const calc_flux_f0 = calc_flux_f0_2DF_2D
    const calc_qf = calc_qf_2DF_2D
end
if DIM === 2
    #=
    2D Z-ordering
    y
    |-----|-----|
    |  3  |  4  |
    |-----|-----|
    |  1  |  2  |
    |-----|-----| x
    =#
    const moment_u = moment_u_2D
    const micro_slope = micro_slope_2D
    const rft = [[2, 4], [1, 4], [3, 2], [1, 3]] # refine face table    
    const rmt = [[-1, -1], [1, -1], [-1, 1], [1, 1]] # refine midpoint table
    const cft = [1, 2, 1, 3]
    const nft = [[1, 3], [2, 4], [1, 2], [3, 4]] # face neighbor inner-sided table, faceid{z-ordering}
end
if NDF === 2
end
