"""
    average_precipiation_per_pw_bin(pw,precipitation,bins,binspacing)
Computes the mean precipitation as a function of the binned precipitable water following Yang, D., 2018: Boundary Layer Height and Buoyancy Determine the Horizontal Scale of Convective Self-Aggregation. J. Atmos. Sci., 75, 469–478, https://doi.org/10.1175/JAS-D-17-0150.1.
It will be happy if bins is quite large (example from 0 to 300 mm).
"""
function average_precipiation_per_pw_bin(pw,precipitation,bins,binspacing)
        probability_of_being_in_bin = zeros(BigInt,length(bins))
        average_precipitation_per_bin = zeros(BigFloat,length(bins))
        number_of_columns = reduce(*,size(pw))
        for column in CartesianIndices(pw)
            binindex = ceil(Int,pw[column]/binspacing)
            average_precipitation_per_bin[binindex] += precipitation[column] 
            probability_of_being_in_bin[binindex] += 1
        end
        average_precipitation_per_bin = average_precipitation_per_bin./probability_of_being_in_bin
        probability_of_being_in_bin = probability_of_being_in_bin./number_of_columns./ (binspacing)
        return probability_of_being_in_bin,average_precipitation_per_bin
end



"""
    velocity_cartesian_to_polar(u,v,index_of_point = (0,0),index_of_center = (0,0))
    Take matrix of U and V components of velocities as well as a matrix of locations and the location of a center.
    Compute the tangential and radial components of velocity around this center.
"""

function velocity_cartesian_to_polar(u,v,index_of_point = (0,0),index_of_center = (0,0))
    tangential = similar(u)
    radial = similar(u)
    for index in CartesianIndices(u)
        t,r = velocity_cartesian_to_polar(u[index],v[index],Tuple(index)[1:2],index_of_center[1:2])
        radial[index] = r
        tangential[index] = t
    end
    return tangential,radial
end


function velocity_cartesian_to_polar(u :: Number,v :: Number,index_of_point = (0,0),index_of_center = (0,0))
    speed = hypot(u,v)
    x,y = index_of_point .- index_of_center
    angle_of_position_vector = atan(y,x)
    angle_of_velocity_vector = atan(v,u)
    return velocity_polar = speed.*sincos( angle_of_velocity_vector - angle_of_position_vector )
end
