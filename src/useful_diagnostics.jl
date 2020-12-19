"""

"""
function average_precipitation_per_pw_bin_dayang(pw,precipitation,max_precip_bin,binspacing)
    precip_bin = Float64[]
    pw_bins = 0:binspacing:max_precip_bin
    n_pw_interval = length(pw_bins)
    for i in 1:n_pw_interval
        pw_i_low = pw_bins[1] + binspacing*(i-1)
        pw_i_up  = pw_bins[1] + binspacing*i
        ind_i    = findall(x -> (x>pw_i_low)&(x<pw_i_up),pw)
        if isempty(ind_i) 
            push!(precip_bin,0.0) 
        else 
            push!(precip_bin,mean(precipitation[ind_i]))
        end
    end
    return precip_bin
end

"""
    average_precipitation_per_pw_bin(pw,precipitation,bins,binspacing)
Computes the mean precipitation as a function of the binned precipitable water following Yang, D., 2018: Boundary Layer Height and Buoyancy Determine the Horizontal Scale of Convective Self-Aggregation. J. Atmos. Sci., 75, 469–478, https://doi.org/10.1175/JAS-D-17-0150.1.
It will be happy if bins is quite large (example from 0 to 300 mm).
"""
function average_precipitation_per_pw_bin(pw,precipitation,pw_bins,binspacing)
        probability_of_being_in_bin = zeros(length(pw_bins))
        average_precipitation_per_bin = zeros(length(pw_bins))
        number_of_columns = reduce(*,size(pw))
        for column in CartesianIndices(pw)
            binindex = ceil(Int,pw[column]/binspacing)
            #@show binindex
            average_precipitation_per_bin[binindex] += precipitation[column] 
            probability_of_being_in_bin[binindex] += 1
        end
        #@show probability_of_being_in_bin
        for (ind,prob) in enumerate(probability_of_being_in_bin)
            prob == 0.0 && continue
            average_precipitation_per_bin[ind] /= prob
            probability_of_being_in_bin[ind] /= (number_of_columns * binspacing)
        end
        return probability_of_being_in_bin,average_precipitation_per_bin
end


"""
    velocity_topolar(u,v,index,center)
Take a velocity vector, an origin of said vector and a center and return the tangential and azimuthal velocities with respect to that center.
"""
function velocity_topolar(u,v,index,center)
    pos = index .- center
    theta1 = atan(v,u)
    theta2 = atan(pos[2],pos[1])
    return -hypot(u,v) * cos(theta1 + theta2), hypot(u,v)*sin(theta1 + theta2)
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
