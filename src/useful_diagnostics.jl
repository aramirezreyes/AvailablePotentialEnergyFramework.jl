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
