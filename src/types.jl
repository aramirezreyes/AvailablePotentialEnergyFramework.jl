struct FrameWithDetectedCyclones{T,U, V, W <: Union{Nothing,P} where P} 
    count :: T
    labels  :: Array{U}
    centers :: Array{CartesianIndex{V}}
    segmented_frame :: W
end
