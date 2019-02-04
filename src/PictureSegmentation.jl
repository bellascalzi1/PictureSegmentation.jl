module PictureSegmentation

using Images

abstract type SegmentationAlgorithm end
struct SLIC_superpixels <: SegmentationAlgorithm end

include("new_slic.jl")

export
	# main functions
	segment_image,
	SLIC_superpixels
end # module
