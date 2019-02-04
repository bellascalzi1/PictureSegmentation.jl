"""
Documentation goes here
"""
function segment_image(algorithm::SLIC_superpixels,imgIn, k, maxiter)
    img = Lab{Float64}.(imgIn)
    nrow, ncol = size(img)
    N = length(img)
    s = floor(Int, sqrt(N / k))
    error = 0.0

    # Initalise arrays for the first time.
    counts = zeros(k)
    totals = zeros(k, 5)
    centers = Array{CartesianIndex}(undef, k)
    center_colors = Array{Lab{Float64}}(undef, k)
    pixel_associactions = zeros(Int, nrow, ncol)
    distances = fill(Inf, size(img))

    # Initalise the cluster centers.
    initialise_centers!(centers, s, k, nrow, ncol)
    seed_locations!(centers, img)

    # Run the segmentation algorithm to generate the cluster centers.
    for i = 1:maxiter - 1
        assign_labels!(centers,center_colors, img, pixel_associactions, counts, distances, totals, s)
        centers,center_colors,img = calc_new_centers(centers,center_colors,counts,totals,img)
        fill!(counts,0)
        fill!(totals,0.0)
        fill!(pixel_associactions,0)
        fill!(distances,Inf)
    end

    # Run the segmentation algorithm one more time to retain data.
    assign_labels!(centers, center_colors, img, pixel_associactions, counts, distances, totals,s)
    centers,center_colors,img=calc_new_centers(centers,center_colors,counts,totals,img)
    # Enforce Connectivity
    pixel_associactions = enforce_connectivity(img,s,pixel_associactions, centers, center_colors)

    # Create the superpixels.
    img2 = similar(img)
    for i in CartesianIndices(img2)
        if pixel_associactions[i] !=0
            img2[i] = center_colors[pixel_associactions[i]]
        end
    end

    return img2

end

function calculate_gradient_magnitude(img::AbstractArray)
    img = RGB{Float64}.(img)
    magnitudes = zeros(Float64, size(img))
    for i in CartesianIndices(img)
        RGBr = img[min(i[1]+1,1),i[2]]-img[max(i[1]-1,first(size(img))),i[2]]
        RGBc = img[i[1],min(i[2]+1,1)]-img[i[1],max(i[2]-1,last(size(img)))]
        magnitudes[i] = (RGBr.r^2+RGBr.g^2+RGBr.b^2) + (RGBc.r^2+RGBc.g^2+RGBc.b^2)
    end
    return magnitudes
end

# Calclates colour distance between 2 pixels.
function d_lab(𝑙ₖ::Float64, 𝑙ᵢ::Float64, 𝑎ₖ::Float64, 𝑎ᵢ::Float64, 𝑏ₖ::Float64, 𝑏ᵢ::Float64)
    return sqrt( (𝑙ₖ - 𝑙ᵢ)^2 + (𝑎ₖ - 𝑎ᵢ)^2 + (𝑏ₖ - 𝑏ᵢ)^2 )
end
# Calculates xy distance between 2 pixels.
function d_xy(𝑥ₖ::Real, 𝑥ᵢ::Real, 𝑦ₖ::Real, 𝑦ᵢ::Real)
    return sqrt( (𝑥ₖ - 𝑥ᵢ)^2 + (𝑦ₖ - 𝑦ᵢ)^2 )
end
#=Calculates the top-left and bottom-right bounds of a neighbourhood
of width w x w around the pixel passed in.
=#
function get_neighbourhood(center::CartesianIndex, w::Integer, img::AbstractArray)
    nrow, ncol = size(img)
    r, c  = center.I
    r₀ = max(r - w, 1)
    r₁ = min(r + w, nrow)
    c₀ = max(c - w, 1)
    c₁ = min(c + w, ncol)
    CartesianIndices((r₀:r₁, c₀:c₁))
end

#= Assigns each pixel with its nearest cluster center, and add the pixel's
labxy value to a running total.
=#
function assign_labels!(centers, center_colors, img, pixel_associactions, counts, distances, totals, s)
    for i = 1:length(centers)
        if centers[i] != CartesianIndex(0,0)
            a = center_colors[i]
            c₁, c₂ = centers[i].I
            search_area = get_neighbourhood(centers[i], s, img)
            for pixel in search_area
                b = img[pixel]
                dlab = sqrt((a.l - b.l)^2 + (a.a - b.a)^2 + (a.b - b.b)^2)
                dxy = sqrt((c₂ - pixel[2])^2 + (c₁ - pixel[1])^2)
                temp_dist =  dlab + (10/s) * dxy
                if temp_dist < distances[pixel]
                    r, c = pixel.I
                    # If the pixel was assigned to another center, remove it.
                    if pixel_associactions[pixel] != 0
                        remove_from_centroid!(img[pixel], r, c, pixel_associactions[pixel],  counts, totals)
                    end
                    # Associate pixel with new center.
                    pixel_associactions[pixel] = i
                    # Update shortest distance.
                    distances[pixel] = temp_dist
                    add_to_centroid!(img[pixel], r, c, i, counts, totals)
                end
            end
        end
    end
    # Remove centers that have overlapped each-other.
    for i in eachindex(counts)
        if counts[i]<=0
            centers[i]=CartesianIndex(0,0)
        end
    end
    return pixel_associactions,counts,distances,totals
end

# Remove a pixel from a centroid's running total.
function remove_from_centroid!(pixel, r, c, index,  counts, totals)
    counts[index] -= 1
    totals[index,1] -= pixel.l
    totals[index,2] -= pixel.a
    totals[index,3] -= pixel.b
    totals[index,4] -= r
    totals[index,5] -= c
end

function add_to_centroid!(pixel, r, c, index,  counts, totals)
    counts[index] += 1
    totals[index,1] += pixel.l
    totals[index,2] += pixel.a
    totals[index,3] += pixel.b
    totals[index,4] += r
    totals[index,5] += c
end

function calc_new_centers(centers,center_colors,counts,totals,img)
    # Compute new centers
    for i in eachindex(centers)
        if centers[i] != CartesianIndex(0,0)
            # Assign new centers to their respective arrays.
            center_colors[i] = Lab{Float64}(totals[i,1]/counts[i],totals[i,2]/counts[i],totals[i,3]/counts[i])
            centers[i] = CartesianIndex(round(Int,totals[i,4]/counts[i]),round(Int,totals[i,5]/counts[i]))
        end
    end
    return centers, center_colors, img
end

#= Generate starting locations for centers, at regularly
spaced intervals.
=#
function initialise_centers!(centers, s, K, nrow, ncol)
    k = 0
    for c =2:s:(ncol-1)
        for r = 2:s:(nrow-1)
            #@show k
            if k != K
                k += 1
                centers[k] = CartesianIndex(r,c)
            end
        end
    end
end

function seed_locations!(centers::Array{CartesianIndex}, img::AbstractArray)
    gradient_magnitude = calculate_gradient_magnitude(img)
    nrow, ncol = size(img)
    for i in eachindex(centers)
        r, c  = centers[i].I
        r₀ = max(r - 1, 1)
        r₁ = min(r + 1, nrow)
        c₀ = max(c - 1, 1)
        c₁ = min(c + 1, ncol)
        minval, pos = findmin(gradient_magnitude[r₀:r₁, c₀:c₁])
        centers[i] += (CartesianIndex(pos) - CartesianIndex(r₁ - r₀, c₁ - c₀))
    end
end

function flood_connections(pixel_associactions::Array{Int},pixels::Array{Int})
    changed = true
    while changed == true
        changed = false
        for i in CartesianIndices(pixels)
            if pixels[i]!=0
                if i[1]+1 <= size(pixel_associactions)[1] && pixel_associactions[i[1]+1,i[2]]==pixels[i] && pixels[i[1]+1,i[2]]==0
                    pixels[i[1]+1,i[2]]=pixels[i]
                    changed = true
                end
                if i[1]-1 >= 1 && pixel_associactions[i[1]-1,i[2]]==pixels[i] && pixels[i[1]-1,i[2]]==0
                    pixels[i[1]-1,i[2]]=pixels[i]
                    changed = true
                end
                if i[2]+1 <= size(pixel_associactions)[2] && pixel_associactions[i[1],i[2]+1]==pixels[i] && pixels[i[1],i[2]+1]==0
                    pixels[i[1],i[2]+1]=pixels[i]
                    changed = true
                end
                if i[2]-1 >= 1 && pixel_associactions[i[1],i[2]-1]==pixels[i] && pixels[i[1],i[2]-1]==0
                    pixels[i[1],i[2]-1]=pixels[i]
                    changed = true
                end
            end
        end
    end
    return pixels
end
function assign_unconnected(img::AbstractArray, s::Int64, pixel_associactions::Array{Int}, pixels::Array{Int}, centers::Array{CartesianIndex}, center_colors::Array{Lab{Float64}})
    changed = true
    while changed == true
        changed = false
        for i in CartesianIndices(pixels)
            if pixels[i] == 0
                dxy = 0
                Dmin=10000.0
                minPos=i
                d=100000.0
                unassigned_pixel_color = img[i]
                label = 0
                if i[1]+1 <= size(pixel_associactions)[1] && pixels[i[1]+1,i[2]]!=0
                    label = pixels[i[1]+1,i[2]]
                    b, c = centers[label].I
                    center_color = center_colors[label]
                    dlab = sqrt((unassigned_pixel_color.l - center_color.l)^2 + (unassigned_pixel_color.a - center_color.a)^2 + (unassigned_pixel_color.b - center_color.b)^2)
                    dxy = sqrt((i[2] - c)^2 + (i[1] - c)^2)
                    d =  dlab + (10/s) * dxy
                    if d<Dmin
                        Dmin = d
                        minPos = CartesianIndex(i[1]+1,i[2])
                        changed = true
                    end
                end
                if i[1]-1 >= 1 && pixels[i[1]-1,i[2]]!=0
                    label = pixels[i[1]-1,i[2]]
                    center_color = center_colors[label]
                    dlab = sqrt((unassigned_pixel_color.l - center_color.l)^2 + (unassigned_pixel_color.a - center_color.a)^2 + (unassigned_pixel_color.b - center_color.b)^2)
                    dxy = sqrt((i[2] - centers[label][2])^2 + (i[1] - centers[label][2])^2)
                    d =  dlab + (10/s) * dxy
                    if d<Dmin
                        Dmin = d
                        minPos = CartesianIndex(i[1]-1,i[2])
                        changed = true
                    end
                end
                if i[2]+1 <= size(pixel_associactions)[2] && pixels[i[1],i[2]+1]!=0
                    label = pixels[i[1],i[2]+1]
                    center_color = center_colors[label]
                    dlab = sqrt((unassigned_pixel_color.l - center_color.l)^2 + (unassigned_pixel_color.a - center_color.a)^2 + (unassigned_pixel_color.b - center_color.b)^2)
                    dxy = sqrt((i[2] - centers[label][2])^2 + (i[1] - centers[label][2])^2)
                    d =  dlab + (10/s) * dxy
                    if d<Dmin
                        Dmin = d
                        minPos = CartesianIndex(i[1],i[2]+1)
                        changed = true
                    end
                end
                if i[2]-1 >= 1 && pixels[i[1],i[2]-1]!=0
                    label = pixels[i[1],i[2]-1]
                    center_color = center_colors[label]
                    dlab = sqrt((unassigned_pixel_color.l - center_color.l)^2 + (unassigned_pixel_color.a - center_color.a)^2 + (unassigned_pixel_color.b - center_color.b)^2)
                    dxy = sqrt((i[2] - centers[label][2])^2 + (i[1] - centers[label][2])^2)
                    d =  dlab + (10/s) * dxy
                    if d<Dmin
                        Dmin = d
                        minPos = CartesianIndex(i[1],i[2]-1)
                        changed = true
                    end
                end
                pixels[i]=pixels[minPos]
            end
        end
    end
    return pixels
end
function enforce_connectivity(img::AbstractArray, s::Int64, pixel_associactions::Array{Int}, centers::Array{CartesianIndex}, center_colors::Array{Lab{Float64}})
    ax = size(pixel_associactions)
    pixels = zeros(Int, ax)
    for i in eachindex(centers)
        current_ctr_row,current_ctr_col=centers[i].I
        if centers[i] != CartesianIndex(0,0)
             pixels[current_ctr_row,current_ctr_col] = i
        end
    end
    pixels = flood_connections(pixel_associactions,pixels)
    pixels = assign_unconnected(img, s, pixel_associactions, pixels, centers, center_colors)
    return pixels
end
