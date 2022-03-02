# # Tomography

#-

# Tomography is the process of reconstructing a density distribution from given
# integrals over sections of the distribution. In our example, we will
# work with tomography on black and white images.
# Suppose $x$ be the vector of $n$ pixel densities, with $x_j$
# denoting how white pixel $j$ is.
# Let $y$ be the vector of $m$ line integrals over the image, with $y_i$
# denoting the integral for line $i$.
# We can define a matrix $A$ to describe the geometry of the lines. Entry
# $A_{ij}$ describes how much of pixel $j$ is intersected by line $i$.
# Assuming our measurements of the line integrals are perfect, we have the relationship that
#
# $$
#   y = Ax
# $$
#
# However, anytime we have measurements, there are usually small errors that occur.
# Therefore it makes sense to try to minimize
#
# $$
#  \|y - Ax\|_2^2.
# $$
#
# This is simply an unconstrained least squares problem; something we can readily solve!

using Convex, ECOS, DelimitedFiles, SparseArrays
aux(str) = joinpath(@__DIR__, "aux_files", str) # path to auxiliary files
line_mat_x = readdlm(aux("tux_sparse_x.txt"))
summary(line_mat_x)

line_mat_y = readdlm(aux("tux_sparse_y.txt"))
summary(line_mat_y)

line_mat_val = readdlm(aux("tux_sparse_val.txt"))
summary(line_mat_val)

line_vals = readdlm(aux("tux_sparse_lines.txt"))
summary(line_vals)

# Form the sparse matrix from the data
# Image is 50 x 50
img_size = 50

# The number of pixels in the image
num_pixels = img_size * img_size

line_mat = spzeros(length(line_vals), num_pixels)

num_vals = length(line_mat_val)

for i in 1:num_vals
    x = Int(line_mat_x[i])
    y = Int(line_mat_y[i])
    line_mat[x+1, y+1] = line_mat_val[i]
end

pixel_colors = Variable(num_pixels)
## line_mat * pixel_colors should be close to the line_integral_values
## to reflect that, we minimize a norm
objective = sumsquares(line_mat * pixel_colors - line_vals)
problem = minimize(objective)
solve!(problem, ECOS.Optimizer; silent_solver = true)

rows = zeros(img_size * img_size)
cols = zeros(img_size * img_size)
for i in 1:img_size
    for j in 1:img_size
        rows[(i-1)*img_size+j] = i
        cols[(i-1)*img_size+j] = img_size + 1 - j
    end
end

# Plot the image using the pixel values obtained!

using Plots
image = reshape(evaluate(pixel_colors), img_size, img_size)
heatmap(
    image,
    yflip = true,
    aspect_ratio = 1,
    colorbar = nothing,
    color = :grays,
)
