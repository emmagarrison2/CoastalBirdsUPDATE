#### Color Palettes Using Colorspace Package ####

library(colorspace)

# so much useful information about this package and picking colors here: https://colorspace.r-forge.r-project.org/index.html

# general idea: there are 3 main types of color palettes for data visualization

# 1) Qualitative: Designed for coding categorical information, i.e., where no particular ordering of categories is available and every color should receive the same perceptual weight.
# 2) Sequential: Designed for coding ordered/numeric information, i.e., going from high to low (or vice versa).
# 3) Diverging: Designed for coding ordered/numeric information around a central neutral value, i.e., where colors diverge from neutral to two extremes.

# colorspace has a number of pre-built palettes of each of these types

# see the entire list of palettes:
hcl_palettes()

##### qualitative palettes #####
# plot all the qualitative palettes:
hcl_palettes("qualitative", plot=T)

##### sequential palettes #####
# for sequential palettes, there are single hue and multi hue

# plot all the single-hu sequential palettes:
hcl_palettes("sequential (single-hue)", plot = TRUE)
# by default, these plot 5 shades, increase the number using n=
hcl_palettes("sequential (single-hue)", n=7, plot = TRUE) # this example shows 7 shades per palette

# plot all the multi-hu sequential palettes:
hcl_palettes("sequential (multi-hue)", plot = TRUE)
# by default, these plot 5 shades, increase the number using n=
hcl_palettes("sequential (multi-hue)", n=7, plot = TRUE) # this example shows 7 shades per palette

##### diverging palettes #####
# plot all the diverging palettes:
hcl_palettes("diverging", plot = TRUE)

# show them with 9 color each:
hcl_palettes("diverging", n=9, plot = TRUE)
# NOTE: diverging palettes typically have an ODD number of colors


### How to Use the Palettes ####

# Lets say you want to use the sequential multi-hue palette BluYl

# you can plot it on its own with the number of colors you require like this:
specplot(sequential_hcl(6, "BluYl"), plot=T) # this shows it with 6 colors. Change the number to plot more or less
swatchplot()

# depending on the type of palette you want to look at in more detail, the function you should use here changes
# sequential_hcl() for sequential palettes
# diverging_hcl() for diverging palettes
# qualitative_hcl() for the qualitative palettes

# to pull the hex codes for the colors, use this:
sequential_hcl(6, "BluYl")
# you can use these hex codes in ggplot2

#### demo plots ####
# there are some fun ways to quickly demo the palette to see how it might look
# these examples show "typical" types of plots where each kind of palette might be used

# for sequential
demoplot(sequential_hcl(50, "BluYl"), type = "heatmap")
demoplot(sequential_hcl(50, "BluYl"), type = "perspective")
demoplot(sequential_hcl( 4, "BluYl"), type = "spine")

# for qualitative
demoplot(qualitative_hcl(4, "Set 2"), type = "pie")
demoplot(qualitative_hcl(4, "Set 2"),    type = "scatter")
demoplot(qualitative_hcl(4, "Set 2"),   type = "lines")

# for diverging
demoplot(diverging_hcl(50, "Purple-Green", power = 2.5), type = "map")
demoplot(diverging_hcl( 11, "Purple-Green"), type = "mosaic")
demoplot(diverging_hcl( 7, "Purple-Green"), type = "bar")


# there are also more sophisticated ways to use the palettes in ggplot that are shown with examples here:
https://colorspace.r-forge.r-project.org/articles/ggplot2_color_scales.html

  