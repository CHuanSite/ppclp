# ppclp

A R package to implement the probabilistic principal curve with length penalty and fixed starting and ending point. This package includes both the 2D and 3D version. 

## Shiny
A shiny app has been built for this package, please refer to [ppclpShiny](https://ppclp.shinyapps.io/ppclpShiny/)

## Functions
```ppclp2D```:This function implements the probabilistic principal curve with length penalty(ppclp) algorithm in 2D image

```ppclp3D```:This function implements the probabilistic principal curve with length penalty(ppclp) algorithm in 3D image

##  Data
```spectExample```: A list, containing data of one 2D "3" number from MNIST dataset, including the x, y coordinates, start point, end point.
```threeExample```: A list, containing data from one 3D SPECT on colon image, including the x, y, z coordinates, start point, end point.

## Example
### 2D MINST example
Load the package into R
```
library(devtools)
install_github("CHuanSite/ppclp")
```
Install the dependent packages
```
library(tidyverse)
library(splines)
library(splines2)
```
Load MNIST data into R from the package **ppclp**
```
data("threeExample")
```
Compute the principal curve 
```
tmpCurve = ppclp2D(threeExample$x, threeExample$y, threeExample$xFix, threeExample$yFix)
```
Plot the final result
```
plot(threeExample$x, threeExample$y, xlim = c(0,1), ylim = c(0,1), pch = 16, cex = 0.8)
lines(tmpCurve$xFit, tmpCurve$yFit, type = "l", col = "red", lwd = 5)
```

### 3D SPECT example
Load the package into R
```
library(devtools)
install_github("CHuanSite/ppclp")
```
Install the dependent packages
```
library(tidyverse)
library(plotly)
library(splines)
library(splines2)
```
Load spect data into R from the package **ppclp**
```
data("spectExample")
```
Compute the principal curve 
```
tmpCurve = ppclp3D(spectExample$x, spectExample$y, spectExample$z, spectExample$xFix, spectExample$yFix, spectExample$zFix)
```
Plot the final result
```
plot_ly() %>% add_trace(x = spectExample$x, y = spectExample$y, z = spectExample$z, type = "scatter3d", mode = "markers", name = 'points', marker = list(size = 1, color = 'rgba(0, 0, 0, .9)', opacity = 0.4)) %>%
add_trace(x = spectExample$xFix[1], y = spectExample$yFix[1], z = spectExample$zFix[1], type = "scatter3d", mode = "markers", name = 'A', marker = list(size = 10, color = 'rgba(0, 255, 0, .9)', opacity = 1)) %>%
add_trace(x = spectExample$xFix[2], y = spectExample$yFix[2], z = spectExample$zFix[2], type = "scatter3d", mode = "markers", name = 'M', marker = list(size = 10, color = 'rgba(0, 0, 255, .9)', opacity = 1)) %>%
add_trace(x = as.vector(tmpCurve$xFit), y = as.vector(tmpCurve$yFit), z = as.vector(tmpCurve$zFit), type = "scatter3d", mode = "lines", name = "theoretical line", line = list(width = 5, color = 'rgba(255, 0, 0, .9)'))
```
