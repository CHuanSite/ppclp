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
library(ppclp)
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
