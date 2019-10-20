library(oro.nifti)
library(brainR)
library(rgl)

Img = readANALYZE("/Users/huanchen/Downloads/DREAM-Amer_Sabrine/J001_Product\ A/J001_111516_2hr/1.2.840.113704.3.1.32161115.10110.83/J001_111516_2\ hr_ECTHd1_IRAC001_DS.img")
ctImg = readNIfTI("/Users/huanchen/Downloads/DREAM-Amer_Sabrine/J001_Product\ A/J001_111516_2hr/1.2.840.113704.3.1.32161115.10110.83/J001_111516_2\ hr_ECTHd1_IRAC001_DS.nii.gz")


nColLevels = 11
col = c(gray(0), heat.colors(nColLevels - 1))
colCT = gray(seq(0, 1, length = nColLevels))

withinThresh = Img[Img >= 10]
withinThreshCT = ctImg[ctImg >= 10]
breaks = quantile(withinThresh, (0 : nColLevels) / nColLevels)
breaksCT = c(0, quantile(withinThreshCT, (0 : (nColLevels - 1)) / (nColLevels - 1)))

XYZ = c(70, 70, 32)

image(z = ctImg[,,XYZ[3]],
      1 : dim(Img)[1],
      1 : dim(Img)[2],
      col = colCT,
      asp = 1,
      axes = FALSE,
      #breaks = breaksCT,
      xlab= "",
      ylab = ""
)


image(z = Img[, ,XYZ[3]],
      1 : dim(Img)[1],
      1 : dim(Img)[2],
      col = col,
      breaks = breaks,
      asp = 1,
      axes = FALSE,
      add = TRUE
)

title("Z slice")
abline(v = XYZ[1], col = "white", lwd = 3)
abline(h = XYZ[2], col = "white", lwd = 3)


spectExample = list(x = spectExample$x, y = spectExample$y, z = spectExample$z, xFix = spectExample$xFix, yFix = spectExample$yFix, zFix = spectExample$zFix, Img = Img, ctImg = ctImg)

save(spectExample, file = "data/spectExample.RData")

