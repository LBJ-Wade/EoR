mkdir Images/
#cw Added the following mkdir so it works regardless of their prior existence
mkdir Images/cs_psf
mkdir Images/fg_psf
mkdir Images/noise_psf
mkdir Images/cs_images
mkdir Images/fg_images
mkdir Images/noise_images

mv -f IMG_*/*cs*.psf.fits    Images/cs_psf/.
mv -f IMG_*/*fg*.psf.fits    Images/fg_psf/.
mv -f IMG_*/*noise*.psf.fits Images/noise_psf/.

mv -f IMG_*/*cs*.image.fits    Images/cs_images/.
mv -f IMG_*/*fg*.image.fits    Images/fg_images/.
mv -f IMG_*/*noise*.image.fits Images/noise_images/.
mv -f IMG_*/exportfits* Images/. #cw addition
rm -rf IMG_*
