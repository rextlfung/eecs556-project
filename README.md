# eecs556-project
# Progress Report
March 27, 2024  
We attempted to reproduce the results in the original paper by Jacobs et al. using the soruce code linked but could not get it to run as there were missing functions in the package. We emailed the original author but have yet to receive a reply. As a backup plan we decided to use Professor Jeff Fessler's NUFFT implementation in his Michigan Image Reconstruction Toolbox (MIRT), which performs min-max interpolation using a Kaiser-Bessel kernel. As preliminary experiments, we retrospectively undersampled non-cartesian k-space data generated by taking the Fourier transform of Matlab's built-in 'mristack' dataset. We tested MIRT's NUFFT on radial and spiral sampled data, which produced images with aliasing artifacts similar to when naively reconstructing via zero-insertion and IFFT. We believe that these results makes sense as the purpose of NUFFT is to reconstruct a cartesian image from non-cartesian k-space data, thus it should have little to no improvement in image quality when applied on undersampled data. That being said, there is perhaps a need for a better way to define fully-sampled vs undersampled data in non-cartesian sampling schemes...  
  
For our next steps, we plan to apply the MIRT NUFFT reconstruction to the data in the paper by Jacobs et al. and quantitatively evaluate the reconstructed images using the following metrics: NRMSE, PSNR, SSIM. We think that it will be interesting to see if the images will be better quality than those reported in the paper.  

Ideas for what to do next:  
- Implement NUFFT from scratch using a simple interpolation kernel described in the course material. It will likely perform worse but I think it'll be interesting to see how far off it is.
