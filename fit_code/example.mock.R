# Galaxy parameters
set.seed(3) # random seed for noise of galaxy image
opt_missing_data = T # if true, some data will be missing
psf = 2 # [pixel] standard deviation of point-spread function (=0.4*FWHM)
nx = 42 # [pixel] horizontal size of image
ny = 40 # [pixel] vertical size of image
noise = 0.6
x0 = 25 # [pixel] x-coordinate of galaxy centre
y0 = 20 # [pixel] y-coordinate of galaxy centre
pa = 50/180*pi # [rad] position angle
inclination = 55/180*pi # [rad] inclination
voffset = -1 # [km/s] velocity offset
rflat = 10 # [pixel] kinematic radius
vflat = 170 # [km/s] asymptotic rotation velocity
soffset = 0.02
rdisk = 7.4 # [pixel] exponential scale radius disk
fbulge = 0.4 # [-] mass fraction of bulge
rbulge = 2 # [pixel] exponential scale radius bulge
kpc_per_pixel = 0.25

# Compute exact specific angular momentum of true galaxy
xd = rdisk/rflat
xb = rbulge/rflat
jd = 2*((1+xd)^3-1)/(1+xd)^3*vflat*rdisk
jb = 2*((1+xb)^3-1)/(1+xb)^3*vflat*rbulge
sd = xd^2
sb = fbulge*xb^2
jref = (jd*sd+jb*sb)/(sd+sb)*kpc_per_pixel # [kpc km/s]

# Make mock galaxy

# initiate galaxy
object = list()
object$name = 'Mock Galaxy'
object$data$height = ny
object$data$width = nx
object$data$psf_in_pixel = psf
object$data$kpc_per_pixel = kpc_per_pixel

# Make density field and velocity field
xrange = seq(0.5,nx-0.5) # pixel center
yrange = seq(0.5,ny-0.5) # pixel center
grid = meshgrid(yrange,xrange)
x = grid$Y
y = grid$X
dep = disk.deprojection(x,y,x0,y0,inclination,pa)
object$data$density = exp(-dep$radius/rdisk)+fbulge*exp(-dep$radius/rbulge)+soffset
object$data$velocity = vflat*(1-exp(-dep$radius/rflat))*dep$vfactor+voffset

# Apply PSF blurring
object$data$density = my.blur(object$data$density,sigma=object$data$psf_in_pixel)
object$data$velocity = my.blur(object$data$velocity,sigma=object$data$psf_in_pixel)

# Make data noisy
if (opt_missing_data) {
  object$data$velocity[(runif(nx*ny)>0.8)&((abs(x-x0)>7)|(abs(y-y0)>7))] = NA
  object$data$velocity[object$data$density<0.08] = NA
  object$data$density[1:10,1:5] = NA
  object$data$density[nx-5,] = NA
  object$data$density[10,(ny/2):ny] = NA
  object$data$density[object$data$density<0.05] = NA
  object$data$velocity[is.na(object$data$density)] = NA
}
object$data$velocity = object$data$velocity+array(rnorm(nx*ny),c(nx,ny))*vflat*0.04*noise
object$data$density = object$data$density+array(rnorm(nx*ny),c(nx,ny))*0.04*noise

# fit model
object = disk.fit(object, show.fit = T, show.details = T, iterations = 0)

# output
cat(sprintf('Center x = %4.1f pixel (%4.1f/%4.1f)\n',object$fit$deprojection$parameter$x0-x0,object$fit$deprojection$parameter$x0,x0))
cat(sprintf('Center y = %4.1f pixel (%4.1f/%4.1f)\n',object$fit$deprojection$parameter$y0-y0,object$fit$deprojection$parameter$y0,y0))
cat(sprintf('Incl.    = %4.1f deg (%4.1f/%4.1f)\n',(object$fit$deprojection$parameter$inclination-inclination)/pi*180,object$fit$deprojection$parameter$inclination/pi*180,inclination/pi*180))
cat(sprintf('PA       = %4.1f deg (%4.1f/%4.1f)\n',Arg(exp((object$fit$deprojection$parameter$positionangle-pa)*1i))/pi*180,object$fit$deprojection$parameter$positionangle/pi*180,pa/pi*180))
cat(sprintf('Voffset  = %4.1f km/s (%3.1f/%3.1f)\n',object$data$velocity_offset-voffset,object$data$velocity_offset,voffset))
cat(sprintf('Vflat    = %4.1f km/s (%4.1f/%4.1f)\n',object$fit$velocity_model$parameter$vflat-vflat,object$fit$velocity_model$parameter$vflat,vflat))
cat(sprintf('Rflat    = %4.1f pixel (%4.1f/%4.1f)\n',object$fit$velocity_model$parameter$rflat-rflat,object$fit$velocity_model$parameter$rflat,rflat))
cat(sprintf('Rdisk    = %4.1f pixel (%4.1f/%4.1f)\n',object$fit$density_model$parameter$double_exponential_r1-rdisk,object$fit$density_model$parameter$double_exponential_r1,rdisk))
cat(sprintf('Rbulge   = %4.1f pixel (%4.1f/%4.1f)\n',object$fit$density_model$parameter$double_exponential_r2-rbulge,object$fit$density_model$parameter$double_exponential_r2,rbulge))
cat(sprintf('Soffset  = %4.1f (%4.2f/%4.2f)\n',object$data$density_offset-soffset,object$data$density_offset,soffset))
cat(sprintf('Scenter  = %4.1f (%4.2f/%4.2f)\n',object$fit$density_model$parameter$single_exponential_s-1,object$fit$density_model$parameter$single_exponential_s,1))
cat(sprintf('j        = %4.1f percent (%3.0f/%3.0f)\n',(object$angularmomentum$j_extrapolated-jref)/jref*100,object$angularmomentum$j_extrapolated,jref))