# Danail Obreschkow, Jan 2015 - May 2016

library(astro)
library(pracma)
library(magicaxis)
library(imager)

disk.fit = function(object,
                    vmax= 1000,
                    iterations = 0,
                    show.fit = T,
                    show.details = F) {
  
  # REQUIRES
  # object$data$density: 2D aray with surface densities in arbitrary linear units
  # object$data$velocity: 2D array with velocities in km/s
  # object$data$width
  # object$data$height
  # object$data$psf_in_pixel
  # object$data$kpc_per_pixel
  
  # MAKES
  # object$fit
  # object$angularmomentum
  
  # Coordinate system:
  # bottom left corner has coordinates (0,0)
  # center of the bottom left pixel has coordinates (0.5,0.5)
  # top right corner has coordinates (nx,ny), where c(nx,ny) = dim(object$data$density)
  # centre of top right pixel has coordinates (nx-0.5,ny-0.5)
  
  object = disk.fit.center(object,F)
  object = disk.fit.orientation(object,plot=show.details,vmax=vmax)
  object = disk.fit.rotationcurve(object,plot=show.details,vmax=vmax)
  object = disk.fit.singe_exponential_density_profile(object,plot=F)
  
  if (iterations>0) {
    for (i in seq(iterations)) {
      object = disk.fit.center(object,F,refit=T)
      object = disk.fit.orientation(object,plot=show.details,refit=T)
      object = disk.fit.singe_exponential_density_profile(object,plot=F)
      object = disk.fit.rotationcurve(object,plot=show.details,vmax)
    }
  }
  
  object = disk.fit.double_exponential_density_profile(object,plot=show.details)
  object = disk.angularmomentum(object)
  
  if (show.fit) {plot.fit(object)}
  
  return(object)
  
}

disk.deprojection = function(x,y,x0,y0,i,pa) {
  # x [pixels] = position along the 1st dimension on the sky (e.g. RA)
  # y [pixels] = position along the 2nd dimension on the sky (e.g. decl)
  # x0,y0 [pixels] = center of rotation
  # i [rad] = inclination relative to line of sight in radians
  # pa [rad] = position angle relative to x-axis
  
  xfo = ((x-x0)*cos(pa)+(y-y0)*sin(pa))/(cos(pa)^2+sin(pa)^2)+1e-7 # face-on x-coordinate
  yfo = ((y-y0)*cos(pa)-(x-x0)*sin(pa))/(cos(i)*(cos(pa)^2+sin(pa)^2))+1e-7 # face-on y-coordinate
  r = sqrt(xfo^2+yfo^2) # [pixel] face-on radius from origin
  f = sin(i)/sqrt(1+(yfo/xfo)^2)*sign(xfo) # [-] velocity projection factor: v_radial = f*v_circular
  return(list(radius=r,vfactor=f))
}

disk.fit.center <- function(object,plot=T,densityGamma=2,refit=F) {
  
  # REQUIRES
  # object$data$velocity: 2D array with velocities in km/s
  # object$data$density: 2D aray with surface densities in arbitrary linear units
  
  # MAKES
  # object$fit$deprojection$parameter$x0 [pixel]
  # object$fit$deprojection$parameter$y0 [pixel]
  
  # Coordinate system:
  # bottom left corner has coordinates (0,0)
  # center of the bottom left pixel has coordinates (0.5,0.5)
  # top right corner has coordinates (nx,ny), where c(nx,ny) = dim(object$data$density)
  # centre of top right pixel has coordinates (nx-0.5,ny-0.5)
  
  # initiate
  if (refit) {
    densityField = object$data$density-object$fit$density_model$parameter$offset
    velocityField = object$data$velocity-object$fit$velocity_model$parameter$voffset
  } else {
    densityField = object$data$density
    velocityField = object$data$velocity
  }
  densityPositive = densityField
  densityPositive[densityPositive<0] = 0
  nx = object$data$width
  ny = object$data$height

  # make weighted velocity field
  v = velocityField*densityPositive^densityGamma
  v[is.na(v)] = 0

  # compute centre by minimizing the convolution with the 180deg-rotated image
  superv = array(0,dim(v)*2)
  superv[1:nx,1:ny] = v
  conv = stats::convolve(flipud(fliplr(superv)),superv,type='c')
  index = which.min(conv)
  ix = (index-1)%%(dim(conv)[1])+1
  iy = ceiling(index/(dim(conv)[1]))
  ixreal = ix+sum(c(-1,0,1)%*%conv[(ix-1):(ix+1),(iy-1):(iy+1)])/sum(conv[(ix-1):(ix+1),(iy-1):(iy+1)])
  iyreal = iy+sum(conv[(ix-1):(ix+1),(iy-1):(iy+1)]%*%c(-1,0,1))/sum(conv[(ix-1):(ix+1),(iy-1):(iy+1)])
  bestshift = c((ixreal-nx-1)/2,(iyreal-ny-1)/2)

  # make final output
  object$fit$deprojection$parameter$x0 = nx/2-bestshift[1]
  object$fit$deprojection$parameter$y0 = ny/2-bestshift[2]
  if (!is.null(object$data[["forced_x_offset"]])) {
    object$fit$deprojection$parameter$x0 = object$fit$deprojection$parameter$x0+object$data$forced_x_offset
  }
  if (!is.null(object$data[["forced_y_offset"]])) {
    object$fit$deprojection$parameter$y0 = object$fit$deprojection$parameter$y0+object$data$forced_y_offset
  }
  
  # plot
  if (plot) {
    my.plotNew(xlim=c(0,nx),ylim=c(0,ny),xlab='[pixel]',ylab='[pixel]',title=paste0('Verification of center (',object$name,')'),aspect=1)
    rgb = my.velocityColor(vel=velocityField,velmax=max(velocityField,na.rm=T),lum=densityPositive/max(densityPositive,na.rm=T),gamma=0.4)
    list = is.na(velocityField)
    if (!is.null(object$data[["pixelInsideFrame"]])) {
      object$data$pixelInsideFrame = array(T,c(object$data$width,object$data$height))
    }
    my.plotMatrix(rgb,xleft=0,xright=nx,ybottom=0,ytop=ny,ticks=F,filter=object$data$pixelInsideFrame)
    abline(v=object$fit$deprojection$parameter$x0,col='white')
    abline(h=object$fit$deprojection$parameter$y0,col='white')
  }
  
  return(object)
}

disk.fit.orientation <- function(object,plot=T,
                                 relative_velocity_weight=0.3, # [0..1] weight of velocity relative to (velocity+density) in disk fitting
                                 refit=F,vmax=300) {
  
  # REQUIRES
  # object$data$velocity: 2D array with velocities in km/s
  # object$data$density: 2D aray with surface densities in arbitrary linear units
  # object$fit$deprojection$parameter$x0
  # object$fit$deprojection$parameter$y0
  
  # MAKES
  # object$fit$deprojection$parameter$inclination
  # object$fit$deprojection$parameter$positionangle
  # object$fit$deprojection$field$radius [pixel]: deprojected galactocentric radius
  # object$fit$deprojection$field$velocity_factor = v_line_of_sight/v_circular
  
  # Coordinate system:
  # bottom left corner has coordinates (0,0)
  # center of the bottom left pixel has coordinates (0.5,0.5)
  # top right corner has coordinates (nx,ny), where c(nx,ny) = dim(object$data$density)
  # centre of top right pixel has coordinates (nx-0.5,ny-0.5)
  
  # initiate
  if (refit) {
    densityField = object$data$density-object$fit$density_model$parameter$offset
    velocityField = object$data$velocity-object$fit$velocity_model$parameter$voffset
  } else {
    densityField = object$data$density
    velocityField = object$data$velocity
  }

  isnb_densityField = !is.na(densityField)
  isnb_velocityField = !is.na(velocityField)
  sub_densityField = densityField[isnb_densityField]
  sub_velocityField = velocityField[isnb_velocityField]
  
  # make mesh-grid
  nx = object$data$width
  ny = object$data$height
  xrange = seq(0.5,nx-0.5) # pixel center
  yrange = seq(0.5,ny-0.5) # pixel center
  grid = meshgrid(yrange,xrange)
  x = grid$Y
  y = grid$X
  
  field_norm = c(mean(densityField,na.rm=T),
                 mean(abs(velocityField),na.rm=T)*2) # factor 2 is roughly the average deprojection factor
  
  min_fctA = function(p) {
    
    # initiation
    dep = disk.deprojection(x,y,object$fit$deprojection$parameter$x0,object$fit$deprojection$parameter$y0,p[1],p[2])
    stdev = list(density=0,velocity=0,total=0)
    
    # computation of deviations in regions of similar radius
    for (field_type in seq(2)) {
      if (field_type==1) {
        field = sub_densityField # density field
        radius_list = sort.int(array(dep$radius[isnb_densityField]),index.return = T)$ix
        weight = array(1,length(field))
      } else {
        field = sub_velocityField/pmax(1e-6,abs(dep$vfactor[isnb_velocityField]))*sign(dep$vfactor[isnb_velocityField]) # deprojected velocity field
        radius_list = sort.int(array(dep$radius[isnb_velocityField]),index.return = T)$ix
        weight = abs(dep$vfactor[isnb_velocityField])
      }
      std = 0
      count = 0
      npt = 50
      for (i in seq(1,length(radius_list)-npt,by=npt)) {
        bin_indices = radius_list[i:(i+npt)]
        avg = sum(weight[bin_indices]*field[bin_indices])/sum(weight[bin_indices])
        std = std+sqrt(sum(weight[bin_indices]*(field[bin_indices]-avg)^2)/sum(weight[bin_indices]))
        count = count+1
      }
      stdev[[field_type]] = std/count/field_norm[field_type]
    }
    
    # final result
    stdev$total = (1-relative_velocity_weight)*stdev$density+relative_velocity_weight*stdev$velocity
    return(stdev$total)
  }
  
  # optimization
  val = 1e99
  for (inclination in c(0.1,0.2,0.4,0.8)) {
  #for (inclination in c(0.1,pi/4)) {
    for (pa in c(0,1,2,3)*pi/4) {
      outTmp = optim(c(inclination,pa),min_fctA,control = list(reltol=1e-10))
      if (outTmp$value<val) {
        outA=outTmp
        val=outTmp$value
      }
    }
  }
  if (plot) cat(sprintf('disk.fit.orientation (Goodness of fit = %1.2f)\n',outA$value))
  
  # check if P.A. needs to be reversed
  dep = disk.deprojection(x,y,object$fit$deprojection$parameter$x0,object$fit$deprojection$parameter$y0,outA$par[1],outA$par[2])
  used_pixels = densityField>0
  sign_check = sum(velocityField[used_pixels]*dep$vfactor[used_pixels],na.rm=T)
  if (sum(sign_check)<0) {outA$par[2]=outA$par[2]+pi}
  
  # recast incl and P.A. into desiried intervals
  outA$par[1] = outA$par[1]%%(2*pi)
  if (outA$par[1]>=pi) {
    outA$par[1] = outA$par[1]-pi
    outA$par[2] = outA$par[2]+pi
  }
  if (outA$par[1]>=pi/2) {
    outA$par[1] = pi-outA$par[1]
  }
  outA$par[2] = outA$par[2]%%(2*pi)
  
  if (!is.null(object$data[["forced_inclination"]])) {
    outA$par[1] = object$data$forced_inclination
  }
  
  # make final deprojection matrix
  deprojection = disk.deprojection(x,y,object$fit$deprojection$parameter$x0,object$fit$deprojection$parameter$y0,outA$par[1],outA$par[2])
  
  # make true position angle (= measured from north celestial pole, rather than from horizontal)
  if (is.null(object$data[["rotation"]])) {
    truePA = NULL
  } else {
    truePA = (outA$par[2]+object$data$rotation/180*pi)%%(2*pi)
  }
  
  # create list
  object$fit$deprojection$parameter$inclination=outA$par[1]
  #object$fit$deprojection$parameter$inclination_err=sigma[1]
  object$fit$deprojection$parameter$positionangle=outA$par[2]
  #object$fit$deprojection$parameter$positionangle_err=sigma[2]
  object$fit$deprojection$parameter$positionangle_from_NCP=truePA
  object$fit$deprojection$field$radius=deprojection$radius
  object$fit$deprojection$field$velocity_factor=deprojection$vfactor
  
  # plot
  if (plot) {
    my.plotNew(xlim=c(0,nx),ylim=c(0,ny),xlab='[pixel]',ylab='[pixel]',title=paste0('Verification of center and orientation (',object$name,')'),aspect=1)
    lum = densityField/max(densityField,na.rm=T)
    rgb = my.velocityColor(vel=velocityField,velmax=max(velocityField,na.rm=T),lum=lum,gamma=0.3)
    list = is.na(velocityField) & !is.na(lum)
    rgb[list] = my.RGBarray(lum[list],lum[list],lum[list],Rmin=0,Gmin=0,Bmin=0,Rmax=1.0,Gmax=1.0,Bmax=1.0,gamma=0.3)
    if (!is.null(object$data[["pixelInsideFrame"]])) {
      object$data$pixelInsideFrame = array(T,c(object$data$width,object$data$height))
    }
    my.plotMatrix(rgb,xleft=0,xright=nx,ybottom=0,ytop=ny,ticks=F,filter=object$data$pixelInsideFrame)
    points(object$fit$deprojection$parameter$x0,object$fit$deprojection$parameter$y0,col='white',cex=1,pch=20)
    xrange = c(0,nx)
    lines(xrange,(xrange-object$fit$deprojection$parameter$x0)*tan(object$fit$deprojection$parameter$positionangle)+object$fit$deprojection$parameter$y0,col='white',lwd=1.5,lty=2)
    lines(xrange,(xrange-object$fit$deprojection$parameter$x0)*tan(object$fit$deprojection$parameter$positionangle+pi/2)+object$fit$deprojection$parameter$y0,col='white',lwd=1.5,lty=2)
    radius = min(c(nx,ny))/2
    draw.ellipse(object$fit$deprojection$parameter$x0,object$fit$deprojection$parameter$y0,a=radius,b=radius*cos(object$fit$deprojection$parameter$inclination),angle=object$fit$deprojection$parameter$positionangle/pi*180,
                 border='white',nv=100,col=NULL,lty=2,lwd=1.5,deg =T)
    radius = min(c(nx,ny))/4
    draw.ellipse(object$fit$deprojection$parameter$x0,object$fit$deprojection$parameter$y0,a=radius,b=radius*cos(object$fit$deprojection$parameter$inclination),angle=object$fit$deprojection$parameter$positionangle/pi*180,
                 border='white',nv=100,col=NULL,lty=2,lwd=1.5,deg =T)
  }
  return(object)  
}

disk.fit.rotationcurve <- function(object,plot=T,vmax=300,gammaDensity=0.0) {
  
  # REQUIRES
  # object$data$velocity: 2D array with velocities in km/s
  # object$data$density: 2D aray with surface densities in arbitrary linear units
  # Output from disk.fit.center
  # Output from disk.fit.orientation
  
  # MAKES
  # diverse fits related to the rotation curve, including a radial profile,
  # circular+radial velocity fields, residual field
  
  # Coordinate system:
  # bottom left corner has coordinates (0,0)
  # center of the bottom left pixel has coordinates (0.5,0.5)
  # top right corner has coordinates (nx,ny), where c(nx,ny) = dim(object$data$density)
  # centre of top right pixel has coordinates (nx-0.5,ny-0.5)
  
  # initiate
  densityPositive = object$data$density
  densityPositive[densityPositive<0] = 0
  
  # kinematic model
  v_circular_model = function(r,vflat,rflat) {
    return(vflat*(1-exp(-r/rflat)))
  }
  v_radial_model = function(voffset,sqrtvflat,sqrtrflat) {
    return(voffset+object$fit$deprojection$field$velocity_factor*v_circular_model(object$fit$deprojection$field$radius,sqrtvflat^2,sqrtrflat^2))
  }
  
  # function to minimize
  min_fctB = function(p) {
    vm = v_radial_model(p[1],p[2],p[3])
    if (object$data$psf_in_pixel>0) {
      vm = my.blur(vm,sigma=object$data$psf_in_pixel)
    }
    dv = vm-object$data$velocity
    punishment = (erf(p[2]^2-vmax)+1)/2 # punishment to ensure that vflat cannot be much above vmax
    return(sum((dv*densityPositive^gammaDensity)^2,na.rm=T)*(1+punishment))
  }
  
  # initial conditions for optimization
  vmax_guess = 2.5*quantile(abs(object$data$velocity),0.5,na.rm=T,names=F)
  pinitB = c(0,sqrt(vmax_guess),sqrt(1))
  outB = optim(pinitB,min_fctB,control=list(reltol=1e-15))
  
  if (plot) cat(sprintf('disk.fit.rotationcurve (Goodness of fit = %1.2f)\n',outB$value/6e5))
  quality_rotationcurve <<- outB$value
  
  voffset = outB$par[1]
  
  # deproject measured velocity field after subtraction of voffset
  f = (0.01+abs(object$fit$deprojection$field$velocity_factor))*sign(object$fit$deprojection$field$velocity_factor)
  deprojected_velocityMap=(object$data$velocity-voffset)/f
  deprojected_velocityMap[abs(f)<0.15*max(abs(f)) & !is.na(deprojected_velocityMap)] = NA
  deprojected_velocityMap[deprojected_velocityMap>outB$par[2]^2*2] = NA
  deprojected_velocityMap[deprojected_velocityMap<0] = NA
  
  # create list
  object$fit$deprojection$field$circular_velocity = deprojected_velocityMap
  object$fit$velocity_model$parameter = list(voffset=voffset,
                                             vflat=outB$par[2]^2,
                                             rflat=outB$par[3]^2)

  v_circular_model_best = function(r) {
    # r is the DEPROECTED galactocentric radius in pixels
    return(v_circular_model(r,object$fit$velocity_model$parameter$vflat,object$fit$velocity_model$parameter$rflat))
  }
  vm = v_radial_model(0,outB$par[2],outB$par[3])
  if (object$data$psf_in_pixel>0) {
    vm = my.blur(vm,sigma=object$data$psf_in_pixel)
  }
  object$fit$velocity_model$field = list(circular_velocity=v_circular_model_best(object$fit$deprojection$field$radius),
                                         projected_velocity=v_radial_model(0,outB$par[2],outB$par[3]),
                                         residual_velocity=object$data$velocity-object$fit$velocity_model$parameter$voffset-vm)
  
  # plot fit
  if (plot) {
    rmax = max(object$fit$deprojection$field$radius[!is.na(object$data$velocity)])
    my.plotNew(xlim=c(0,rmax),ylim=c(-2*outB$par[2]^2,2*outB$par[2]^2),xlab='Deprojected radius [pixel]',ylab='Velocity [Units of object$data$velocity]',title=paste0('Verification of rotation curve (',object$name,')'))
    list = !is.na(object$data$velocity)
    threshold = mean(abs(object$fit$deprojection$field$velocity_factor[list]))
    list = list & (abs(object$fit$deprojection$field$velocity_factor)>threshold)
    threshold = mean(abs(object$fit$deprojection$field$velocity_factor[list]))
    list = list & (abs(object$fit$deprojection$field$velocity_factor)>threshold)
    color = my.velocityColor(vel=array(object$data$velocity[list]),velmax=1.2*outB$par[2]^2*sin(object$fit$deprojection$parameter$inclination),lum=1)
    vm = object$data$velocity-object$fit$velocity_model$parameter$voffset
    points(array(object$fit$deprojection$field$radius[list]),array(vm[list]/object$fit$deprojection$field$velocity_factor[list]),pch=20,cex=0.5,col=color)
    r = seq(0,rmax,length.out=200)
    lines(r,v_circular_model(r,outB$par[2]^2,outB$par[3]^2),col='red',lwd=1)
    
    vm = object$fit$velocity_model$field$projected_velocity
    if (object$data$psf_in_pixel>0) {
      vm = my.blur(vm,sigma=object$data$psf_in_pixel)
    }
    vm = vm/object$fit$deprojection$field$velocity_factor
    blurred_model = array(0,ceiling(rmax)+1)
    blurred_model_n = array(0,ceiling(rmax)+1)
    nx = object$data$width
    ny = object$data$height
    for (ix in seq(nx)) {
      for (iy in seq(ny)) {
        if (abs(object$fit$deprojection$field$velocity_factor[ix,iy])>threshold) {
          i = object$fit$deprojection$field$radius[ix,iy]+1
          fi = floor(i)
          ci = ceiling(i)
          f = ci-i
          blurred_model_n[fi] =  blurred_model_n[fi]+f
          blurred_model_n[ci] =  blurred_model_n[ci]+(1-f)
          blurred_model[fi] =  blurred_model[fi]+f*vm[ix,iy]
          blurred_model[ci] =  blurred_model[ci]+(1-f)*vm[ix,iy]
        }
      }
    }
    blurred_model = blurred_model/blurred_model_n
    lines(seq(ceiling(rmax)),blurred_model[2:(ceiling(rmax)+1)],col='red',lwd=3)
  }
  #stop()
  
  return(object)
  
}

disk.fit.singe_exponential_density_profile <- function(object,plot=T,rexp_scale_factor=1.0) {
  
  # ! rexp_scale_factor is a fudge factor to correct the bias of inclination and pa errors
  
  # REQUIRES
  # object$data$velocity: 2D array with velocities in km/s
  # object$data$density: 2D aray with surface densities in arbitrary linear units
  # Output from disk.fit.center
  # Output from disk.fit.orientation
  
  # MAKES
  # diverse fits related to the surface density, including a radial profile,
  # model field, resitual field, morphology estimate, half-mass radius
  
  # Coordinate system:
  # bottom left corner has coordinates (0,0)
  # center of the bottom left pixel has coordinates (0.5,0.5)
  # top right corner has coordinates (nx,ny), where c(nx,ny) = dim(object$data$density)
  # centre of top right pixel has coordinates (nx-0.5,ny-0.5)
  
  # initiate
  densityPositive = object$data$density
  densityPositive[densityPositive<0] = 0
  nx = object$data$width
  ny = object$data$height
  xrange = seq(0.5,nx-0.5) # pixel center
  yrange = seq(0.5,ny-0.5) # pixel center
  grid = meshgrid(yrange,xrange)
  x = grid$Y
  y = grid$X
  
  # surface density model
  sigma_model = function(r,para) {
    # r is the DEPROECTED galactocentric radius in pixels
    sd = para[1] # projected disk normalization
    rd = para[2] # [pixel] disk scale length
    return(sd*exp(-r/max(1e-10,rd)))
  }
  
  # disk optimization
  min_fct1 = function(p) {
    sig = sigma_model(object$fit$deprojection$field$radius,exp(p))
    if (object$data$psf_in_pixel>0) {
      sig = my.blur(sig,sigma=object$data$psf_in_pixel)
    }
    dev = object$data$density-p[3]-sig
    return(sum(dev^2,na.rm=T))
  }
  par1 = c(log(c(max(densityPositive,na.rm=T),(nx+ny)/5)),0)
  out1 = optim(par1,min_fct1,hessian=FALSE)
  par1 = out1$par
  
  min_fct2 = function(p) {
     sig = sigma_model(object$fit$deprojection$field$radius,exp(p))
     if (object$data$psf_in_pixel>0) {
       sig = my.blur(sig,sigma=object$data$psf_in_pixel)
     }
    dev = object$data$density-p[3]-sig
    if (!is.null(object$data[["forced_rexp_fraction"]])) {
      dev[object$fit$deprojection$field$radius<object$data$forced_rexp_fraction*exp(par1[2])] = 0
    } else {
      dev[object$fit$deprojection$field$radius<exp(par1[2])] = 0
    }
    return(sum(dev^2,na.rm=T))
  }
  out2 = optim(par1,min_fct2,hessian=FALSE)
  par2 = out2$par
  
  if (plot) cat(sprintf('disk.fit.singe_exponential_density_profile (Goodness of fit = %1.2f)\n',out2$value))
  quality_exponentialdensityprofile <<- out2$value
  
  # create list
  object$fit$density_model$parameter = list(single_exponential_s=exp(par2[1]),
                                            single_exponential_r=exp(par2[2])*rexp_scale_factor,
                                            offset=par2[3])
  
  # plot fit
  if (plot) {
    list = !is.na(object$data$density)
    rmax = max(object$fit$deprojection$field$radius[list])
    my.plotNew(xlim=c(0,rmax),ylim=c(-0.1,1),xlab='Deprojected radius [pixel]',ylab='Density [Units of object$data$density]',title=paste0('Verification of single exponential density profile (',object$name,')'))
    points(array(object$fit$deprojection$field$radius[list]),array(object$data$density[list]),pch=21,cex=0.5)
    r = seq(0,rmax,length.out=200)
    #lines(r,sigma_model(r,exp(par2[1:2])),col='blue',lwd=1)
    #lines(r,array(par2[3],length(r)),col='black',lwd=1)
    lines(r,sigma_model(r,exp(par2[1:2]))+object$fit$density_model$parameter$offset,col='red',lwd=1)
    
    sig = sigma_model(object$fit$deprojection$field$radius,exp(par2[1:2]))
    if (object$data$psf_in_pixel>0) {
      sig = my.blur(sig,sigma=object$data$psf_in_pixel)
    }
    blurred_model = array(0,ceiling(rmax)+1)
    blurred_model_n = array(0,ceiling(rmax)+1)
    for (ix in seq(nx)) {
      for (iy in seq(ny)) {
        if (!is.na(object$data$density[ix,iy])) {
          i = object$fit$deprojection$field$radius[ix,iy]+1
          fi = floor(i)
          ci = ceiling(i)
          f = ci-i
          blurred_model_n[fi] =  blurred_model_n[fi]+f
          blurred_model_n[ci] =  blurred_model_n[ci]+(1-f)
          blurred_model[fi] =  blurred_model[fi]+f*sig[ix,iy]
          blurred_model[ci] =  blurred_model[ci]+(1-f)*sig[ix,iy]
        }
      }
    }
    blurred_model = blurred_model/blurred_model_n
    lines(seq(ceiling(rmax)),blurred_model[2:(ceiling(rmax)+1)]+object$fit$density_model$parameter$offset,col='red',lwd=3)
  }
  
  return(object)
}

disk.fit.double_exponential_density_profile <- function(object,plot=T) {
  
  # REQUIRES
  # object$data$velocity: 2D array with velocities in km/s
  # object$data$density: 2D aray with surface densities in arbitrary linear units
  # Output from disk.fit.center
  # Output from disk.fit.orientation
  
  # MAKES
  # diverse fits related to the surface density, including a radial profile,
  # model field, resitual field, morphology estimate, half-mass radius
  
  # Coordinate system:
  # bottom left corner has coordinates (0,0)
  # center of the bottom left pixel has coordinates (0.5,0.5)
  # top right corner has coordinates (nx,ny), where c(nx,ny) = dim(object$data$density)
  # centre of top right pixel has coordinates (nx-0.5,ny-0.5)
  
  # initiate  
  densityPositive = object$data$density
  densityPositive[densityPositive<0] = 0
  nx = object$data$width
  ny = object$data$height
  xrange = seq(0.5,nx-0.5) # pixel center
  yrange = seq(0.5,ny-0.5) # pixel center
  grid = meshgrid(yrange,xrange)
  x = grid$Y
  y = grid$X
  
  # surface density model
  sigma_model = function(r,para) {
    # r is the DEPROECTED galactocentric radius in pixels
    # we assume that the bulge is exponential and 'flat' in the sense
    # that it also needs to be deprojected
    sd = para[1] # projected disk normalization
    rd = para[2] # [pixel] disk scale length
    sb = para[3] # projected bluge normalization
    rb = para[4] # [pixel] bluge scale length
    return(sd*exp(-r/max(1e-10,rd))+sb*exp(-r/max(1e-10,rb)))
  }
  
  # disk+bulge optimization
  dm = object$data$density-object$fit$density_model$parameter$offset
  min_fct2 = function(p) {
    sig = sigma_model(object$fit$deprojection$field$radius,exp(p))
    if (object$data$psf_in_pixel>0) {
      sig = my.blur(sig,sigma=object$data$psf_in_pixel)
    }
    dev = dm-sig
    punishment1 = erf((exp(p[2])-1.3*object$fit$density_model$parameter$single_exponential_r)*2)+1 # punishment if disk radius larger than 2-times old disk radius
    punishment2 = erf((p[4]-p[2])*2)+1 # punishment if bulge radius larger than disk radius
    return(sum(dev^2,na.rm=T)*(1+punishment1+punishment2))
  }
  par2 = log(c(object$fit$density_model$parameter$single_exponential_s,
               object$fit$density_model$parameter$single_exponential_r,
               object$fit$density_model$parameter$single_exponential_s,
               object$fit$density_model$parameter$single_exponential_r/2))
  out2 = optim(par2,min_fct2,hessian=FALSE)
  par2 = out2$par
  if (plot) cat(sprintf('disk.fit.densityprofile (Goodness of fit = %1.2f)\n',out2$value/5))
  
  # compute bulge to total ratio
  Md = exp(par2[1])*exp(par2[2])^2
  Mb = exp(par2[3])*exp(par2[4])^2
  B_to_T = Mb/(Md+Mb)
  
  # compute half mass radius
  rmax = max(c(exp(par2[2]),exp(par2[4])))*20
  rlist = seq(0,rmax,by=0.01) # pixel
  cs = cumsum(rlist*sigma_model(rlist,exp(par2)))
  rhalf = (rlist[cs>=max(cs)/2])[1] # pixel
  
  # create list
  object$fit$density_model$parameter = c(object$fit$density_model$parameter,
                                         list(double_exponential_s1=exp(par2[1]),
                                              double_exponential_r1=exp(par2[2]),
                                              double_exponential_s2=exp(par2[3]),
                                              double_exponential_r2=exp(par2[4]),
                                              r_halfmass=rhalf,
                                              B_to_T=B_to_T))
  
  # make density model
  densityModel = sigma_model(object$fit$deprojection$field$radius,exp(par2[1:4]))
  if (object$data$psf_in_pixel>0) {
    sig = my.blur(densityModel,sigma=object$data$psf_in_pixel)
  } else {
    sig = densityModel
  }
  object$fit$density_model$field = list(projected_density_offset_corrected=densityModel,
                                        residual_density=object$data$density-object$fit$density_model$parameter$offset-sig)
  
  # plot fit
  if (plot) {
    list = !is.na(object$data$density)
    rmax = max(object$fit$deprojection$field$radius[list])
    my.plotNew(xlim=c(0,rmax),ylim=c(-0.1,1),xlab='Deprojected radius [pixel]',ylab='Density [Units of object$data$density]',title=paste0('Verification of double exponential density profile (',object$name,')'))
    points(array(object$fit$deprojection$field$radius[list]),array(object$data$density[list]),pch=21,cex=0.5)
    r = seq(0,rmax,length.out=200)
    #lines(r,sigma_model(r,c(exp(par2[1:2]),0,0)),col='blue',lwd=1)
    #lines(r,sigma_model(r,c(0,0,exp(par2[3:4]))),col='red',lwd=1)
    lines(r,sigma_model(r,exp(par2[1:4]))+object$fit$density_model$parameter$offset,col='red',lwd=1)
    
    blurred_model = array(0,ceiling(rmax)+1)
    blurred_model_n = array(0,ceiling(rmax)+1)
    for (ix in seq(nx)) {
      for (iy in seq(ny)) {
        if (!is.na(object$data$density[ix,iy])) {
          i = object$fit$deprojection$field$radius[ix,iy]+1
          fi = floor(i)
          ci = ceiling(i)
          f = ci-i
          blurred_model_n[fi] =  blurred_model_n[fi]+f
          blurred_model_n[ci] =  blurred_model_n[ci]+(1-f)
          blurred_model[fi] =  blurred_model[fi]+f*sig[ix,iy]
          blurred_model[ci] =  blurred_model[ci]+(1-f)*sig[ix,iy]
        }
      }
    }
    blurred_model = blurred_model/blurred_model_n
    lines(seq(ceiling(rmax)),blurred_model[2:(ceiling(rmax)+1)]+object$fit$density_model$parameter$offset,col='red',lwd=3)
  }
  
  return(object)
}

disk.angularmomentum = function(object) {
  
  # REQUIRES
  # objects processed via
  # disk.fit.rotationcurve
  # disk.fit.densityprofile
  
  # MAKES
  # object$angular_momentum$j_observed
  # object$angular_momentum$j_model
  # object$angular_momentum$j_extrapolated (=observed where available, model otherwise)
  
  # j_model
  xd = object$fit$density_model$parameter$single_exponential_r/object$fit$velocity_model$parameter$rflat
  jd = 2*((1+xd)^3-1)/(1+xd)^3*object$fit$velocity_model$parameter$vflat*object$fit$density_model$parameter$single_exponential_r
  j_model = jd*object$data$kpc_per_pixel # [kpc km/s]
  
  # j_observed
  list = !is.na(object$data$density+object$fit$deprojection$field$circular_velocity)
  dm = object$data$density-object$fit$density_model$parameter$offset
  J_observed = sum(object$fit$deprojection$field$radius[list]*dm[list]*object$fit$deprojection$field$circular_velocity[list])
  M_observed = sum(dm[list])
  j_observed = J_observed/M_observed*object$data$kpc_per_pixel # [kpc km/s]
  
  # j_extrapolated (= observed and model, where no obs available)
  J_model_part = sum(object$fit$density_model$parameter$single_exponential_s*exp(-object$fit$deprojection$field$radius[list]/object$fit$density_model$parameter$single_exponential_r)*object$fit$deprojection$field$radius[list]*object$fit$velocity_model$parameter$vflat*(1-exp(-object$fit$deprojection$field$radius[list]/object$fit$velocity_model$parameter$rflat)))
  M_model_part = sum(object$fit$density_model$parameter$single_exponential_s*exp(-object$fit$deprojection$field$radius[list]/object$fit$density_model$parameter$single_exponential_r))
  drortho = cos(object$fit$deprojection$parameter$inclination)
  M_model = drortho*2*pi*(object$fit$density_model$parameter$single_exponential_s*object$fit$density_model$parameter$single_exponential_r^2)
  J_model = jd*M_model
  J_extrapolated = J_model-J_model_part+J_observed
  M_extrapolated = M_model-M_model_part+M_observed
  j_extrapolated = max(j_observed,J_extrapolated/M_extrapolated*object$data$kpc_per_pixel) # [kpc km/s]
  fobs = max(0,min(1,J_observed/max(J_extrapolated,J_model)))
  gobs = max(0,min(1,j_observed/max(j_extrapolated,j_model)))
  j_extrapolated_err = (min(fobs,gobs)*0.03+(1-min(fobs,gobs))*0.2)*j_extrapolated
  
  # the model j out to r=infinity is calculated analytically
  xd = object$fit$density_model$parameter$double_exponential_r1/object$fit$velocity_model$parameter$rflat
  xb = object$fit$density_model$parameter$double_exponential_r2/object$fit$velocity_model$parameter$rflat
  jd = 2*((1+xd)^3-1)/(1+xd)^3*object$fit$velocity_model$parameter$vflat*object$fit$density_model$parameter$double_exponential_r1
  jb = 2*((1+xb)^3-1)/(1+xb)^3*object$fit$velocity_model$parameter$vflat*object$fit$density_model$parameter$double_exponential_r2
  sd = object$fit$density_model$parameter$double_exponential_s1*xd^2
  sb = object$fit$density_model$parameter$double_exponential_s2*xb^2
  j_model = (jd*sd+jb*sb)/(sd+sb)*object$data$kpc_per_pixel # [kpc km/s]
  
  # create list
  object$angularmomentum = list(j_observed=j_observed,
                                j_model=j_model,
                                j_extrapolated=j_extrapolated,
                                j_extrapolated_err=j_extrapolated_err)
  
  return(object)
}

merge.objects = function(master_object,slave_object,slave_weight_density=0.35,plot=F) {
  
  # REQUIRES
  # ..._object$data$density: 2D aray with surface densities in arbitrary linear units
  # ..._object$data$velocity: 2D array with velocities in km/s
  # ..._object$data$rotation [deg]
  # ..._object$fit$deprojection$parameter$x0 [pixel]
  # ..._object$fit$deprojection$parameter$y0 [pixel]
  
  # MAKES
  # new object where density and velocity data is taken from master_object if available
  # and from slave_object otherwise
  
  # NOTES
  # slave_object will be rotated positively (= counterclockwise) by
  # slave_object$data$rotation-master_object$data$rotation degrees
  
  new_object = list()
  new_object$name = master_object$name
  if (!is.null(slave_object$data[["forced_inclination"]])) {
    new_object$data$forced_inclination = slave_object$data$forced_inclination
  }
  if (!is.null(master_object$data[["forced_inclination"]])) {
    new_object$data$forced_inclination = master_object$data$forced_inclination
  }
  new_object$data$psf_in_pixel = mean(master_object$data$psf_in_pixel,slave_object$data$psf_in_pixel)
  if (master_object$data$kpc_per_pixel==slave_object$data$kpc_per_pixel) {
    new_object$data$kpc_per_pixel = master_object$data$kpc_per_pixel
  } else {
    stop('master and slave properties do not match')
  }
  new_object$data$rotation = master_object$data$rotation
  
  alpha = (slave_object$data$rotation-master_object$data$rotation)/180*pi
  co = cos(alpha)
  si = sin(alpha)
  Rs2m = array(c(co,si,-si,co),c(2,2))
  Rm2s = array(c(co,-si,si,co),c(2,2))
  
  s2m = function(xposs,yposs) {
    # converts slave coordinates to master coordinates
    vects = c(xposs-slave_object$fit$deprojection$parameter$x0,
              yposs-slave_object$fit$deprojection$parameter$y0)
    vectm = Rs2m%*%vects
    xposm = vectm[1]+master_object$fit$deprojection$parameter$x0
    yposm = vectm[2]+master_object$fit$deprojection$parameter$y0
    return(list(x=xposm,y=yposm))
  }
  
  m2s = function(xposm,yposm) {
    # converts slave coordinates to master coordinates
    vectm = c(xposm-master_object$fit$deprojection$parameter$x0,
              yposm-master_object$fit$deprojection$parameter$y0)
    vects = Rm2s%*%vectm
    xposs = vects[1]+slave_object$fit$deprojection$parameter$x0
    yposs = vects[2]+slave_object$fit$deprojection$parameter$y0
    return(list(x=xposs,y=yposs))
  }
  
  nxm = dim(master_object$data$density)[1]
  nym = dim(master_object$data$density)[2]
  xrangem = seq(0.5,nxm-0.5) # pixel center
  yrangem = seq(0.5,nym-0.5) # pixel center
  gridm = meshgrid(yrangem,xrangem)
  xm = gridm$Y
  ym = gridm$X
  
  nxs = dim(slave_object$data$density)[1]
  nys = dim(slave_object$data$density)[2]
  xranges = seq(0.5,nxs-0.5) # pixel center
  yranges = seq(0.5,nys-0.5) # pixel center
  grids = meshgrid(yranges,xranges)
  xs = gridm$Y
  ys = gridm$X
  
  # make new array, containing both master and slave
  xmmin = min(xm)
  xmmax = max(xm)
  ymmin = min(ym)
  ymmax = max(ym)
  for (ix in seq(0,1)) {
    for (iy in seq(0,1)) {
      cornerxs = 0.5+ix*(nxs-1) # this is either the min or max of xs
      cornerys = 0.5+iy*(nys-1) # this is either the min or max of ys
      coordm = s2m(cornerxs,cornerys)
      xmmin = min(xmmin,coordm$x)
      xmmax = max(xmmax,coordm$x)
      ymmin = min(ymmin,coordm$y)
      ymmax = max(ymmax,coordm$y)
    }
  }
  ixmmin = ceiling(xmmin)
  ixmmax = ceiling(xmmax)
  iymmin = ceiling(ymmin)
  iymmax = ceiling(ymmax)
  nxnew = ixmmax-ixmmin+1
  nynew = iymmax-iymmin+1
  new_object$data$width = nxnew
  new_object$data$height = nynew
  
  for (mode in seq(4)) { # 1=density map, 2=velocity map, 3=pixelInsideFrame
    if (mode == 1) {
      mapm = master_object$data$density-master_object$fit$density_model$parameter$offset
      maps = slave_weight_density*slave_object$data$density # the slave does not use the offset corrected data, because the density fit and offset are very uncertain due to central beam smearing
    }
    if (mode==2) {
      mapm = master_object$data$velocity-master_object$fit$velocity_model$parameter$voffset
      maps = slave_object$data$velocity-slave_object$fit$velocity_model$parameter$voffset
    }
    if (mode==3) {
      if (!is.null(master_object$data[["pixelInsideFrame"]])) {
        master_object$data$pixelInsideFrame = array(T,c(master_object$data$width,master_object$data$height))
      }
      mapm = master_object$data$pixelInsideFrame
      mapm[mapm] = 1
      mapm[mapm==0] = NA
      if (!is.null(slave_object$data[["pixelInsideFrame"]])) {
        slave_object$data$pixelInsideFrame = array(T,c(slave_object$data$width,slave_object$data$height))
      }
      maps = slave_object$data$pixelInsideFrame
      maps[maps] = 1
      maps[maps==0] = NA
    }
    mapnew = array(NA,c(nxnew,nynew))
    for (ixnew in seq(nxnew)) {
      for (iynew in seq(nynew)) {
        ixm = ixnew+ixmmin-1
        iym = iynew+iymmin-1
        if ((ixm>=1)&(ixm<=nxm)&(iym>=1)&(iym<=nym)) {
          if (!is.na(mapm[ixm,iym])) {
            mapnew[ixnew,iynew]=mapm[ixm,iym]
          }
        }
        if (is.na(mapnew[ixnew,iynew])) {
          coords = m2s(ixm-0.5,iym-0.5)
          ixs = ceiling(coords$x)
          iys = ceiling(coords$y)
          if ((ixs>=1)&(ixs<=nxs)&(iys>=1)&(iys<=nys)) {
            if (!is.na(maps[ixs,iys])) {
              mapnew[ixnew,iynew]=maps[ixs,iys]
            }
          }
        }
      }
    }
    if (mode == 1) new_object$data$density = mapnew
    if (mode == 2) new_object$data$velocity = mapnew
    if (mode == 3) {
      mapnew = as.logical(mapnew)
      mapnew[is.na(mapnew)] = F
      new_object$data$pixelInsideFrame = mapnew
    }
  }
  
  if (plot) {
    densityPositive = new_object$data$density
    densityPositive[densityPositive<0] = 0
    my.plotNew(xlim=c(0,nxnew),ylim=c(0,nynew),xlab='[pixel]',ylab='[pixel]',title=paste0('Verification of superposition (',new_object$name,')'),aspect=1)
    lum = densityPositive/max(densityPositive,na.rm=T)
    rgb = my.velocityColor(vel=new_object$data$velocity,velmax=max(new_object$data$velocity,na.rm=T),lum=lum,gamma=0.3)
    list = is.na(new_object$data$velocity) & !is.na(lum)
    rgb[list] = my.RGBarray(lum[list],lum[list],lum[list],Rmin=0,Gmin=0,Bmin=0,Rmax=1.0,Gmax=1.0,Bmax=1.0,gamma=0.3)
    my.plotMatrix(rgb,xleft=0,xright=nxnew,ybottom=0,ytop=nynew,ticks=F,filter=new_object$data$pixelInsideFrame)
  }
  
  return(new_object)
}

plot.fit <- function(object, residual_scale = 0.25) {
  
  if (!is.null(object$data[["pixelInsideFrame"]])) {
    object$data$pixelInsideFrame = array(T,c(object$data$width,object$data$height))
  }
  
  # make plots
  NAcolor = '#555555'
  posx = c(0,1,0,1)
  posy = c(1,1,0,0)
  w = object$data$width
  h = object$data$height
  vmax_plot = 2.5*quantile(abs(object$data$velocity),0.5,na.rm=T,names=F)
  my.plotNew(xlim=w*c(min(posx),max(posx+1)),ylim=h*c(min(posy),max(posy+1)),title=object$name,
             aspect=T,side=c(),ticks=c(),equal=T,showborder=F)
  
  # velocity field
  lum1 = pmax(0,pmin(1,object$data$density-object$fit$density_model$parameter$offset))
  rgb = my.velocityColor(vel=object$data$velocity,velmax=vmax_plot,lum=lum1,gamma=0.5)
  list = is.na(object$data$velocity-object$fit$velocity_model$parameter$voffset) & !is.na(lum1)
  rgb[list] = my.RGBarray(lum1[list],lum1[list],lum1[list],Rmin=0,Gmin=0,Bmin=0,Rmax=1.0,Gmax=1.0,Bmax=1.0,gamma=0.5)
  my.plotMatrix(rgb,filter=object$data$pixelInsideFrame,NAcolorRGB=NAcolor,
                xleft=posx[1]*w,xright=posx[1]*w+w,ybottom=posy[1]*h,ytop=posy[1]*h+h,ticks=F)
  
  # fitted velocity field
  lum2 = pmax(0,pmin(1,object$fit$density_model$field$projected_density_offset_corrected))
  rgb = my.velocityColor(vel=object$fit$velocity_model$field$projected_velocity,velmax=vmax_plot,lum=lum2,gamma=0.5)
  my.plotMatrix(rgb,filter=object$data$pixelInsideFrame,NAcolorRGB=NAcolor,
                xleft=posx[2]*w,xright=posx[2]*w+w,ybottom=posy[2]*h,ytop=posy[2]*h+h,ticks=F)
  
  # residual velocity field
  rgb = my.velocityColor(vel=object$fit$velocity_model$field$residual_velocity,velmax=residual_scale*vmax_plot,lum=1)
  rgb[is.na(object$fit$velocity_model$field$residual_velocity)] = my.velocityColor(vel=0,lum=1)
  my.plotMatrix(rgb,filter=object$data$pixelInsideFrame,NAcolorRGB=NAcolor,
                xleft=posx[3]*w,xright=posx[3]*w+w,ybottom=posy[3]*h,ytop=posy[3]*h+h,ticks=F)
  
  # residual density field
  rgb = my.velocityColor(vel=object$fit$density_model$field$residual_density,velmax=residual_scale,lum=1,gamma=0.5)
  rgb[is.na(object$fit$density_model$field$residual_density)] = my.velocityColor(vel=0,lum=1)
  my.plotMatrix(rgb,filter=object$data$pixelInsideFrame,NAcolorRGB=NAcolor,
                xleft=posx[4]*w,xright=posx[4]*w+w,ybottom=posy[4]*h,ytop=posy[4]*h+h,ticks=F)
  
  abline(v=posx[4]*w,col='white',lwd=2)
  abline(h=posy[4]*h+h,col='white',lwd=2)
  
  #abline(v=c(posx[1]*w,posx[2]*w+w),col='black',lwd=1.5)
  #abline(h=c(posy[1]*h+h,posy[4]*h),col='black',lwd=1.5)
  
  text(posx[1]*w,posy[1]*h+0.9*h,' Data',pos=4,family='Arial',col='white',cex=1.0)
  text(posx[2]*w,posy[2]*h+0.9*h,' Model',pos=4,family='Arial',col='white',cex=1.0)
  text(posx[3]*w,posy[3]*h+0.9*h,' Velocity residual',pos=4,family='Arial',col='black',cex=1.0)
  text(posx[4]*w,posy[4]*h+0.9*h,' Density residual',pos=4,family='Arial',col='black',cex=1.0)
}

check.deprojection_formula = function() {
  i = 36/180*pi
  alpha = 70/180*pi
  x = 3
  y = 1.84
  d = disk.deprojection(3,1.84,0,0,i,alpha)
  arg = atan2(y,x)
  phi = arg-alpha
  f = sqrt(cos(phi)^2+sin(phi)^2/cos(i)^2)
  g = f/cos(phi)/sin(i)
  print(c(f,d$radius/sqrt(x^2+y^2)))
  print(c(g,1/d$vfactor))
}

my.plotNew <- function(side=c(1,2), ticks=c(1,2,3,4), log='', equal=FALSE, aspect=NULL, xlim=NULL, ylim=NULL, xlab='', ylab='',
                       title=NULL, family='serif',lwd=1, majorn=5, minorn=5,showborder=T,col='black') {
  # Can deal with complicated labels, e.g.
  # xlab="Example:"~bolditalic("V")["max"]~"[km s"^"-1"*"]"
  # ylab="Example: "~frac(sigma,bolditalic("V")["max"])
  if (equal) {par(pty='s')} else {par(pty='m')}
  if (showborder) {setbty='o'} else {setbty='n'}
  plot(x=1,y=NULL,type='n',log=log,xlim=xlim,ylim=ylim,xaxs='i',yaxs='i',xaxt='n',yaxt='n',
       ann=FALSE,asp=aspect,bty=setbty,col.lab=col)
  for (s in seq(4)) {
    if (any(side==s) & any(ticks==s)) {
      magaxis(side=s,xlab=bquote(bold(.(xlab))),ylab=bquote(bold(.(ylab))),family=family,col=col,col.lab=col,lwd=lwd,majorn=majorn,minorn=minorn)
    }
    if (!any(side==s) & any(ticks==s)) {
      magaxis(side=s,labels=FALSE,lwd=lwd,majorn=majorn,minorn=minorn,col=col,col.lab=col)
    }
    if (any(side==s) & !any(ticks==s)) {
      magaxis(side=s,xlab=bquote(bold(.(xlab))),ylab=bquote(bold(.(ylab))),family=family,col=col,col.lab=col,lwd=lwd,majorn=majorn,minorn=minorn,tcl=0)
    }
  }
  title(main=title,family=family,col=col)
  par(family=family,col=col)
}

my.velocityColor <- function(vel,velmax=1,lum=1,lummax=1,gamma=0.65) {
  # vel = 0..velmax
  # lum = 0..lummax
  matrix = (!is.null(dim(lum)))|(!is.null(dim(vel)))
  if (matrix) {
    if (!is.null(dim(lum))) {
      tmpdim = dim(lum)
    } else {
      tmpdim = dim(vel)
    }
  }
  lum = pmax(0,array(lum))
  vel = array(vel)
  if (length(lum)==1) lum = array(lum,length(vel))
  if (length(vel)==1) vel = array(vel,length(lum))
  badData =  is.na(vel+lum)
  vel[badData] = 0
  lum[badData] = 0
  lum = (lum/lummax)^gamma
  col = 1-(pmin(1,pmax(-1,vel/velmax))/2+0.5) # 0..1
  r = (lum*(tanh(2.5-6*col)+1)/2)
  g = lum*(tanh(2-8*abs(col-0.5))+1)/2
  b = (lum*(tanh(-3.5+6*col)+1)/2)^0.8
  g = g+(1-g)*b*0.3
  if (matrix) {
    return(array(rgb(r,g,b),tmpdim))
  } else {
    return (rgb(r,g,b))
  }
}

my.RGBarray <- function(R,G,B,
                        Rmin=min(R,na.rm=T),Gmin=min(G,na.rm=T),Bmin=min(B,na.rm=T),
                        Rmax=max(R,na.rm=T),Gmax=max(G,na.rm=T),Bmax=max(B,na.rm=T),
                        gamma=1,alpha=1,returnNA=F) {
  outputDim = dim(R+G+B)
  Rmax = max(Rmax,Rmin,na.rm=T)
  Gmax = max(Gmax,Gmin,na.rm=T)
  Bmax = max(Bmax,Bmin,na.rm=T)
  R[is.na(R)] = Rmin
  G[is.na(G)] = Gmin
  B[is.na(B)] = Bmin
  gamma[is.na(gamma)] = 1
  alpha[is.na(alpha)] = 1
  R = pmin(1,pmax(0,(R-Rmin)/(Rmax-Rmin+1e-99)))^gamma
  G = pmin(1,pmax(0,(G-Gmin)/(Gmax-Gmin+1e-99)))^gamma
  B = pmin(1,pmax(0,(B-Bmin)/(Bmax-Bmin+1e-99)))^gamma
  if (alpha<1) {
    x = rgb(R,G,B,alpha)
  } else {
    x = rgb(R,G,B)
  }
  dim(x) = outputDim
  if (returnNA) {
    badData =  is.na(R+G+B+gamma+alpha)
    x[badData] = NA
  }
  return(x)
}

my.plotMatrix <- function(M,colorscale=gray.colors(256,start=0,end=1,gamma=1),
                          swapXY=F,reverseX=F,reverseY=F,
                          interpolateRGB=F,filterRGB=NULL,NAcolorRGB='grey',
                          ticks=TRUE,add=TRUE,contour=FALSE,
                          raster=FALSE,lwd=1,lty=1,
                          xleft=NULL,ybottom=NULL,xright=NULL,ytop=NULL,...) {
  
  # M is a 2D matrix M[x,y] real values or colors with x increasing from left to right, and y increasing from bottom to top,
  # the 'colorscale' can be specified
  
  xlogold = par('xlog')
  ylogold = par('ylog')
  usr_old = par('usr')
  par('xlog'=FALSE)
  par('ylog'=FALSE)
  usr = par('usr') # cannot use usr_ol directly, because it messes up with log axis
  if (is.null(xleft)) {xleft=usr[1]}
  if (is.null(xright)) {xright=usr[2]}
  if (is.null(ybottom)) {ybottom=usr[3]}
  if (is.null(ytop)) {ytop=usr[4]}
  if ((!is.numeric(M[1])) & (!is.null(filterRGB))) {M[!filterRGB]=NAcolorRGB}
  if (swapXY) {M=t(M)}
  nx = dim(M)[1]
  ny = dim(M)[2]
  if (reverseX) {M=M[nx:1,]}
  if (reverseY) {M=M[,ny:1]}
  xrange = sort(seq(1/2/nx,1-1/2/nx,length.out=nx)*(xright-xleft)+xleft) # sort because needs to be increasing for image()
  yrange = sort(seq(1/2/ny,1-1/2/ny,length.out=ny)*(ytop-ybottom)+ybottom)
  if (is.numeric(M[1])) {
    if (exists('zlim')) {M=array(pmax(zlim[1]+(zlim[2]-zlim[1])*1e-5,pmin(zlim[2],M)),dim(M))}
    if ((usr[1]<usr[2]) & (usr[3]<usr[4])) {
      par(usr=usr) # needed not to mess up log-plots
    } else if ((usr[1]>usr[2]) & (usr[3]<usr[4])) {
      par(usr=c(usr[2],usr[1],usr[3],usr[4]))
    } else if ((usr[1]<usr[2]) & (usr[3]>usr[4])) {
      par(usr=c(usr[1],usr[2],usr[4],usr[3]))
    } else {
      par(usr=c(usr[2],usr[1],usr[4],usr[3]))
    }
    if (contour) {
      contour(x=xrange,y=yrange,z=M,axes=FALSE,add=add,...)
    } else {
      image(x=xrange,y=yrange,z=M,axes=FALSE,add=add,col=colorscale,xlim=c(xleft,xright),ylim=c(ybottom,ytop))
    }
  } else {
    rasterImage(t(M[,ny:1]),xleft=xleft,ybottom=ybottom,xright=xright,ytop=ytop,interpolate=interpolateRGB)
  }
  if (raster) {
    xlist= seq(1/nx,1-1/nx,length.out=nx-1)
    ylist= seq(1/ny,1-1/ny,length.out=ny-1)
    for (i in seq(length(xlist))) {lines(xlist[i]*c(1,1),c(0,1),lwd=lwd,lty=lty)}
    for (i in seq(length(ylist))) {lines(c(0,1),ylist[i]*c(1,1),lwd=lwd,lty=lty)}
  }
  par('xlog'=xlogold)
  par('ylog'=ylogold)
  par(usr=usr_old)
  if (ticks) {magaxis(side=seq(4),labels=FALSE)}
}

my.hist = function(x,breaks) {
  b = c(-1e99,breaks,1e99)
  counts = hist(x,breaks=b,plot=F)$counts[2:(length(b)-2)]
  center = (breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2
  x = array(rbind(array(breaks),array(breaks)))
  y = c(0,array(rbind(array(counts),array(counts))),0)
  return(list(counts=counts,center=center,x=x,y=y))
}

my.plotNewSub <- function(xleft=0.1,xright=0.3,ybottom=0.1,ytop=0.3) {
  # coordinates relative to full plot
  pomdglobal <<- par(no.readonly = T)#par('omd')
  par(omd=c(0,1,0,1))
  xmarg = sum(par()$mai[c(2,4)])
  xplot = par()$pin[1]
  ymarg = sum(par()$mai[c(1,3)])
  yplot = par()$pin[2]
  #print(c(xleft*xplot/(xplot+xmarg),(xright*xplot+xmarg)/(xplot+xmarg),ybottom*yplot/(yplot+ymarg),(ytop*yplot+ymarg)/(yplot+ymarg)))
  par(new=T,omd=c(xleft*xplot/(xplot+xmarg),(xright*xplot+xmarg)/(xplot+xmarg),ybottom*yplot/(yplot+ymarg),(ytop*yplot+ymarg)/(yplot+ymarg)))
}

my.blur <- function(A,sigma) {
  # A is a 2D matrix to be smoothed by a Gaussian of std sigma
  nx = dim(A)[1]
  ny = dim(A)[2]
  A = array(A,c(nx,ny,1,1))
  A = isoblur(A,sigma)
  A = array(A,c(nx,ny))
  return(A)
}
