dynamo.load = function(index) {
  
  # REQUIRES:
  # index: numerical index of galaxy to be loaded
  
  # MAKES:
  # object$name
  # object$data$density [any normalization] 2D density field, NA if not observed
  # object$data$velocity [km/s] 2D velocity field, NA if not observed
  # object$data$lineIntensity [any normalization] strength of line used for velocity map, NA if not observed
  # object$data$pixelInsideFrame T/F indicates whether pixel is inside observation
  # object$data$rotation [deg] clockwise from north right, east up
  # object$data$width [pixel]
  # object$data$height [pixel]
  # object$data$kpc_per_pixel
  
  # OBJECT OVERVIEW:
  # 001: Dynamo C00_1 @ OSIRIS
  # 002: Dynamo C22_2 @ OSIRIS
  # 003: Dynamo D13_5 @ OSIRIS
  # 004: Dynamo G04_1 @ OSIRIS
  # 005: Dynamo G20_2 @ OSIRIS
  # 006: Dynamo G04_1 @ GMOS
  # 007: Dynamo G20_2 @ GMOS
  # 008: THINGS NGC 3198
  # 009: Dynamo H10_2 @ OSIRIS
  # 010: Dynamo H10_2 @ GMOS
  
  H0 = 70 # km/s/Mpc
  
  object = list()
  
  master_path = '/Users/gsavorgnan/angularmomentum/DanailData/'
  identifier = c('C001','C222','D135','G041','G202','G041','G202','ThingsN3198','H102','H102')[index]
  redshift = c(0.06083,0.07116,0.07535,0.12981,0.14113,0.12981,0.14113,0.14907,0.14907)[index]
  
  # make object name
  object$name = identifier
  
  if ((index<=5)|(index==9)) { # OSIRIS data
    
    path = paste0(master_path,'OSIRIS/')
    
    # make filenames
    filename = {}
    filename[1] = paste0(identifier,'/lum/',identifier,'_Pa_lum.fits')
    filename[2] = paste0(identifier,'/lum/',identifier,'_continuum_lum.fits')
    filename[3] = paste0(identifier,'/kine/',identifier,'_Pa_velocity.fits')
    
    # read maps
    rawmap = list(Pa=array(),continuum=array(),velocity=array())
    for (i in seq(length(rawmap))) {
      fn = file.path(path,filename[i])
      if (file.exists(fn)) {
        fits = read.fits(fn)
        rawmap[[i]] = t(fits$dat[[1]]) # transpose image
        rawmap[[i]] = rawmap[[i]][,dim(rawmap[[i]])[2]:1] # revert y-axis
        # interpret header
        id = grep('axes rotated by',fits$hdr[[1]][,3])[1]
        object$data$rotation = as.numeric(substr(fits$hdr[[1]][id,3],25,31)) # [deg] clockwise from north right, east up
        id = grep('Pixel scale is',fits$hdr[[1]][,3])[1]
        arcsec_per_pixel = as.numeric(substr(fits$hdr[[1]][id,3],16,23))
        Da = 2.998e5*redshift/H0/(1+redshift) # [Mpc] approximate angular diameter distance
        object$data$kpc_per_pixel = arcsec_per_pixel/3600/180*pi*Da*1000
        object$data$psf_in_pixel = 2.4*0.4
      }
    }
  
    # extract dinemsions of map
    object$data$width = dim(rawmap[[1]])[1] # [pixel]
    object$data$height = dim(rawmap[[1]])[2] # [pixel]
    
    continuumDataExists = length(rawmap$continuum)>1
    
    # make matrix showing which pixels are inside the observation
    if (continuumDataExists) {
      object$data$pixelInsideFrame = rawmap$continuum!=0
    } else {
      object$data$pixelInsideFrame = rawmap$Pa>-1000
    }
    
    # make normalised surface density map
    if (continuumDataExists) {
      object$data$density = rawmap$continuum/max(rawmap$continuum,na.rm=T)
      object$data$density[rawmap$continuum==0] = NA
    } else {
      object$data$density = rawmap$Pa/max(rawmap$Pa,na.rm=T)
      object$data$density[(rawmap$Pa==-1000)|is.na(rawmap$Pa)] = NA
    }
    if (index==3) {
      object$data$density[object$data$density<(-0.15)] = NA
    }
    
    # make velocity maps
    object$data$velocity = rawmap$velocity
    object$data$velocity[(rawmap$Pa==-1000)|is.na(rawmap$Pa)] = NA
    vmax = 250 # [km/s]
    object$data$velocity[abs(object$data$velocity)>vmax] = NA # exclude velocity points that seem outside the galaxy velocity
    if (index==4) object$data$velocity[object$data$velocity>100] = NA
    
    # make lineIntensity
    object$data$lineIntensity = rawmap$Pa/max(rawmap$Pa,na.rm=T)
    object$data$lineIntensity[is.na(object$data$velocity)] = NA
    
  }
  
  if ((index==6)|(index==7)) { # GMOS data
    
    path = paste0(master_path,'GMOS/')
    
    # make filenames
    filename = {}
    filename[1] = paste0(identifier,'/',identifier,'_hb.fits')
    filename[2] = paste0(identifier,'/',identifier,'_ct.fits')
    filename[3] = paste0(identifier,'/',identifier,'_vm.fits')
    
    # read maps
    rawmap <<- list(Hbeta=array(),continuum=array(),velocity=array())
    for (i in seq(length(rawmap))) {
      fn = file.path(path,filename[i])
      fits <<- read.fits(fn)
      rawmap[[i]] <<- t(fits$dat[[1]]) # transpose image
      rawmap[[i]] <<- rawmap[[i]][,dim(rawmap[[i]])[2]:1] # revert y-axis
    }
    
    # make some basic parameters
    if (index==6) object$data$rotation = 0 # [deg] clockwise from north right, east up
    if (index==7) object$data$rotation = 225 # [deg] clockwise from north right, east up
    arcsec_per_pixel = 0.1 # [arcsec/pixel]
    object$data$psf_in_pixel = 5*0.4
    Da = 2.998e5*redshift/H0/(1+redshift) # [Mpc] approximate angular diameter distance
    object$data$kpc_per_pixel = arcsec_per_pixel/3600/180*pi*Da*1000
    object$data$width = dim(rawmap[[1]])[1] # [pixel]
    object$data$height = dim(rawmap[[1]])[2] # [pixel]
    
    # make matrix showing which pixels are inside the observation
    object$data$pixelInsideFrame = rawmap$Hbeta>(-100)
    
    # make normalised surface density map
    object$data$density = rawmap$continuum/max(rawmap$continuum,na.rm=T)
    object$data$density[rawmap$continuum<=(-20)] = NA
    object$data$density[rawmap$Hbeta==(-100)] = NA
    
    # make velocity maps
    object$data$velocity = rawmap$velocity#*(1+redshift)
    object$data$velocity[rawmap$velocity<=(-200)] = NA
    object$data$velocity[rawmap$Hbeta==(-100)] = NA
    
    # make lineIntensity
    object$data$lineIntensity = array(pmax(0,rawmap$Hbeta/max(rawmap$Hbeta,na.rm=T)),dim(rawmap$Hbeta))
    object$data$lineIntensity[rawmap$Hbeta==(-100)] = NA
    
  }
  
  if (index==10) { # GMOS data
    
    path = paste0(master_path,'GMOS/')
    
    # make filenames
    filename = {}
    filename[1] = paste0(identifier,'/greg_hb.fits')
    filename[2] = paste0(identifier,'/greg_hb.fits') # ct
    filename[3] = paste0(identifier,'/greg_vm.fits')
    
    # read maps
    rawmap <<- list(Hbeta=array(),continuum=array(),velocity=array())
    for (i in seq(length(rawmap))) {
      fn = file.path(path,filename[i])
      fits <<- read.fits(fn)
      rawmap[[i]] <<- t(fits$dat[[1]]) # transpose image
      rawmap[[i]] <<- rawmap[[i]][,dim(rawmap[[i]])[2]:1] # revert y-axis
    }
    
    # make some basic parameters
    object$data$rotation = 270 # [deg] clockwise from north right, east up
    arcsec_per_pixel = 0.1 # [arcsec/pixel]
    object$data$psf_in_pixel = 5*0.4
    Da = 2.998e5*redshift/H0/(1+redshift) # [Mpc] approximate angular diameter distance
    object$data$kpc_per_pixel = arcsec_per_pixel/3600/180*pi*Da*1000
    object$data$width = dim(rawmap[[1]])[1] # [pixel]
    object$data$height = dim(rawmap[[1]])[2] # [pixel]
    
    # make matrix showing which pixels are inside the observation
    #object$data$pixelInsideFrame = rawmap$Hbeta>(-100)
    
    # make normalised surface density map
    object$data$density = rawmap$continuum/max(rawmap$continuum,na.rm=T)
    object$data$density[rawmap$continuum<=(-200)] = NA
    #object$data$density[rawmap$Hbeta==(-100)] = NA
    
    # make velocity maps
    object$data$velocity = rawmap$velocity
    #object$data$velocity[rawmap$velocity<=(-200)] = NA
    #if (index==10) {
    #  object$data$velocity[rawmap$velocity==(-100)] = NA
    #  object$data$velocity[1:40,] = NA
    #}
    #object$data$velocity[rawmap$Hbeta==(-100)] = NA
    
    # make lineIntensity
    object$data$lineIntensity = rawmap$Hbeta/max(rawmap$Hbeta,na.rm=T)
    #object$data$lineIntensity = array(pmax(0,rawmap$Hbeta/max(rawmap$Hbeta,na.rm=T)),dim(rawmap$Hbeta))
    #object$data$lineIntensity[rawmap$Hbeta==(-100)] = NA
    
  }
  
  if (index==8) { # THINGS data
    
    path = '/Users/do/Documents/Science/Projects/Observational Data/THINGS/'
    
    # make filenames
    filename = {}
    filename[1] = '/mom0/NGC_3198_NA_MOM0_THINGS.FITS'
    filename[2] = 'doc/StellarSurfaceDensity_NGC3198.img'
    filename[3] = '/mom1/NGC_3198_NA_MOM1_THINGS.FITS'
    
    # read maps
    rawmap = list(HI=array(),continuum=array(),velocity=array())
    smallrawmap = list(HI=array(),continuum=array(),velocity=array())
    for (i in seq(length(rawmap))) {
      fn = file.path(path,filename[i])
      if (i==1) {
        fits = read.fits(fn)
        rawmap[[i]] = array(fits$dat[[1]],dim(fits$dat[[1]])[1:2])
      }
      if (i==2) {
        load(file = fn)
        rawmap[[i]] = StellarSurfaceDensity_NGC3198
      }
      if (i==3) {
        fits = read.fits(fn)
        rawmap[[3]] = array(fits$dat[[1]],dim(fits$dat[[1]])[1:2])
        list = !is.na(rawmap[[3]])
        mean_velocity = sum(rawmap[[3]][list]*rawmap[[1]][list])/sum(rawmap[[1]][list])
        rawmap[[3]][list] = (rawmap[[3]][list]-mean_velocity)/1000 # [km/s]
      }
      #smallrawmap[[i]] = rawmap[[i]][257:768,257:768]
      dx = 128-68
      dy = 120-64
      smallrawmap[[i]] = rawmap[[i]][(449-dx):(576+dx),(449-dy):(576+dy)]
    }
    rawmap = smallrawmap
    
    # reduce size
    niterations = 2
    for (iterations in seq(niterations)) {
      nx = dim(rawmap[[1]])[1]
      ny = dim(rawmap[[1]])[2]
      for (i in seq(3)) {
        rawmap[[i]] = 0.25*(rawmap[[i]][seq(1,nx,by=2),seq(1,ny,by=2)]+
                            rawmap[[i]][seq(1,nx,by=2),seq(2,ny,by=2)]+
                            rawmap[[i]][seq(2,nx,by=2),seq(1,ny,by=2)]+
                            rawmap[[i]][seq(2,nx,by=2),seq(2,ny,by=2)])
      }
    }
    
    nx = dim(rawmap[[1]])[1]
    ny = dim(rawmap[[1]])[2]
    for (i in seq(3)) {
      rawmap[[i]] = as.matrix(blur(as.im(rawmap[[i]]),sigma=1.5))
      rawmap[[i]] = rawmap[[i]][,3:ny-2]
    }
    
    # make matrix showing which pixels are inside the observation
    nx = dim(rawmap[[1]])[1]
    ny = dim(rawmap[[1]])[2]
    object$data$pixelInsideFrame = array(T,dim(rawmap[[1]]))
    for (ix in seq(nx)) {
      for (iy in seq(ny)) {
        if (iy>(-ix*5+280)) {object$data$pixelInsideFrame[ix,iy] = F}
        if (iy<(-ix*5+90)) {object$data$pixelInsideFrame[ix,iy] = F}
        if (iy>(ix*0.2+40)) {object$data$pixelInsideFrame[ix,iy] = F}
        if (iy<(ix*0.2+5)) {object$data$pixelInsideFrame[ix,iy] = F}
      }
    }
    
    # cut data
    nx = dim(rawmap[[1]])[1]
    ny = dim(rawmap[[1]])[2]
    set.seed(1)
    list = rawmap[[2]]<0.04*array(runif(nx*ny),c(nx,ny))*max(rawmap[[2]])
    for (i in seq(3)) {
      rawmap[[i]][list] = NA
      rawmap[[i]][!object$data$pixelInsideFrame] = NA
    }
    
    # make some basic parameters
    object$data$rotation = 270 # [deg] clockwise from north right, east up
    arcsec_per_pixel = 1.5*2^niterations # [arcsec/pixel]
    object$data$psf_in_pixel = 0
    distance = 13.8 # [Mpc] from Walter 2008
    object$data$kpc_per_pixel = arcsec_per_pixel/3600/180*pi*distance*1000
    object$data$width = dim(rawmap[[1]])[1] # [pixel]
    object$data$height = dim(rawmap[[1]])[2] # [pixel]
  
    
    # make normalised surface density map
    object$data$density = rawmap$continuum/max(rawmap$continuum,na.rm=T)
    
    # make velocity maps
    object$data$velocity = rawmap$velocity
    
    # make lineIntensity
    object$data$lineIntensity = array(pmax(0,rawmap$HI/max(rawmap$HI,na.rm=T)),dim(rawmap$HI))
    
  }
  
  if (index==2) {
    object$data$forced_x_offset = -1
    object$data$forced_y_offset = +0.0
  }
  if (index==3) {
    object$data$forced_x_offset = 0
    object$data$forced_y_offset = 0
    object$data$velocity = object$data$velocity-40
    object$data$velocity[1:5,1:25] = NA
    object$data$forced_rexp_fraction = 0.5
  }
  if (index==4) {
    object$data$forced_x_offset = 1
    object$data$forced_y_offset = 0
    #object$data$forced_inclination = 40/180*pi
  }
  
  return(object)
}