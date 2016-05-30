# load, process, fit, and plot galaxies

#galaxy = 'G041'
#galaxy = 'C222'
#galaxy = 'D135'
#galaxy = 'G202'
#galaxy = 'SDSS013527'
#galaxy = 'SDSS234657'


main = function(galaxy, id, osiris, gmos) {

if (osiris &! gmos) {observation = 'osiris'}
if (gmos &! osiris) {observation = 'gmos'}
if (osiris & gmos) {observation = 'osiris+gmos'}

slave_weight_density_value = 1.
if (galaxy=='G041') {slave_weight_density_value = 0.35} 
if (galaxy=='C222') {slave_weight_density_value = 0.9} #0.9} 
if (galaxy=='D135') {slave_weight_density_value = 0.5} 
if (galaxy=='SDSS013527') {slave_weight_density_value = 0.4} 
if (galaxy=='SDSS234657') {slave_weight_density_value = 0.5} 

rd_fixed = NULL # fixed disk scale length, in pixels
if (galaxy=='G041') {rd_fixed = 9.7} 
if (galaxy=='G202') {rd_fixed = 4.1} 
if (galaxy=='D135') {rd_fixed = 20.1}#17.0} 

forced_inclination = NULL
if (galaxy=='G041') {forced_inclination = 35} 
if (galaxy=='G202') {forced_inclination = 46} #25} 
if (galaxy=='D135') {forced_inclination = 38}#57} 



if (osiris & gmos) {
  osiris_object = dynamo.load(galaxy, 'osiris') 
  gmos_object = dynamo.load(galaxy, 'gmos') 
  osiris_object = disk.fit(osiris_object, rd_fixed, show.fit = F, show.details = F)
  gmos_object = disk.fit(gmos_object, rd_fixed, show.fit = F, show.details = F)
  
  # merge osiris and gmos: treats first as master and second as slave
  object = merge.objects(osiris_object,gmos_object,slave_weight_density=slave_weight_density_value)
} else {
  if (osiris) {
    object = dynamo.load(galaxy, 'osiris') 
  } else if (gmos) {
    object = dynamo.load(galaxy, 'gmos')
  } else {
    object = NULL
    print('No data supplied.')
  }
}

if (galaxy=='C222' & osiris & gmos) {
  object$data$forced_x_offset = -4.7
  object$data$forced_y_offset = +0
} else if (galaxy=='D135'& osiris & gmos) {
  object$data$forced_x_offset = +0.5
  object$data$forced_y_offset = -0.5
} else if (galaxy=='SDSS013527'& osiris & gmos) {
  object$data$forced_x_offset = +0.5
  object$data$forced_y_offset = -0.5
} else if (galaxy=='SDSS234657'& osiris & gmos) {
  object$data$forced_x_offset = -0.5
  object$data$forced_y_offset = +0.5
}

if (!is.null(forced_inclination)) {object$data$forced_inclination = forced_inclination/180*pi }
object = disk.fit(object, rd_fixed, show.fit = T, show.details = T, iterations = 0, vmax=250) # fit disk


# print text
cat('===============================================================\n')
cat(sprintf('Galaxy            = %s\n',galaxy))
cat(sprintf('Observation       = %s\n',observation))
if (is.null(object$data$forced_inclination)) { 
  cat(sprintf('inclination       = %6.1f deg\n',object$fit$deprojection$parameter$inclination*180/pi)) 
} else {
  cat(sprintf('inclination       = %6.1f deg [forced]\n',object$fit$deprojection$parameter$inclination*180/pi)) 
}
cat(sprintf('PA from north     = %6.1f deg\n',object$fit$deprojection$parameter$positionangle_from_NCP*180/pi))
cat(sprintf('rdisk             = %6.3f kpc\n',object$fit$density_model$parameter$single_exponential_r*object$data$kpc_per_pixel))
cat(sprintf('rdisk             = %6.3f arcsec\n',object$fit$density_model$parameter$single_exponential_r*0.1))
cat(sprintf('rflat             = %6.1f kpc\n',object$fit$velocity_model$parameter$rflat*object$data$kpc_per_pixel))
cat(sprintf('vflat             = %6.1f km/s\n',round(object$fit$velocity_model$parameter$vflat)))
cat(sprintf('j_observed        = %6.1f kpc km/s\n',object$angularmomentum$j_observed))
cat(sprintf('j_model           = %6.1f kpc km/s\n',object$angularmomentum$j_model))
cat(sprintf('j_total           = %6.1f kpc km/s\n',object$angularmomentum$j_extrapolated))
cat(sprintf('j_total_err       = %6.1f kpc km/s\n',object$angularmomentum$j_extrapolated_err))
cat(sprintf('extrapol fraction = %6.1f%s\n',(object$angularmomentum$j_extrapolated-object$angularmomentum$j_observed)/object$angularmomentum$j_extrapolated*100,'%'))

result_row = c(id, galaxy, observation, (osiris&!gmos), (gmos&!osiris), (osiris&gmos), 
               format(round(object$fit$deprojection$parameter$inclination*180/pi, 2), nsmall = 2), 
               !is.null(object$data$forced_inclination), 
               format(round(object$fit$density_model$parameter$single_exponential_r*object$data$kpc_per_pixel, 2), nsmall = 2), 
               format(round(object$fit$velocity_model$parameter$rflat*object$data$kpc_per_pixel, 2), nsmall = 2),
               round(object$fit$velocity_model$parameter$vflat), 
               format(round(object$angularmomentum$j_observed, 2), nsmall = 2), 
               format(round(object$angularmomentum$j_model, 2), nsmall = 2), 
               format(round(object$angularmomentum$j_extrapolated, 2), nsmall = 2),
               format(round(object$angularmomentum$j_extrapolated_err, 2), nsmall = 2),
               format(round((object$angularmomentum$j_extrapolated-object$angularmomentum$j_observed)/object$angularmomentum$j_extrapolated, 2), nsmall = 2) )

result_file = '/Users/gsavorgnan/angularmomentum/results.dat'


if (!file.exists(result_file)) {
  cat('# ID   Gal   Obs   OSI   GMOS   OSI+GMOS   incl   incl_forced   h_disk   r_flat   v_flat   j_obs   j_mod   j_tot   err_j_tot   f_extr   \n', file=result_file)
}
cat(result_row, file=result_file, sep='   ', append=T)  
cat('\n', file=result_file, append=T)
}

#galaxy_list = c('G041', 'C222', 'D135', 'G202', 'SDSS013527', 'SDSS234657')
galaxy_list = c('D135')


id = 1

for (galaxy in galaxy_list) {
  
  main(galaxy, id, T, F)
  main(galaxy, id, F, T)
  main(galaxy, id, T, T)
  
  id = id+1
}  
