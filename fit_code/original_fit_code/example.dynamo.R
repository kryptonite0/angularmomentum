# load, process, fit, and plot galaxies
index = 46#46
if (index>10) {
  slave_index = index%%10
  master_index = (index-slave_index)/10
  master_object = dynamo.load(master_index) # load fits
  slave_object = dynamo.load(slave_index) # load fits
  master_object = disk.fit(master_object, show.fit = F, show.details = F)
  slave_object = disk.fit(slave_object, show.fit = F, show.details = F)
  object = merge.objects(master_object,slave_object)
} else {
  object = dynamo.load(index)
}
object = disk.fit(object, show.fit = T, show.details = T, iterations = 0) # fit disk
  
# print text
cat('===============================================================\n')
cat(sprintf('inclination       = %6.1f deg\n',object$fit$deprojection$parameter$inclination*180/pi))
cat(sprintf('PA from north     = %6.1f deg\n',object$fit$deprojection$parameter$positionangle_from_NCP*180/pi))
cat(sprintf('rdisk             = %6.1f kpc\n',object$fit$density_model$parameter$single_exponential_r*object$data$kpc_per_pixel))
cat(sprintf('rflat             = %6.1f kpc\n',object$fit$velocity_model$parameter$rflat*object$data$kpc_per_pixel))
cat(sprintf('vflat             = %6.1f km/s\n',round(object$fit$velocity_model$parameter$vflat)))
cat(sprintf('j_observed        = %6.1f kpc km/s\n',object$angularmomentum$j_observed))
cat(sprintf('j_total           = %6.1f kpc km/s\n',object$angularmomentum$j_extrapolated))
cat(sprintf('j_total_err       = %6.1f kpc km/s\n',object$angularmomentum$j_extrapolated_err))
cat(sprintf('extrapol fraction = %6.1f%s\n',(object$angularmomentum$j_extrapolated-object$angularmomentum$j_observed)/object$angularmomentum$j_extrapolated*100,'%'))
