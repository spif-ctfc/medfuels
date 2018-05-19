#' Draw allometric equation
#'
#' @param species the specie or functional group to be drawn 
#' @param type either 'total'  (total fuel) of 'fine' (fine fuels)
#'
#' @return
#' @export
#'
#' @examples
#' 
drawallometry<-function(species, type="total", log = "") {
  if(type=="total") {
    sp_params = get("sp_params_total")
    group_params = get("group_params_total")
  } else if(type=="fine") {
    sp_params = get("sp_params_fine")
    group_params = get("group_params_fine")
  }
  vol = seq(sp_params[species,"minVol"], sp_params[species,"maxVol"],length.out = 100)
  weight = sp_params[species,"a"]*vol^sp_params[species,"b"]
  vars = (weight^2)*sp_params[species,"gamma_disp"]
  rate = weight/vars
  shape = weight * rate
  q05 = qgamma(0.05, shape, rate)
  q95 = qgamma(0.95, shape, rate)
  df = data.frame(volume = vol, weight = weight, q05 = q05, q95 = q95)
  h <- ggplot(df, aes(volume, weight)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), fill="grey70", alpha=0.5) +
    geom_line() +
    xlab(label = expression(paste("Volume  ", (m^3))))+
    ylab(label = expression(paste("Dry weight  ", (kg))))
  if(log=="x" || log=="xy") h <- h+ scale_x_continuous(trans="log10")
  if(log=="y" || log=="xy") h <- h+ scale_y_continuous(trans="log10")
  h
}