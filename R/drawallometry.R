#' Draw allometric equation
#'
#' @param species the specie or functional group to be drawn 
#' @param type 'total'  (total fuel), 'fine' (fine fuels), or 'both'
#'
#' @return
#' @export
#'
#' @examples
#' drawallometry("Cistus albidus")
drawallometry<-function(species, type="total", log = "") {
  type = match.arg(type, c("total","fine","both"))
  if(type=="total" || type=="both") {
    sp_params = get("sp_params_total")
    group_params = get("group_params_total")
    vol_total = seq(sp_params[species,"minVol"], sp_params[species,"maxVol"],length.out = 100)
    weight_total = sp_params[species,"a"]*vol_total^sp_params[species,"b"]
    vars_total = (weight_total^2)*sp_params[species,"gamma_disp"]
  } 
  if(type=="fine" || type=="both") {
    sp_params = get("sp_params_fine")
    group_params = get("group_params_fine")
    vol_fine = seq(sp_params[species,"minVol"], sp_params[species,"maxVol"],length.out = 100)
    weight_fine = sp_params[species,"a"]*vol_fine^sp_params[species,"b"]
    vars_fine = (weight_fine^2)*sp_params[species,"gamma_disp"]
  } 
  if(type=="total") {
    df = data.frame(volume = vol_total, weight = weight_total, vars = vars_total,
                    fuel="total")
  } else if(type=="fine") {
    df = data.frame(volume = vol_fine, weight = weight_fine, vars = vars_fine, 
                    fuel="fine")
  } else {
    df = data.frame(volume = c(vol_total,vol_fine), 
                    weight = c(weight_total,weight_fine),
                    vars = c(vars_total,vars_fine),
                    fuel = c(rep("total", length(vol_total)), 
                             rep("fine",length(vol_fine))))
  }
  rate = df$weight/df$vars
  shape = df$weight * rate
  df$q05 = qgamma(0.05, shape, rate)
  df$q95 = qgamma(0.95, shape, rate)
  h <- ggplot(df, aes(volume, weight))
  if(type=="both") {
    h<-h+ geom_ribbon(aes(ymin = q05, ymax = q95, fill = fuel), alpha=0.5)+
      geom_line(aes(color=fuel)) + labs(color="Fraction", fill="Fraction")
  } else {
    h<-h+ geom_ribbon(aes(ymin = q05, ymax = q95), fill="grey70", alpha=0.5) +
      geom_line()
  }
  h<- h + xlab(label = expression(paste("Volume  ", (m^3))))+
          ylab(label = expression(paste("Dry weight  ", (kg))))
  if(log=="x" || log=="xy") h <- h+ scale_x_continuous(trans="log10")
  if(log=="y" || log=="xy") h <- h+ scale_y_continuous(trans="log10")
  h<-h+labs(title=species)
  h
}