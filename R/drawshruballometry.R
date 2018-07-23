#' Draw shrub allometric equation
#'
#' @param species The shrub species or functional group to be drawn
#' @param type Allometry to draw: 'total'  (total fuel vs. volume), 'fine' (fine fuels vs. volume), 'both' (total and fine fuels vs. volume) or 'area' (area vs. height)
#' @param log Specifies logarithmic axis: 'x', 'y' or 'xy'
#'
#' @return
#' @export
#'
#' @examples
#' drawshruballometry("Cistus albidus")
drawshruballometry<-function(species, type="total", log = "") {
  type = match.arg(type, c("total","fine","both", "area"))
  if(type=="total" || type=="both") {
    sp_params = get("sp_params_total")
    group_params = get("group_params_total")
    if(species %in% row.names(sp_params)) {
      vol_total = seq(sp_params[species,"minVol"], sp_params[species,"maxVol"],length.out = 100)
      weight_total = sp_params[species,"a"]*vol_total^sp_params[species,"b"]
      vars_total = (weight_total^2)*sp_params[species,"gamma_disp"]
    } else if(species %in% row.names(group_params)) {
      vol_total = seq(group_params[species,"minVol"], group_params[species,"maxVol"],length.out = 100)
      weight_total = group_params[species,"a"]*vol_total^group_params[species,"b"]
      vars_total = (weight_total^2)*group_params[species,"gamma_disp"]
    } else {
      stop("Not a valid species or group name")
    }
  }
  if(type=="fine" || type=="both") {
    sp_params = get("sp_params_fine")
    group_params = get("group_params_fine")
    if(species %in% row.names(sp_params)) {
      vol_fine = seq(sp_params[species,"minVol"], sp_params[species,"maxVol"],length.out = 100)
      weight_fine = sp_params[species,"a"]*vol_fine^sp_params[species,"b"]
      vars_fine = (weight_fine^2)*sp_params[species,"gamma_disp"]
    } else if(species %in% row.names(group_params)) {
      vol_fine = seq(group_params[species,"minVol"], group_params[species,"maxVol"],length.out = 100)
      weight_fine = group_params[species,"a"]*vol_fine^group_params[species,"b"]
      vars_fine = (weight_fine^2)*group_params[species,"gamma_disp"]
    } else {
      stop("Not a valid species or group name")
    }
  }
  if(type=="total") {
    df = data.frame(volume = vol_total, weight = weight_total, vars = vars_total,
                    fuel="total")
  } else if(type=="fine") {
    df = data.frame(volume = vol_fine, weight = weight_fine, vars = vars_fine,
                    fuel="fine")
  } else if(type=="both"){
    df = data.frame(volume = c(vol_total,vol_fine),
                    weight = c(weight_total,weight_fine),
                    vars = c(vars_total,vars_fine),
                    fuel = c(rep("total", length(vol_total)),
                             rep("fine",length(vol_fine))))
  }
  if(type=="area") {
    sp_params = get("sp_params_area")
    group_params = get("group_params_area")
    if(species %in% row.names(sp_params)) {
      H = seq(sp_params[species,"minH"], sp_params[species,"maxH"],length.out = 100)
      area = sp_params[species,"a"]*H^sp_params[species,"b"]
      vars = (area^2)*sp_params[species,"gamma_disp"]
    } else if(species %in% row.names(group_params)) {
      H = seq(group_params[species,"minH"], group_params[species,"maxH"],length.out = 100)
      area = group_params[species,"a"]*H^group_params[species,"b"]
      vars = (area^2)*group_params[species,"gamma_disp"]
    } else {
      stop("Not a valid species or group name")
    }
    df = data.frame(height = H, area = area, vars = vars)
    rate = df$area/df$vars
    shape = df$area * rate
    df$q05 = qgamma(0.05, shape, rate)
    df$q95 = qgamma(0.95, shape, rate)
    h <- ggplot(df, aes(df$height, df$area))
    h<-h+ geom_ribbon(aes(ymin = df$q05, ymax = df$q95), fill="grey70", alpha=0.5) +
        geom_line()
    h<- h + xlab(label = expression(paste("Height  ", (cm))))+
      ylab(label = expression(paste("Area  ", (cm^2))))
    if(log=="x" || log=="xy") h <- h+ scale_x_continuous(trans="log10")
    if(log=="y" || log=="xy") h <- h+ scale_y_continuous(trans="log10")
    h<-h+labs(title=species)
  } else {
    rate = df$weight/df$vars
    shape = df$weight * rate
    df$q05 = qgamma(0.05, shape, rate)
    df$q95 = qgamma(0.95, shape, rate)
    h <- ggplot(df, aes(df$volume, df$weight))
    if(type=="both") {
      h<-h+ geom_ribbon(aes(ymin = df$q05, ymax = df$q95, fill = df$fuel), alpha=0.5)+
        geom_line(aes(color=df$fuel)) + labs(color="Fraction", fill="Fraction")
    } else {
      h<-h+ geom_ribbon(aes(ymin = df$q05, ymax = df$q95), fill="grey70", alpha=0.5) +
        geom_line()
    }
    h<- h + xlab(label = expression(paste("Volume  ", (m^3))))+
      ylab(label = expression(paste("Dry weight  ", (kg))))
    if(log=="x" || log=="xy") h <- h+ scale_x_continuous(trans="log10")
    if(log=="y" || log=="xy") h <- h+ scale_y_continuous(trans="log10")
    h<-h+labs(title=species)
  }
  h
}
