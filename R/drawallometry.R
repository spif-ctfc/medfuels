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
drawallometry<-function(species, type="total") {
  if(type=="total") {
    sp_params = get("sp_params_total")
    group_params = get("group_params_total")
  } else if(type=="fine") {
    sp_params = get("sp_params_fine")
    group_params = get("group_params_fine")
  }
  vol = seq(sp_params[species,"minVol"], sp_params[species,"maxVol"],length.out = 100)
  weight = sp_params[species,"a"]*vol^sp_params[species,"b"]
  df = data.frame(volume = vol, weight = weight)
  ggplot(df, aes(volume, weight)) +
    geom_line()
}