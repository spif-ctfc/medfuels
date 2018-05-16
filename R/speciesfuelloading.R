#' Species fuel loading
#'
#' Calculates loading (kg/m2) of total or fine fuels corresponding to species data
#'
#' @param x data frame with columns 'plot', 'species', 'H' (mean height in cm), 'C' (cover in percent)
#' @param type either 'total'  (total fuel) or 'fine' (fine fuels)
#' @param agg aggregation of results. Either 'none' or 'plot'
#' @param customParams custom allometry parameter table (for species not in default params)
#' @param na.rm whether to exclude missing values when aggregating loading
#'
#' @return a vector with loading values (kg/m2)
#' @export
#'
#' @examples
#' plot = c(1,1,2,2,2)
#' species = c("Erica arborea","Cistus albidus", "Erica arborea", "Chamaerops humilis", "Unknown")
#' H = c(60,200,100,250,100)
#' C = c(10,100,30, 50,25)
#' x = data.frame(plot, species, H, C)
speciesfuelloading <- function(x, type= "total", agg = "none", customParams = NULL, na.rm = TRUE) {
  type = match.arg(type, c("total","fine"))
  agg = match.arg(agg, c("none", "species", "plot"))
  x = as.data.frame(x)
  vars = names(x)
  if(!("plot" %in% vars)) stop("Variable 'plot' needed in 'x'")
  if(!("species" %in% vars)) stop("Variable 'species' needed in 'x'")
  if(!("H" %in% vars)) stop("Variable 'H' needed in 'x'")
  if(!("C" %in% vars)) stop("Variable 'C' needed in 'x'")

  hm = x$H #hm in cm

  data("sp_params_area")
  data("group_params_area")
  sp_params = sp_params_area
  group_params = group_params_area

  sp_list = row.names(sp_params)
  nrec = nrow(x)
  area = rep(NA,nrec)
  for(i in 1:nrec) {
    sp = as.character(x$species[i])
    if(sp %in% sp_list) {
      area[i] = sp_params[sp,"a"]*hm[i]^2 # area in cm2
    } else {
      gr = .getSpeciesGroup(sp)
      if(!is.na(gr)) {
        area[i] = group_params[gr,"a"]*hm[i]^2
      } else {
        warning(paste0("Species '", sp,"' not found in parameter file for area!"))
      }
    }
  }

  N = x$C*100/area #Density in ind/m2
  vol = (area*hm)/(10^6) #vol in m3

  if(type=="total") {
    data("sp_params_total")
    data("group_params_total")
    sp_params = sp_params_total
    group_params = group_params_total
  } else if(type=="fine") {
    data("sp_params_fine")
    data("group_params_fine")
    sp_params = sp_params_fine
    group_params = group_params_fine
  }

  sp_list = row.names(sp_params)
  weight = rep(NA,nrec)
  for(i in 1:nrec) {
    sp = as.character(x$species[i])
    if(sp %in% sp_list) {
      weight[i] = sp_params[sp,"a"]*vol[i]^sp_params[sp,"b"]
    } else {
      gr = .getSpeciesGroup(sp)
      if(!is.na(gr)) {
        weight[i] = group_params[gr,"a"]*vol[i]^group_params[gr,"b"]
      } else {
        warning(paste0("Species '", sp,"' not found in parameter file for biomass!"))
      }
    }
  }

  load = weight*N

  if(agg=="plot") {
    load = tapply(load, x$plot, FUN = sum, na.rm=na.rm)
  } else {
    names(load) = row.names(x)
  }
  return(load)
}
