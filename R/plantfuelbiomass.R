#' Plant fuel biomass
#'
#' Calculates dry weight (biomass, in kg) of total or fine fuels corresponding to plant data
#'
#' @param x data frame with columns 'plot', 'individual', 'species', 'H' (height in cm), 'D1' and 'D2' (in cm)
#' @param type either 'total'  (total fuel) or 'fine' (fine fuels)
#' @param agg aggregation of results. Either 'none', 'species' or 'plot'
#' @param customParams custom allometry parameter table (for species not in default params)
#' @param na.rm whether to exclude missing values when aggregating biomass
#'
#' @return a vector of dry weight (biomass) in kg
#' @export
#'
#' @examples
#' plot = c(1,1,2,2,2)
#' individual = c(1,2,1,2,3)
#' species = c("Erica arborea","Cistus albidus", "Erica arborea", "Chamaerops humilis", "Unknown")
#' H = c(60,200,100,250,100)
#' D1 = c(10,100,30, 50,25)
#' D2 = D1
#' x = data.frame(plot, individual, species, H, D1, D2)
#'
#' plantfuelbiomass(x)

plantfuelbiomass <- function(x, type= "total", agg = "none", customParams = NULL, na.rm = TRUE) {
  type = match.arg(type, c("total","fine"))
  agg = match.arg(agg, c("none", "species", "plot"))
  x = as.data.frame(x)
  vars = names(x)
  if(!("plot" %in% vars)) stop("Variable 'plot' needed in 'x'")
  if(!("individual" %in% vars)) stop("Variable 'individual' needed in 'x'")
  if(!("species" %in% vars)) stop("Variable 'species' needed in 'x'")
  if(!("H" %in% vars)) stop("Variable 'H' needed in 'x'")
  if(!("D1" %in% vars)) stop("Variable 'D1' needed in 'x'")
  if(!("D2" %in% vars)) stop("Variable 'D2' needed in 'x'")
  h = x$H/100 #from cm to m
  d1 = x$D1/100
  d2 = x$D2/100
  vol = h*pi*(d1/2)*(d2/2)

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
  nind = nrow(x)
  load = rep(NA,nind)
  for(i in 1:nind) {
    sp = as.character(x$species[i])
    if(sp %in% sp_list) {
      load[i] = sp_params[sp,"a"]*vol[i]^sp_params[sp,"b"]
    } else {
      gr = .getSpeciesGroup(sp)
      if(!is.na(gr)) {
        load[i] = group_params[gr,"a"]*vol[i]^group_params[gr,"b"]
      } else {
        warning(paste0("Species '", sp,"' not found in parameter file!"))
      }
    }
  }

  if(agg=="species") {
    load = tapply(load, x$species, FUN = sum, na.rm=na.rm)
  } else if(agg=="plot") {
    load = tapply(load, x$plot, FUN = sum, na.rm=na.rm)
  } else {
    names(load) = row.names(x)
  }
  return(load)
}

.getSpeciesGroup<-function(sp) {
  data("species_groups")
  if(sp %in% species_groups$Name) {
    return(species_groups[which(species_groups$Name==sp),"shrub type"])
  }
  return(NA)
}
