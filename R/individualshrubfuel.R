#' Shrub fuel biomass
#'
#' Calculates dry weight (biomass, in kg) of total or fine fuels corresponding to individual shrub data
#'
#' @param x data frame with columns 'plot', 'species', 'H' (height in cm), 'D1' and 'D2' (in cm)
#' @param type either 'total'  (total fuel) of 'fine' (fine fuels)
#' @param allometric wether to use allometric equations or bulk density estimates
#' @param excludeSSP excludes subspecies information for species matching
#' @param var a flag to indicate that variance of estimates is desired
#' @param agg aggregation of results. Either 'none', 'species', 'speciesplot', 'plotspecies' or 'plot'
#' @param customParams custom allometry parameter table (for species not in default params)
#' @param na.rm whether to exclude missing values when aggregating biomass
#'
#' @return a vector of dry weight (biomass) in kg
#' @export
#'
#' @examples
#' plot = c(1,1,2,2,2)
#' species = c("Erica arborea","Cistus albidus", "Erica arborea", "Chamaerops humilis", "Unknown")
#' H = c(60,200,100,250,100)
#' D1 = c(10,100,30, 50,25)
#' D2 = D1
#' x = data.frame(plot, species, H, D1, D2)
#'
#' individualshrubfuel(x)

individualshrubfuel <- function(x, type= "total",  allometric = TRUE, excludeSSP = TRUE, var = FALSE,
                             agg = "none",  customParams = NULL, na.rm = TRUE) {
  type = match.arg(type, c("total","fine"))
  agg = match.arg(agg, c("none", "species", "plot", "plotspecies", "speciesplot"))
  x = as.data.frame(x)
  vars = names(x)
  if(!("plot" %in% vars)) stop("Variable 'plot' needed in 'x'")
  if(!("species" %in% vars)) stop("Variable 'species' needed in 'x'")
  if(!("H" %in% vars)) stop("Variable 'H' needed in 'x'")
  if(!("D1" %in% vars)) stop("Variable 'D1' needed in 'x'")
  if(!("D2" %in% vars)) stop("Variable 'D2' needed in 'x'")
  h = x$H/100 #from cm to m
  d1 = x$D1/100
  d2 = x$D2/100
  vol = h*pi*(d1/2)*(d2/2)

  sp = as.character(x$species)
  if(excludeSSP) {
    s = strsplit(as.character(sp), " ")
    sp = unlist(lapply(s, function(x) {paste(x[1:min(2,length(x))], collapse=" ")}))
  }

  if(type=="total") {
    sp_params = get("sp_params_total")
    group_params = get("group_params_total")
  } else if(type=="fine") {
    sp_params = get("sp_params_fine")
    group_params = get("group_params_fine")
  }

  sp_list = row.names(sp_params)
  nind = nrow(x)
  weight = rep(NA,nind)
  vars =  rep(NA,nind)
  for(i in 1:nind) {
    if(sp[i] %in% sp_list) {
      if(allometric) {
        if(vol[i] > sp_params[sp[i],"maxVol"]) warning(paste0("Volume '", vol[i],"' outside the calibration range for '", sp[i],"'"))
        weight[i] = sp_params[sp[i],"a"]*vol[i]^sp_params[sp[i],"b"]
        vars[i] = (weight[i]^2)*sp_params[sp[i],"gamma_disp"]
      }
      else {
        weight[i] = sp_params[sp[i],"BD"]*vol[i]
        vars[i] = NA
      }
    } else {
      gr = .getSpeciesGroup(sp[i])
      if(!is.na(gr)) {
        if(allometric) {
          if(vol[i] > group_params[gr,"maxVol"]) warning(paste0("Volume '", vol[i],"' outside the calibration range for '", gr[i],"'"))
          weight[i] = group_params[gr,"a"]*vol[i]^group_params[gr,"b"]
          vars[i] = (weight[i]^2)*group_params[gr,"gamma_disp"]
        }
        else {
          weight[i] = group_params[gr,"BD"]*vol[i]
          vars[i] = NA
        }
      } else {
        warning(paste0("Species '", sp[i],"' not found in parameter file for biomass!"))
      }
    }
  }

  if(agg=="species") {
    weight = tapply(weight, x$species, FUN = sum, na.rm=na.rm)
  } else if(agg=="plot") {
    weight = tapply(weight, x$plot, FUN = sum, na.rm=na.rm)
  } else if(agg=="plotspecies") {
    weight = aggregate(weight, by=list(x$species, x$plot), FUN = sum, na.rm=na.rm, drop=TRUE, simplify=TRUE)
    weight = weight[,c(2,1,3)]
    names(weight)<-c("plot","species", "biomass")
  } else if(agg=="speciesplot") {
    weight = aggregate(weight, by=list(x$plot, x$species), FUN = sum, na.rm=na.rm, drop=TRUE, simplify=TRUE)
    weight = weight[,c(2,1,3)]
    names(weight)<-c("species","plot", "biomass")
  } else {
    names(weight) = row.names(x)
  }
  if(!var) {
    res = weight
  } else {
    if(agg=="species") {
      vars = tapply(vars, x$species, FUN = sum, na.rm=na.rm)
      res = data.frame(biomass= weight, var = vars)
    } else if(agg=="plot") {
      vars = tapply(vars, x$plot, FUN = sum, na.rm=na.rm)
      res = data.frame(biomass= weight, var = vars)
    } else if(agg=="plotspecies") {
      vars = aggregate(vars, by=list(x$species, x$plot), FUN = sum, na.rm=na.rm, drop=TRUE, simplify=TRUE)
      vars = vars[,c(3)]
      res = weight
      res$var = vars
    } else if(agg=="speciesplot") {
      vars = aggregate(vars, by=list(x$species, x$plot), FUN = sum, na.rm=na.rm, drop=TRUE, simplify=TRUE)
      vars = vars[,c(3)]
      res = weight
      res$var = vars
    } else {
      res = data.frame(biomass= weight, var = vars)
    }
  }
  return(res)
}

.getSpeciesGroup<-function(sp) {
  species_groups = get("species_groups")
  if(sp %in% species_groups$Name) {
    return(species_groups[which(species_groups$Name==sp),"shrub type"])
  }
  return(NA)
}
