#' Shrub fuel biomass
#'
#' Calculates dry weight (biomass, in kg) of total or fine fuels corresponding to individual shrub data
#'
#' @param x data frame with columns 'plot', 'species', 'H' (height in cm), 'D1' and 'D2' (in cm). If 'D2' is ommitted then
#'          shrub crowns are assumed to be circular (i.e. D2 = D1). Column 'group' may be used to specify a functional group to be used when
#'          the species has not a specific allometry (see details).
#' @param type either 'total'  (total fuel) of 'fine' (fine fuels)
#' @param allometric wether to use allometric equations or bulk density estimates
#' @param excludeSSP excludes subspecies information for species matching
#' @param var a flag to indicate that variance of estimates is desired
#' @param agg aggregation of results. Either 'none', 'species', 'speciesplot', 'plotspecies' or 'plot'
#' @param customParams custom allometry parameter table (for species not in default params)
#' @param na.rm whether to exclude missing values when aggregating biomass
#'
#' @details The function determines the allometry to be applied using the following rules, sequentially:
#' \enumerate{
#'   \item{If the species name is included in the list of species with parameter values, it takes the parameters from that species.}
#'   \item{If the species name is a synonym for a species included the list of species with parameter values, it takes the parameters from that species (and gives a warning).}
#'   \item{If the user has specified a column 'group' and the value is included in the list of valid groups, it takes the parameters from that group.}
#'   \item{If the species is listed within the checklist in 'species_groups', it takes the parameters from the group specified in that file.}
#'   \item{Otherwise, it gives a warning and no allometry is applied.}
#' }
#'
#' @return a vector of dry weight (biomass) in kg
#' @export
#'
#' @examples
#' plot = c(1,1,2,2,2)
#' species = c("Erica arborea","Genista cinerea", "Erica arborea", "Chamaerops humilis", "Genista balansae")
#' H = c(60,70,200,100,10)
#' D1 = c(60,40,100, 100,25)
#' D2 = D1
#' x = data.frame(plot, species, H, D1, D2, stringsAsFactors=FALSE)
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

  h = x$H/100 #from cm to m
  d1 = x$D1/100
  if(!("D2" %in% vars)) {
    warning("Assuming D2 = D1 (circular crown projection)")
    d2 = d1
  } else {
    d2 = x$D2/100
  }
  vol = h*pi*(d1/2)*(d2/2)

  hasGroup = ("group" %in% vars)
  if(hasGroup) gr = as.character(x$group)
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
  gr_list = row.names(group_params)
  nind = nrow(x)
  weight = rep(NA,nind)
  vars =  rep(NA,nind)
  for(i in 1:nind) {
    spi = NULL
    if(sp[i] %in% sp_list) {
      spi = sp[i]
    } else { #try translating the synonym
      spi = .translateSynonym(sp[i])
      if(!is.null(spi)) warning(paste0("Input species '", sp[i],"' translated to '",spi,"'."))
    }
    if(!is.null(spi)) {
      if(allometric) {
        if(vol[i] > sp_params[spi,"maxVol"]) warning(paste0("Volume '", vol[i],"' outside the calibration range for '", spi,"'"))
        weight[i] = sp_params[spi,"a"]*vol[i]^sp_params[spi,"b"]
        vars[i] = (weight[i]^2)*sp_params[spi,"gamma_disp"]
      }
      else {
        weight[i] = sp_params[spi,"BD"]*vol[i]
        vars[i] = (sp_params[spi,"BD.sd"]^2)*vol[i]
      }
    } else {
      gri = NULL
      if(hasGroup) {
        if(gr[i] %in% gr_list) gri = gr[i]
      }
      if(is.null(gri)) {
        gri = .getSpeciesGroup(sp[i])
      }
      if(!is.null(gri)) {
        if(allometric) {
          if(vol[i] > group_params[gri,"maxVol"]) warning(paste0("Volume '", vol[i],"' outside the calibration range for '", gri,"'"))
          weight[i] = group_params[gri,"a"]*vol[i]^group_params[gri,"b"]
          vars[i] = (weight[i]^2)*group_params[gri,"gamma_disp"]
        }
        else {
          weight[i] = group_params[gri,"BD"]*vol[i]
          vars[i] = (group_params[gri,"BD.sd"]^2)*vol[i]
        }
      } else {
        warning(paste0("Input species '", sp[i],"' not found in parameter file for biomass!"))
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

