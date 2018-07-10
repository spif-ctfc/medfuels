#' Shrub species fuel loading
#'
#' Calculates loading (kg/m2) of total or fine fuels corresponding to shrub species data
#'
#' @param x data frame with columns 'plot', 'species', 'H' (mean height in cm), 'C' (cover in percent). Column 'group' may be used to specify a functional group to be used when
#'          the species has not a specific allometry (see details).
#' @param type either 'total'  (total fuel) or 'fine' (fine fuels)
#' @param allometric wether to use allometric equations or bulk density estimates
#' @param excludeSSP excludes subspecies information for species matching
#' @param var a flag to indicate that variance of estimates is desired
#' @param agg aggregation of results. Either 'none' or 'plot'
#' @param customParams custom allometry parameter table (for species not in default params)
#' @param na.rm whether to exclude missing values when aggregating loading
#'
#' @return a vector with species loading values (kg/m2)
#' @export
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
#' @examples
#' plot = c(1,1,2,2,2)
#' species = c("Erica arborea","Genista cinerea", "Erica arborea", "Chamaerops humilis", "Sarothamnus scoparius")
#' H = c(60,70,200,100,10)
#' C = c(10,10,30, 20,5)
#' x = data.frame(plot, species, H, C)
#'
#' shrubspeciesfuel(x)
shrubspeciesfuel <- function(x, type= "total", allometric = TRUE, excludeSSP = TRUE, var = FALSE,
                               agg = "none", customParams = NULL, na.rm = TRUE) {
  type = match.arg(type, c("total","fine"))
  agg = match.arg(agg, c("none", "species", "plot"))
  x = as.data.frame(x)
  vars = names(x)
  if(!("plot" %in% vars)) stop("Variable 'plot' needed in 'x'")
  if(!("species" %in% vars)) stop("Variable 'species' needed in 'x'")
  if(!("H" %in% vars)) stop("Variable 'H' needed in 'x'")
  if(!("C" %in% vars)) stop("Variable 'C' needed in 'x'")

  hasGroup = ("group" %in% vars)
  if(hasGroup) gr = as.character(x$group)
  sp = as.character(x$species)
  if(excludeSSP) {
    s = strsplit(as.character(sp), " ")
    sp = unlist(lapply(s, function(x) {paste(x[1:min(2,length(x))], collapse=" ")}))
  }

  hm = x$H #hm in cm

  sp_params = get("sp_params_area")
  group_params = get("group_params_area")
  sp_list = row.names(sp_params)
  gr_list = row.names(group_params)

  nrec = nrow(x)
  if(allometric) {
    area = rep(NA,nrec)
    for(i in 1:nrec) {
      spi = NULL
      if(sp[i] %in% sp_list) {
        spi = sp[i]
      } else { #try translating the synonym
        spi = .translateSynonym(sp[i])
        if(!is.null(spi)) warning(paste0("Input species '", sp[i],"' translated to '",spi,"'."))
      }
      if(!is.null(spi)) {
        if(hm[i] > sp_params[spi,"maxH"]) warning(paste0("Height '", hm[i],"' outside the calibration range for '", spi,"'"))
        area[i] = sp_params[spi,"a"]*hm[i]^2 # area in cm2
      } else {
        gri = NULL
        if(hasGroup) {
          if(gr[i] %in% gr_list) gri = gr[i]
        }
        if(is.null(gri)) {
          gri = .getSpeciesGroup(sp[i])
        }
        if(!is.null(gri)) {
          if(hm[i] > group_params[gri,"maxH"]) warning(paste0("Height '", hm[i],"' outside the calibration range for '", gri,"'"))
          area[i] = group_params[gri,"a"]*hm[i]^2
        } else {
          warning(paste0("Input species '", sp[i],"' not found in parameter file for area!"))
        }
      }
    }

    N = x$C*100/area #Density in ind/m2
    vol = (area*hm)/(10^6) #vol in m3
  } else {
    vol = (x$C/100)*(hm/100) #vol in m3/m2
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
  weight = rep(NA,nrec)
  vars =  rep(NA,nrec)
  for(i in 1:nrec) {
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
      } else {
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

  if(allometric) {
    load = weight*N
    varload = vars*N
  }
  else {
    load = weight
    varload = vars
  }

  if(agg=="plot") {
    load = tapply(load, x$plot, FUN = sum, na.rm=na.rm)
    varload = tapply(varload, x$plot, FUN = sum, na.rm=na.rm)
  } else {
    names(load) = row.names(x)
  }

  if(!var) {
    res = load
  } else {
    res = data.frame(loading= load, var = varload)
  }
  return(res)
}
