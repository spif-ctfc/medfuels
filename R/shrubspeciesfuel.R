#' Shrub species fuel loading
#'
#' Calculates loading (kg/m2) of total or fine fuels corresponding to shrub species data
#'
#' @param x data frame with columns 'plot', 'species', 'H' (mean height in cm), 'C' (cover in percent)
#' @param type either 'total'  (total fuel) or 'fine' (fine fuels)
#' @param allometric wether to use allometric equations or bulk density estimates
#' @param excludeSSP excludes subspecies information for species matching
#' @param var a flag to indicate that variance of estimates is desired
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

  sp = as.character(x$species)
  if(excludeSSP) {
    s = strsplit(as.character(sp), " ")
    sp = unlist(lapply(s, function(x) {paste(x[1:min(2,length(x))], collapse=" ")}))
  }

  hm = x$H #hm in cm

  sp_params = get("sp_params_area")
  group_params = get("group_params_area")

  sp_list = row.names(sp_params)
  nrec = nrow(x)
  if(allometric) {
    area = rep(NA,nrec)
    for(i in 1:nrec) {
      if(sp[i] %in% sp_list) {
        if(hm[i] > sp_params[sp[i],"maxH"]) warning(paste0("Height '", hm[i],"' outside the calibration range for '", sp[i],"'"))
        area[i] = sp_params[sp[i],"a"]*hm[i]^2 # area in cm2
      } else {
        gr = .getSpeciesGroup(sp[i])
        if(!is.na(gr)) {
          if(hm[i] > group_params[gr,"maxH"]) warning(paste0("Height '", hm[i],"' outside the calibration range for '", gr,"'"))
          area[i] = group_params[gr,"a"]*hm[i]^2
        } else {
          warning(paste0("Species '", sp[i],"' not found in parameter file for area!"))
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
  weight = rep(NA,nrec)
  vars =  rep(NA,nrec)
  for(i in 1:nrec) {
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
