#' Transform shrub individual data to species data
#'
#' Transforms shrub individual data (plant heights and diameters) into species level data (average height and cover).
#'
#' @param x data frame with columns 'plot', 'species', 'H' (height in cm), 'D1' and 'D2' (in cm)
#' @param sampledarea sampled area in squared meters
#' @param maxcover maximum allowed cover (set to NA to avoid truncation)
#' @param na.rm whether to exclude missing values when averaging heights and adding areas
#'
#' @return data frame with columns 'plot', 'species', 'H' (mean height in cm), 'C' (cover in percent)
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
#' individualshrub2species(x, sampledarea = 2)
individualshrub2species<-function(x, sampledarea = 20, maxcover = 100, na.rm = TRUE) {
  x = as.data.frame(x)
  vars = names(x)
  if(!("plot" %in% vars)) stop("Variable 'plot' needed in 'x'")
  if(!("species" %in% vars)) stop("Variable 'species' needed in 'x'")
  if(!("H" %in% vars)) stop("Variable 'H' needed in 'x'")
  if(!("D1" %in% vars)) stop("Variable 'D1' needed in 'x'")
  if(!("D2" %in% vars)) stop("Variable 'D2' needed in 'x'")


  plots = levels(as.factor(x$plot))
  shspdata = NULL
  for(i in 1:length(plots)) {
    pc = x[as.character(x$plot)==plots[i],]
    area = pi*((pc$D1/200)*(pc$D2/200)) #area in m2
    species = as.factor(as.character(pc$species))
    shpc = data.frame(plot = rep(plots[i], length(levels(species))),
                      species = levels(species),
                      H = tapply(pc$H, species, FUN=mean, na.rm=na.rm),
                      C = 100*tapply(area, species, FUN=sum, na.rm=na.rm)/sampledarea)
    row.names(shpc) = NULL
    if(!is.null(shspdata)) shspdata = rbind(shspdata, shpc)
    else shspdata = shpc
  }
  if(!is.na(maxcover)) shspdata$C = pmin(maxcover, shspdata$C)
  return(shspdata)
}
