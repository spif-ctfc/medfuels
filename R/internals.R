.getSpeciesGroup<-function(sp) {
  species_groups = get("species_groups")
  if(sp %in% species_groups$name) {
    return(species_groups[which(species_groups$name==sp),"shrub type"])
  }
  return(NULL)
}
.translateSynonym<-function(sp) {
  s = strsplit(species_synonyms,"/")
  ind = which(!is.na(unlist(lapply(s, function(x) {x[x==sp]}))))
  if(length(ind)==1) return(names(ind))
  return(NULL)
}
