#'@importFrom Matrix KhatriRao
khatrirao <- function(a,b,c = NULL){
  if(!is.null(c)){
    if(dim(c)[2]==0){c=NULL}}
  if(dim(b)[2]==0){return(a)}

  out = t(as.matrix(KhatriRao(t(a),t(b))))
  if(!is.null(c)){
    out= t(as.matrix(KhatriRao(t(out),t(c))))
  }
  out}