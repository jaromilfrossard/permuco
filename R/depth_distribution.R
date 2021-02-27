depth_distribution <- function(distribution, head_mat = NULL, tail_mat = NULL){
  if(is.null(head_mat)&is.null(tail_mat)){
    stop("You need to specify a head_mat or tail_mat argument")
  }else if((!is.null(head_mat))&(is.null(tail_mat))){
    out <- depth_distribution_head(distribution, head_mat)
  }else if((is.null(head_mat))&(!is.null(tail_mat))){
    out <- depth_distribution_tail(distribution, tail_mat)
  }else if((!is.null(head_mat))&(!is.null(tail_mat))){
    out <- depth_distribution_unique(distribution, head = head_mat, tail = tail_mat)

  }
  return(out)

}