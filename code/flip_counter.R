count <- 0
in_progress <- FALSE
for(i in 1:(length(df$POS_START)-1)){
  if(df$POS_END[i] == df$POS_START[i+1]){
    if(!in_progress){
      count <- count+1
      in_progress <- TRUE
    } else{
      in_progress <- FALSE
    }
  } else{
    in_progress <- FALSE
  }
}
