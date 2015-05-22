
smoothfunc = function(idx, seq, windowsize){
sum(seq[idx:(idx+windowsize-1)])/(windowsize)
}

extendsmooth = function(seq, windowsize=2){
idx = 1:length(seq)
smooth = sapply(idx, smoothfunc, seq=seq, windowsize=windowsize)
keep = !is.na(smooth)
smooth = smooth[keep]
additions = sum(!keep)
if (additions %% 2 == 0){
   smooth = c(rep(smooth[1],additions/2),smooth, rep(smooth[length(smooth)],additions/2))
   } else{
   smooth = c(smooth, smooth[length(smooth)])
   additions = additions-1
   smooth = c(rep(smooth[1],additions/2),smooth, rep(smooth[length(smooth)],additions/2))
   }
smooth
}
