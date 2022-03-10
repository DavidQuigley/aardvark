#' Plot a reversion map showing the physical locations of the reversions in summary sums
#'
#' @param sums Summary of reversions
#' @param xlim x axis limits
#' @param height height of a single row; defaults 1
#' @param xinterval X axis interval between ticks in the label; defaults 250
#' @param exclude_blacklist boolean, if TRUE don't plot reversions marked with the evidence "overlaps_blacklist"; defaults TRUE
#' @param min_freq minimum number of observations to plot a candidate; defaults 2
plot_reversion_summary = function( sums, xlim, height=1,
                                   xinterval=250,
                                   exclude_blacklist=TRUE,
                                   min_freq=2){
    sr = sums$summary
    cigs = sums$cigar_ranges
    if( exclude_blacklist ){
        cigs = cigs[sr$evidence != "overlaps_blacklist"]
        sr = sr[sr$evidence != "overlaps_blacklist",]
    }
    cigs = cigs[ sr$N>=min_freq]
    sr = sr[ sr$N>=min_freq,]
    ord = order(sr$N, decreasing=FALSE)
    N_y = dim(sr)[1]
    plot( -1, -1, xlim=xlim,
          ylim=c( 0, N_y * height ),
          axes=FALSE,
          xlab="", ylab="", yaxs="i", xaxs="i")
    axis(1, seq( from=xlim[1], to=xlim[2], by=xinterval ), las=2 )

    cur_y = 1
    for(i in 1:N_y){
        idx = ord[i]
        if( i %% 2 == 0 ){
            rect( xlim[1], cur_y, xlim[2], (cur_y + height - 1), col="#eeeeee", border="#eeeeee")
        }
        for( rr in 1:dim( cigs[[idx]] )[1]){
            x1=cigs[[idx]]$ref_start[rr]
            x2=cigs[[idx]]$ref_end[rr]
            if( cigs[[idx]]$cigar_code[rr] == "D" ){
                rect( x1, cur_y, x2, (cur_y + height - 1), col="black", border="black")
            }else if( cigs[[idx]]$cigar_code[rr]  == "I" ){
                rect( x1, cur_y, x2, (cur_y + height - 1), col="#1f78b4",border="#1f78b4")
            }
        }
        cur_y = cur_y + height
    }
    box()
}

#' Barplot the frequency of each reversion within a sample
#'
#'
#' @param sums Summary of reversions
#' @param height height of each bar
#' @param xinterval interval of x axis ticks
#' @param exclude_blacklist exclude candidate reversions that intersect a blacklist region
#' @param min_freq minimum number of observations to plot a candidate, default 2
plot_reversion_frequency = function( sums, height=1, xinterval=10,
                                     exclude_blacklist=TRUE,
                                     min_freq=2){
    sr = sums$summary
    if( exclude_blacklist ){
        sr = sr[sr$evidence != "overlaps_blacklist",]
    }
    sr = sr[ sr$N>=min_freq,]
    ord = order(sr$N, decreasing=FALSE)
    N_y = dim(sr)[1]
    xlim=c(0, max(sr$N))
    plot( -1, -1, xlim=xlim,
          ylim=c( 0, N_y * height ),
          axes=FALSE,
          xlab="", ylab="", yaxs="i", xaxs="i")
    axis(1, seq( from=xlim[1], to=xlim[2], by=xinterval ) )
    cur_y = 1
    # 1f78b4
    for(i in 1:N_y){
        idx = ord[i]
        rect( 0, cur_y, sr$N[ idx ], (cur_y + height - 1), col="black", border="black")
        cur_y = cur_y + height
    }
    box()
}
