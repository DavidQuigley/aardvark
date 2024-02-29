#' Create new hash
#'
hsh_new = function(){
    new.env(hash=TRUE, parent=emptyenv())
}

#' Test whether key exists in hash H
#' @param H hash object
#' @param key key to check
hsh_in = function(H, key){
    if( key == "" ){
        FALSE
    }else{
        exists(key, H)
    }
}

#' Return value of key in hash H
#' @param H hash object
#' @param key key to set
#' @param na.if.not.found If key not in H and TRUE, return NA; if FALSE, raise an error
hsh_get = function( H, key, na.if.not.found = FALSE ){
    if( length(key)==1 ){
        if( na.if.not.found ){
            if( exists(key, H) )
                get(key, H)
            else
                NA
        }
        else{
            get(key, H)
        }
    }
    else{
        results = rep(0, length(key) )
        if( !na.if.not.found ){
            for(i in 1:length(key) ){
                if( exists(key[i], H) ){
                    results[i] = get(key[i], H )
                }
                else{
                    results[i] = NA
                }
            }
        }
        else{
            for(i in 1:length(key) ){
                results[i] = get(key[i], H )
            }
        }
        results
    }
}

#' Set key to value in hash H
#' @param H hash object
#' @param key key to set
#' @param value new value of key in H
hsh_set = function( H, key, value ){
    assign(key, value, envir=H)
}


#' split vector all at once
#'
#' @param v vector
#' @param string, split value
#' @param col column to return
#' @param last return last
#' @param first return first
get.split.col = function(v, string, col=0, last=F, first=F){
    if( last & first )
        stop("Cannot request both last and first column")
    if( col==0 & !last & !first)
        stop("Must request either a column by index, first, or last")

    for(i in 1:length(v)){
        x = strsplit( v[i], string, fixed=T)[[1]]
        if(last){
            v[i] = x[length(x)]
        }
        else if(first){
            v[i] = x[1]
        }
        else{
            v[i] = x[col]
        }
    }
    v
}



#' return indices into A and B restricted to perfect matches where idx.A == idx.B for each i in matched pairs
#'
#' @param A vector
#' @param B vector
#' @param allow.multiple.B If true, A->B can be one to many; otherwise return first match
match.idx=function(A, B, allow.multiple.B = FALSE){
    if( allow.multiple.B ){
        idx.B = which(B %in% A)
        idx.A = match(B[idx.B], A)
    }
    else{
        in.both = intersect(A,B)
        idx.A = match(in.both, A)
        idx.B = match(in.both, B)
    }
    C= data.frame(idx.A, idx.B)
    if( sum( A[ C$idx.A ] != B[ C$idx.B] )>0 )
        stop("ERROR! At least one in idx.A not the same as matched item in idx.B")
    C
}


#' Returns reverse complement of a vector of A,C,G,T
#' @param ts_fwd character string of transcript to generate the reverse complement
#' @return reverse complement vector of ts_fwd
#'
reverse_complement = function( ts_fwd ){
    ts_rev = vector(mode="character", length=length( ts_fwd ) )
    ctr=1
    for(i in seq(from=length( ts_fwd ), to=1, by=-1)){
        if( ts_fwd[i]=="A" ){
            ts_rev[ctr] = "T"
        }else if( ts_fwd[i]=="C" ){
            ts_rev[ctr] = "G"
        }else if( ts_fwd[i]=="G" ){
            ts_rev[ctr] = "C"
        }else{
            ts_rev[ctr] = "A"}
        ctr=ctr+1
    }
    ts_rev
}

#' Remove an arbitrary row from a dataframe and return the result
#
#' @param df dataframe to modify
#' @param row_id_to_remove row to remove
remove_row = function( df, row_id_to_remove ){
    if( row_id_to_remove > dim(df)[1]){
        stop( "row_id_to_remove is larger than dimension of df")
    }
    df[ setdiff( 1:dim(df)[1], row_id_to_remove), ]
}

#' replace an existing row with the contents of df_new (may be multiple rows)
#'
#' @param df dataframe to modify
#' @param row_id_to_replace row to replace
#' @param df_new new dataframe, columns must match df
replace_row = function( df, row_id_to_replace, df_new ){
    if( dim(df)[1] == 1 ){
        df_new
    }else if( row_id_to_replace == 1) {
        rbind( df_new, df[2: dim(df)[1],] )
    }else if( row_id_to_replace == dim(df)[1] ){
        rbind( df[1:(dim(df)[1]-1),], df_new )
    }else{
        rbind( df[ 1 : (row_id_to_replace-1),],
               df_new,
               df[ (row_id_to_replace+1) : dim(df)[1] , ])
    }
}


#' convert VCF-style REF and ALT to equal length with dashes for indels
#'
#' e.g. a REF,ALT pair TTCTT,T -> TTCTT,T----
#' e.g. a REF,ALT pair T,TTCTT -> T----, TTCTT
#' e.g. a REF,ALT pair TCT, TAG -> TCT,TAG
#'
#' @param s_ref reference nucleotide string
#' @param s_alt alternate nucleotide string
#' @export
#' @returns list with values REF, ALT for new formatting
convert.VCF.REFALT.to.dash.format = function( s_ref, s_alt ){
    v_ref = strsplit(s_ref,"")[[1]]
    v_alt = strsplit(s_alt,"")[[1]]
    n = max( c( length( v_ref ), length( v_alt ) ) )
    v_ref_final = rep( "-",  n )
    v_alt_final = rep( "-",  n )
    for( i in 1:length( v_ref )){
        v_ref_final[i] = v_ref[i]
    }
    for( i in 1:length( v_alt )){
        v_alt_final[i] = v_alt[i]
    }
    list( REF = paste( v_ref_final, collapse=""), ALT = paste(v_alt_final, collapse="") )
}

#' Given a read, return a simple data frame that write the reference sequence
#' where there is a perfect match and D or I for deletion/insertion.
#' @param read aardvark::Read object under consideration
#' @param alignment_window aardvark::AlignmentWindow object that holds reference sequence search space for realignment
read_to_genome_sequence = function( read, alignment_window ){

    # replaced because last (or first) value can be NA if softclipped
    # nt_width = read$cigar_ranges$ref_end[ dim(read$cigar_ranges)[1] ] - read$cigar_ranges$ref_start[ 1 ]
    ref_max = max( read$cigar_ranges$ref_end, na.rm=TRUE )
    ref_min = min( read$cigar_ranges$ref_start, na.rm=TRUE )
    nt_width = ref_max - ref_min

    if( is.na( ref_max ) | is.na(ref_min) ){
        data.frame(pos=c(), nt=c())
    }else{
        nt = rep(".", nt_width )
        pos = rep(0, nt_width )
        for( i in 1:dim(read$cigar_ranges)[1]){
            if( read$cigar_ranges$cigar_code[i] != "S" ){
                idx = (1 + read$cigar_ranges$ref_start[i] - read$pos) : (1 + read$cigar_ranges$ref_end[i] - read$pos )
                positions = read$cigar_ranges$ref_start[i] : read$cigar_ranges$ref_end[i]
                if( length(positions) != length(idx ) ){
                    print( "debug me")
                }
                pos[ idx ] = positions
                if( read$cigar_ranges$cigar_code[i]=="M" ){
                    nts = strsplit( as.character( AW_seq(alignment_window, read$cigar_ranges$ref_start[i], read$cigar_ranges$ref_end[i] ) ), "")[[1]]
                    if( length( nts ) != length( idx ) ){
                        print( "debug me" )
                    }
                    nt[ idx ] = nts
                }else if( read$cigar_ranges$cigar_code[i]=="D" ){
                    nt[ idx ] = "-"
                }else if( read$cigar_ranges$cigar_code[i]=="I" ){
                    nt[ idx ] = "I"
                }
            }
        }
    }
    data.frame( pos, nt )
}

