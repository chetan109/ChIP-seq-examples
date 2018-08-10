readData = function( dir, fnpOHT, fnmOHT = NULL, type = "bw", as = "RleList" ){
    rv = list();
    message( "pOHT" );
    if( type == "bw" ){
        rv[["pOHT"]] = import.bw( paste0( dir, fnpOHT ), as = as );
        if( !is.null( fnmOHT ) ){
            message( "mOHT" );
            rv[["mOHT"]] = import.bw( paste0( dir, fnmOHT ), as = as );
        }
    }else if( type == "wig" ){
        rv[["pOHT"]] = import.wig( paste0( dir, fnpOHT ), as = as );
        if( !is.null( fnmOHT ) ){
            message( "mOHT" );
            rv[["mOHT"]] = import.wig( paste0( dir, fnmOHT ), as = as  );	
        }
    }else{
        stop( "ERROR : unknown type ", type );
    }
    return( rv );
}

computeProfile = function( bed, wig, w = 20000, span = 200, seqlens ,method="mean"){
    if( class( wig ) != "SimpleRleList" ){
        stop( "ERROR : unknown class of wig, please provide a SimpleRleList" );
    }
    mat = NULL;
    for( i in 1:length( bed ) ){
        message( i, "/", length( bed ) );
        bedi = bed[i, ];
        chr = as.character( seqnames( bedi ) );
        cov = wig[[chr]];
        
        
        center = start( bedi ) + 4;
        stW = center - w;
        edW = center + w;
        
        if( span == 1 ){
            vm = as.numeric( Views( cov, start = stW, end = edW )[[1]] )
        }else{
            sts = seq( stW, edW - span + 1, span );
            eds = seq( stW + span - 1, edW, span );
            v = Views( cov, start = sts, end = eds );
            if(method =="sum"){
                vm =  sum( v );    
            }else {
                vm =  mean( v );
            }
            
            vm[sts >= seqlens[chr] | sts > length( cov )] = 0;
            vm[sts < 1] = 0;
        }
        mat = rbind( mat, vm );
    }
    #rv = colMeans( mat );
    return( mat );
}

compute1ValPerSite = function( bed, wig, w = 20000, seqlens, fun = "sum" ){
    if( class( wig ) != "SimpleRleList" ){
        stop( "ERROR : unknown class of wig, please provide a SimpleRleList" );
    }
    vec = NULL;
    for( i in 1:length( bed ) ){
        if( i %% 100 == 0 ){
            message( i, "/", length( bed ) );
        }
        bedi = bed[i, ];
        chr = as.character( seqnames( bedi ) );
        cov = wig[[chr]];
        center = start( bedi ) + 4;
        stW = center - w;
        edW = center + w;
        stW[stW < 1] = 1;
        edW[edW > seqlens[chr] | edW > length( cov )] = min( length( cov), seqlens[chr] );
        v = Views( cov, start = stW, end = edW );
        if( fun == "sum" ){
            vm = sum( v );
        }else if( fun == "mean" ){
            vm =  mean( v );
        }else{
            stop( "ERROR : unknown function fun = ", fun, " - fun must be 'sum' or 'mean'" );
        }
        vec = c( vec, vm );
    }
    return( vec );
}


computeProfileLog = function( bed, wig, w = 2500, span = 1, seqlens ){
    if( class( wig ) != "SimpleRleList" ){
        stop( "ERROR : unknown class of wig, please provide a SimpleRleList" );
    }
    mat = NULL;
    for( i in 1:length( bed ) ){
        # message( i, "/", length( bed ) );
        bedi = bed[i, ];
        chr = as.character( seqnames( bedi ) );
        cov = wig[[chr]];
        stW = start( bedi ) - w;
        edW = end( bedi ) + w;
        if( span == 1 ){
            stW = start( bedi ) + 3 - w;
            edW = end( bedi ) - 3 + w;
            if( stW > length( cov ) ){
                vm = rep( NA, ( w * 2 ) + 1 );
            }else if( edW > length( cov ) ){
                edW = length( cov );
                vm = c( as.numeric( Views( cov, start = stW, end = edW )[[1]] ), rep( 0, edW - length( cov ) ) )
            }else{
                vm = as.numeric( Views( cov, start = stW, end = edW )[[1]] )
            }
        }else if( span < 8 ){
            stW = start( bedi ) + 3 - w;
            edW = end( bedi ) - 3 + w;
            sts = seq( stW, edW - span + 1, span );
            eds = seq( stW + span - 1, edW, span );
            v = Views( cov, start = sts, end = eds );
            vm =  mean( v );
            vm[sts >= seqlens[chr] | sts > length( cov )] = 0;
            vm[sts < 1] = 0;
        }else{
            stW = start( bedi ) - w;
            edW = end( bedi ) + w;
            sts = seq( stW, edW - span + 1, span );
            eds = seq( stW + span - 1, edW, span );
            v = Views( cov, start = sts, end = eds );
            vm =  mean( v );
            vm[sts >= seqlens[chr] | sts > length( cov )] = 0;
            vm[sts < 1] = 0;
        }
        mat = rbind( mat, vm );
    }
    rv = colMeans( mat );
    return( rv );
}

RenderColor <- function(pval,datalog,inter=0.05,maxi=0.01){
    dicolor = list("neutral" = "#AFABAB", 
                   "sign_pos"="#CC6600",
                   "very_sign_pos"="#C00000",
                   "sign_neg"="#9DC3E6",
                   "very_sign_neg"="#0070C0")
    if(pval >inter){
        res = dicolor[["neutral"]]
    }else {
        if(pval > maxi){
            type = "sign"
        }else {
            type = "very_sign"
        }
        if(datalog > 0){
            res = dicolor[[paste0(type,"_pos")]]
        }else if(datalog < 0){
            res = dicolor[[paste0(type,"_neg")]]
        }else {
            res = dicolor[["neutral"]]
        }
    }
    return(res)
}


computeProfileGenes = function( bed, wig, w = 3000, span = 200, seqlens , fun="sum",get="one"){
    if( class( wig ) != "SimpleRleList" ){
        stop( "ERROR : unknown class of wig, please provide a SimpleRleList" );
    }
    mat = NULL;
    for( i in 1:length( bed ) ){
        if( i %% 100 == 0 ){
            message( i, "/", length( bed ) );
        }
        bedi = bed[i, ];
        len = end(bedi) - start(bedi);
        if( len < 100 ){
            message( "Warning : Gene too small - n° ", i, " - ", bedi$name );
        }else{
            chr = as.character( seqnames( bedi ) );
            if(!chr %in% names(wig)){
                next;
            }
            cov = wig[[chr]];
            stW_st = start(bedi) - w;
            edW_st = start(bedi);
            stW_ed = end(bedi);
            edW_ed = end(bedi) + w;
            sts_st = seq( stW_st, edW_st - span + 1, span );
            eds_st = seq( stW_st + span - 1, edW_st, span );
            sts_ed = seq( stW_ed, edW_ed - span + 1, span );
            eds_ed = seq( stW_ed + span - 1, edW_ed, span );
            div = len / 100;
            ent = trunc( div );
            dec = div - ent;
            if( dec > 0.75 ){
                span_c = ent + 1;
            }else{
                span_c = ent;
            }
            sts_c = seq( start(bedi), end(bedi) - 1, span_c );
            eds_c = seq( start(bedi) + span_c - 1, end(bedi) + span_c - 1, span_c );
            if( length( sts_c ) > 100 ){
                sts_c = sts_c[-c( 101:length( sts_c ) )];
            }
            if( length( eds_c ) > 100 ){
                eds_c = eds_c[-c( 101:length( eds_c ) )];
            }
            if( length( sts_c ) < 100 | length( eds_c ) < 100 ){
                message( "Warning : Less than 100 regions - n° ", i, " - ", bedi$name);
            }else{
                eds_c[100] = end(bedi);
                v_st = Views( cov, start = sts_st, end = eds_st );
                vm_st =  switch(fun,
                                "sum"=sum(v_st),
                                "mean"=mean( v_st ));
                vm_st[sts_st >= seqlens[chr] | sts_st > length( cov )] = 0;
                vm_st[sts_st < 1] = 0;
                v_ed = Views( cov, start = sts_ed, end = eds_ed );
                vm_ed =  switch(fun,
                                "sum"=sum(v_ed),
                                "mean"=mean( v_ed ));
                vm_ed[sts_ed >= seqlens[chr] | sts_ed > length( cov )] = 0;
                vm_ed[sts_ed < 1] = 0;
                v_c = Views( cov, start = sts_c, end = eds_c );
                vm_c =  switch(fun,
                               "sum"=sum(v_c),
                               "mean"=mean( v_c ));
                vm_c[sts_c >= seqlens[chr] | sts_c > length( cov )] = 0;
                vm_c[sts_c < 1] = 0;
                if( as.character( strand(bedi) ) == "+" ){
                    vec = c( vm_st, vm_c, vm_ed );
                }else{
                    vec = c( rev( vm_ed ), rev( vm_c ), rev( vm_st ) );
                }	
                if( length( vec ) == ( 100 + ( 2 * ( w / span ) ) ) ){
                    mat = rbind( mat, vec );
                }else{
                    stop( "ERROR : profiles of different size" );
                }	
            }
        }
    }
    if(get == "one"){
        rv = switch(fun,
                    "sum"=colSums( mat ),
                    "mean"=colMeans( mat ));
    }else{
        rv = mat
    }
    return( rv );
}


compute1ValPerSiteForBLESS <- function( bed, wig, w = 20000, seqlens, fun = "sum" ){
    if( class( wig ) != "SimpleRleList" ){
        stop( "ERROR : unknown class of wig, please provide a SimpleRleList" );
    }
    vec = NULL;
    for( i in 1:length( bed ) ){
        if( i %% 100 == 0 ){
            message( i, "/", length( bed ) );
        }
        bedi = bed[i, ];
        chr = as.character( seqnames( bedi ) );
        cov = wig[[chr]];
        
        stW = start(bedi) - w;
        edW = end(bedi) + w;
        stW[stW < 1] = 1;
        edW[edW > seqlens[chr] | edW > length( cov )] = min( length( cov), seqlens[chr] );
        v = Views( cov, start = stW, end = edW );
        if( fun == "sum" ){
            vm = sum( v );
        }else if( fun == "mean" ){
            vm =  mean( v );
        }else{
            stop( "ERROR : unknown function fun = ", fun, " - fun must be 'sum' or 'mean'" );
        }
        vec = c( vec, vm );
    }
    return( vec );
}