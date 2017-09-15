resf_vc		<- function( y, x = NULL, xconst = NULL, meig, method = "reml" ){

    lik_resf_vc		<- function( par0, ev, M, m, yy, n, nx, nsv, ne, emet ){
    	par		<- par0 ^ 2
    	evSqrt		<- NULL
    	for( i in 1:nsv ){
    		evv  	<- ev ^ par[ nsv + i ] * sum( ev ) / sum( ev ^ par[ nsv + i ] )
    		evSqrt	<- c( evSqrt, par[ i ] * sqrt( evv ) )
    	}
    	M[ -( 1:nx ), -( 1:nx ) ]	<- t(M[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
    	M[ -( 1:nx ),    1:nx   ]	<-   M[ -( 1:nx ),    1:nx   ] * evSqrt
    	M[    1:nx  , -( 1:nx ) ]	<- t(M[ -( 1:nx ),    1:nx   ] )
    	M0		<- M
    	diag( M [ -( 1:nx ), -( 1:nx ) ] ) <- diag( M[ -( 1:nx ), -( 1:nx ) ] ) + 1

    	m[-(1:nx)]	<- m[ -( 1:nx ) ] * evSqrt
    	test		<-try(Minv	<- solve( M, tol = 1e-25 ))
    	if(class(test)=="try-error"){
    		loglik  	<- Inf
    	} else {
    		b		<- Minv %*% m
    		sse		<- yy - 2 * t( b ) %*% m + t( b ) %*% M0 %*% b
    		dd		<- sse + sum( b[ -( 1:nx ) ] ^ 2 )
    		if( emet == "reml" ){
    			term1	<- determinant( M )$modulus
    			term2	<- ( n - nx ) * ( 1 + log( 2 * pi * dd / ( n - nx ) ) )
    		} else if( emet == "ml" ){
    			term1	<- determinant( as.matrix( M[ -( 1:nx ), -( 1:nx ) ] ) )$modulus
    			term2	<- n * ( 1 + log( 2 * pi * dd / n ) )
    		}
    		loglik		<- term1 + term2
    	}
    	return( loglik )
    }

    if( is.null( x ) ){
    	result	<- resf( y = y, x = xconst, meig = meig )
    	b_par	<- result$b
    	sf_par	<- result$v
    	e_stat	<- result$e
    	b_vc	<- NULL
    	bse_vc	<- NULL
    	bt_vc	<- NULL
    	bp_vc	<- NULL
    	pred	<- result$pred
    	resid	<- result$resid

    } else {
    	n     	<- length( y )
    	ev	<- meig$ev
    	X1	<- x
    	if( is.null( X1 ) == F ){
    		X1	<- as.matrix( X1 )
    		if( is.numeric( X1 ) == F ){
    			mode( X1 ) <- "numeric"
    		}
    		xind	<- apply( X1, 2, sd ) != 0
    		X1	<- as.matrix( X1[ , xind ] )
        	nx0	<- sum( xind )
    		if( nx0 == 0 ){
    			xname	<- NULL
    		} else {
    			xname	<- names( as.data.frame( X1 ) )
    		}
    	} else {
    		xname	<- NULL
    	}
        if( nx0 > 5 ) {
        	message( " Warning:" )
        	message( "  At most 4 or 5 covariates in x are prefarable in terns of stability" )
        	message( "  Please consider using xconst" )
        }

    	Xconst	<- xconst
    	if( is.null( xconst ) == F ){
    		Xconst	<- as.matrix( Xconst )
    		if( is.numeric( Xconst ) == F ){
    			mode( Xconst ) <- "numeric"
    		}
    		xfind	<- apply( Xconst, 2, sd ) != 0
    		Xconst	<- as.matrix( Xconst[ , xfind ] )
        	nxf	<- sum( xfind )
    		if( nxf == 0 ){
    			xfname	<- NULL
    		} else {
    			xfname	<- names( as.data.frame( Xconst ) )
    		}
    	} else {
        	xfind	<- 0
    		xfname	<- NULL
        	nxf	<- 0
    	}

        if( ( nx0 > 0 ) && ( nxf > 0 ) ){
        	for( dd  in 1:nx0  ){
        	for( ddf in 1:nxf ){
        		if( sum( X1[ , dd ] != Xconst[ , ddf ] ) == 0 ){
        			stop( " x and xconst cannot have the same column" )
        		}
        	}
        	}
        }

    	nsv	<- ifelse(is.null( X1 ), 1, dim( X1 )[ 2 ] + 1 )
    	X2	<- meig$sf
    	if( nsv >= 2 ){
    		for( i in 1:( nsv - 1 ) ){
   			X2   	<- cbind( X2, X1[, i ] * meig$sf )
    		}
    	}

    	X	<- as.matrix( cbind( 1, Xconst, X1, X2 ) )
    	ne	<- dim( X2 )[ 2 ]
    	nx	<- dim( X )[ 2 ] - ne
    	M   	<- crossprod( X )
    	m	<- crossprod( X, y )
    	yy     	<- sum( y ^ 2 )
    	par0	<- rep( 1, 2 * nsv )
    	res    	<- optim( fn = lik_resf_vc, par0, ev = ev, M = M, m = m, yy = yy,
    			   n = n, nx = nx, nsv = nsv, ne = ne, emet = method )
    	par2   	<- res$par ^ 2
    	loglik 	<- ( -1 / 2 ) * res$value

    	evSqrt	<- NULL
    	evSqrt2	<- NULL
    	for( i in 1:nsv ){
    		evv  	<- ev ^ par2[ nsv + i ] * sum( ev ) / sum( ev ^ par2[ nsv + i ] )
    		evSqrt	<- c( evSqrt, par2[ i ] * sqrt( evv ) )
    		evSqrt2	<- cbind( evSqrt2, par2[ i ] * sqrt( evv ) )
    	}

    	M[ -( 1:nx ), -( 1:nx ) ]	<- t( M[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
    	M[ -( 1:nx ),    1:nx   ]	<-    M[ -( 1:nx ),    1:nx   ] * evSqrt
    	M[    1:nx  , -( 1:nx ) ]	<- t( M[ -( 1:nx ),    1:nx   ] )
    	diag( M [ -( 1:nx ), -( 1:nx ) ] ) <- diag( M[ -( 1:nx ), -( 1:nx ) ] ) + 1
    	Minv		<- solve( M, tol = 1e-25 )
    	m[ -( 1:nx ) ]	<- m[ -( 1:nx ) ] * evSqrt
    	b		<- Minv %*% m
    	b[ -( 1:nx ) ]	<- b[ -( 1:nx ) ] * evSqrt
    	pred		<- X %*% b
    	resid		<- y - pred
    	SSE		<- sum( resid ^ 2 )
    	SSY		<- sum( ( y - mean( y ) ) ^ 2 )
    	sig		<- SSE / ( n - nx )
    	b_cov		<- sig * Minv
    	bse		<- sqrt( diag( b_cov ) )

    	nev		<- length( ev )
    	b_vc		<- matrix( 0, nrow = n, ncol = nsv )
    	bse_vc		<- matrix( 0, nrow = n, ncol = nsv )
    	for( j in 1:nsv ){
        	bid0		<- ifelse( j == 1, 1, nxf + j )
    		bid0_vc		<- ( nx + nev * ( j - 1 ) + 1 ):( nx + nev * j )
        	bid		<- c( bid0, bid0_vc )
        	b_vc  [ , j ]   <- b[ bid0 ] + meig$sf %*% b[ bid0_vc ]

    		b_cov_sub	<- b_cov[ bid, bid ]
        	sf2		<- t( t( meig$sf ) * evSqrt2[ , j ] )
        	x_sf		<- as.matrix( cbind( X[ , bid0 ], sf2 ) )
   		bse_vc[ , j ]	<- sqrt( colSums( t( x_sf ) * ( b_cov_sub %*% t( x_sf ) ) ) )
    	}

    	parV		<- par2[   1:nsv  ]
    	parR		<- par2[ -(1:nsv) ]

    	Xm		  <- X
    	Xm[ , -( 1:nx ) ] <- t( t( Xm[ , -( 1:nx ) ] ) * evSqrt )
    	np		<- nxf + nsv * 2 + 2
    	AIC		<- -2 * loglik + np * 2
    	BIC		<- -2 * loglik + np * log( n )
    	r2_0		<- 1 - SSE / SSY
    	r2		<- 1 - ( 1- r2_0 ) * ( n - 1 ) / ( n - np - 1 )

        if( nxf != 0 ) {
    		b		<- b  [ 2:( nxf + 1 ) ]
    		bse		<- bse[ 2:( nxf + 1 ) ]
    		bt		<- b / bse
    		df		<- sum( t( Xm ) * ( Minv %*% t( Xm ) ) )
        	bp		<- 2 - 2 * pt( abs( bt ), df = n - df )
    		b_par		<- data.frame( Estimate = b, SE = bse, t_value = bt, p_value = bp )
        	rownames( b_par )<- xfname
        } else {
    		df		<- sum( t( Xm ) * ( Minv %*% t( Xm ) ) )
        	b_par		<- NULL
        }

        bt_vc		<- b_vc / bse_vc
        bp_vc		<- 2 - 2 * pt( abs( bt_vc ), df = n - df )
        b_vc		<- data.frame( b_vc )
        bse_vc		<- data.frame( bse_vc )
        bt_vc		<- data.frame( bt_vc )
        bp_vc		<- data.frame( bp_vc )
        names( b_vc )	<- c( "(Intercept)", xname )
        names( bse_vc )	<- c( "(Intercept)", xname )
        names( bt_vc )	<- c( "(Intercept)", xname )
        names( bp_vc )	<- c( "(Intercept)", xname )
        parV		<- parV * sqrt( sig )
    	sf_par		<- data.frame( rbind( parV, parR ) )
    	names( sf_par )	<- c( "(Intercept)", xname )
    	rownames( sf_par )<- c( "shrink_sf_SE", "shrink_sf_alpha" )

    	e_stat		<- data.frame( stat = c( sqrt( sig ), r2, loglik, AIC, BIC ) )
    	if( method == "reml" ){
    		rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "rlogLik", "AIC", "BIC")
    		message( " RE-ESF model fit by REML" )
    	} else if( method == "ml" ){
    		rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "logLik", "AIC", "BIC" )
    		message( " RE-ESF model fit by ML" )
    	}
    }
    return( list( b = b_par, s = sf_par, e = e_stat,
    		  b_vc = b_vc, bse_vc = bse_vc, t_vc = bt_vc, p_vc = bp_vc, pred = pred, resid = resid ) )
}
