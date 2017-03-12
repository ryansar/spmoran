esf		<-function( y, x = NULL, vif = NULL, meig, fn = "r2" ){
    ObjEval	<-function( sfmod, fn ){
    	if( fn == "r2" ){
    		Obj	<- summary( sfmod )$adj.r.squared
    	} else if( fn == "aic" ){
    		Obj	<- -AIC( sfmod )
    	} else if( fn == "bic" ){
    		Obj	<- -BIC( sfmod )
    	}
    	return( Obj )
    }

    n     <- length( y )
    dum   <- 0
    if( is.null( x ) ){
    	sfmod	<- lm( y ~ 1 )
    	X	  <- NULL
    	nx	<- 1
    } else {
    	X00	<- as.matrix( x )
    	if( is.numeric( X00 ) == F ){
    		mode( X00 ) <- "numeric"
    	}
    	X	<- X00[ , apply( X00, 2, sd ) != 0 ]
    	if( dim( X )[ 2 ] == 0 ) {
    		sfmod	<- lm( y ~ 1 )
    		X	<- NULL
    		nx	<- 1
    	} else {
    		sfmod	<- lm( y ~ X )
    		nx	<- dim( X )[ 2 ] + 1
    		vif0	<- max( diag( solve( cor( X ) ) ) )
    		if( is.null( vif ) == F ){
    			if( vif0 >= vif ){
    			message( "Note:" )
    			message( paste( "  max(VIF) of the explanatory variables exceeds", vif ) )
    			message( "  No eigenvectors are selected" )
    			dd	<- data.frame( y, X )
    			dum	<- 1
    			}
    		}
    	}
    }

    sf      	<- as.matrix( meig$sf )
    ev      	<- meig$ev
    ne		<- length( ev )
    sf_list	<- 1:ne

    if( dum == 0 ){
    if( ( fn != "all" ) ){
    	Obj	<- ObjEval( sfmod = sfmod, fn = fn )
    	Obj_old	<- Obj - 1
    	obj	<- Obj
    	i	<- 1
    	while( ( Obj_old < Obj ) & ( i <= ne ) ){
    		Obj_old	<- Obj
    		Obj_list<- NULL
    		for( ii in 1:dim( sf )[ 2 ] ){
    			d	<- data.frame( cbind( y, X, sf[ , ii ] ) )
    			sfmod0	<- lm( y ~ ., d )
    			Obj	<- ObjEval( sfmod = sfmod0, fn = fn )
    			Obj_list<- c( Obj_list, Obj )
    		}
    		Obj	<- max( Obj_list )
    		if( Obj_old < Obj ){
    			obj	<- c( obj, Obj )
    			Obj_ind0<- ( 1:length( Obj_list ) ) * ( Obj_list == Obj )
    			Obj_ind	<- Obj_ind0[ Obj_ind0 != 0 ][ 1 ]
			sf_sel_l<- sf_list[    Obj_ind ]
    			sf_sel	<- sf     [ ,  Obj_ind ]
    			sf_list	<- sf_list[   -Obj_ind ]
    			sf	<- sf     [ , -Obj_ind ]
    			X0	<- cbind( X, sf_sel )
    			if( is.null(vif) == F ){
    				VIF	<- diag( solve( cor( X0 ) ) )
    				if( max( VIF ) < vif ){
    					X	<- data.frame( X0 )
    					names( X )[ dim( X )[ 2 ] ] <- paste( "sf", sf_sel_l, sep = "" )
    					dd	<- data.frame( y, X )
    					sfmod	<- lm( y ~ ., dd )
    				}
    			} else {
    				X	<- data.frame( X0 )
    				names( X )[ dim( X )[ 2 ] ] <- paste( "sf", sf_sel_l, sep = "" )
    				dd	<- data.frame( y, X )
    				sfmod	<- lm( y ~ ., dd )
    			}
    			i	<- i + 1
    		}
    		if( i %% 50 == 0 ) print( i )
    	}
    } else {
    	dd	<- data.frame( y, X )
    	for( ii in sf_list ){
    		X0	<- cbind( X, sf[, ii ] )
    		if( is.null( vif ) == F ){
    			VIF	<- diag( solve( cor( X0 ) ) )
    			if( max( VIF ) < vif ){
    				X	<- data.frame( X0 )
    				names( X )[ dim( X )[ 2 ] ] <- paste( "sf", ii, sep = "" )
    				dd	<- data.frame( y, X )
    			}
    		} else {
    			X	<- data.frame( X0 )
    			names( X )[ dim( X )[ 2 ] ] <- paste( "sf", ii, sep = "" )
    			dd	<- data.frame( y, X )
    		}
    	}
    	sfmod	<- lm( y ~ ., dd )
    }
    }
    pred	<- predict( sfmod )
    resid	<- residuals ( sfmod )
    df		<- dim( X )[ 2 ] + 1
    sig		<- sum( resid ^ 2 ) / ( n - df )
    if( nx == 1 ){
    	b_par	<- summary( sfmod )$coefficients[   1:nx , 1:4 ]
    	names( b_par ) <- c( "(Intercept)", "SE","t_value", "p_value" )
    } else {
    	b_par	<- data.frame( summary( sfmod )$coefficients[   1:nx , 1:4 ] )
    	names( b_par ) <- c( "Estimate", "SE","t_value", "p_value" )
    }

    if( dim( dd )[ 2 ] == ( nx + 1 ) ){
    	r	<- NULL
    	nr	<- 0
    } else {
    	r	<- data.frame( summary( sfmod )$coefficients[ -(1:nx), 1:4 ] )
    	nr	<- dim( r )[ 1 ]
    	names( r )	<- c( "Estimate", "SE","t_value", "p_value" )
    }

    if( is.null( r ) ){
    	SF	<- NULL
    } else {
    	if( nx == 1 ){
    		SF	<- as.matrix( X ) %*% c(r[,1])
    	} else {
    		SF	<- as.matrix( X[ , -( 1:( nx - 1 ) ) ] ) %*% c(r[,1])
    	}
    }
    vif 	<- data.frame(diag( solve( cor( X ) ) ) )
    names(vif)	<- "VIF"

    r2		<- summary( sfmod )$adj.r.squared
    loglik	<- logLik( sfmod )
    AIC		<- AIC( sfmod )
    BIC		<- BIC( sfmod )
    e_stat	<- data.frame( stat = c( sqrt( sig ), r2, loglik, AIC, BIC ) )
    rownames( e_stat ) <- c( "resid_SE", "adjR2", "logLik", "AIC", "BIC")
    if( dum ==0) {
    message( paste( "  ", nr, "/", ne, " eigenvectors are selected", sep = "" ) )
    }
    return( list( b = b_par, r = r, vif = vif, e = e_stat, sf = SF, pred = pred, resid = resid ) )
}






