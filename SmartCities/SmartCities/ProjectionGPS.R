
projection <- function(data,EPSG = 2154) {
        
    # loading libraries
    for (package in c('data.table','proj4','rgdal')) {
        if (!require(package, character.only=T, quietly=T)) {
            install.packages(package)
            library(package, character.only=T)
        }
    }
        
        tempsExec <- Sys.time()
        data <- read.csv(file = "../data/AllDatas_20151103_180302_NmeaTimeLatLongSatnbAtlHSpeedVSpeed_output.csv", header = FALSE, sep = ";")
        
        lat <- data[,3]
        lon <- data [,4]
        
        tempsExec1 <- Sys.time()
        # Creation du dataframe contenant les correspondances entre les codes EPSG et les parametre de projection
        EPSGPro4jMatching <- make_EPSG()
        # On remove les NAs
        EPSGPro4jMatching <- EPSGPro4jMatching[complete.cases(EPSGPro4jMatching),]
        # On regarde dans ce data frame a quelle ligne correspond le code epsg passe en parametre de fonction
        
        ligneParametreProjectionMatch <- match(EPSG, EPSGPro4jMatching$code)
        
        # Creation du dataframe contenant les correspondances entre les codes EPSG et les parametre de projection
        #EPSGPro4jMatching <- make_EPSG()
        # On remove les NAs
        #EPSGPro4jMatching <- EPSGPro4jMatching[complete.cases(EPSGPro4jMatching),]
        # On regarde dans ce data frame a quelle ligne correspond le code epsg passe en parametre de fonction
        #for(val in 1:nrow(EPSGPro4jMatching)){
                #if (EPSGPro4jMatching[val, "code"] == EPSG){
                        #ligne = val
                #}
        #}
        print(Sys.time()-tempsExec1)
        #On recupere le parametre de proj4 correspondant au code epsg passe en parmaetre
        paramProj4EFonctionEPSG <- EPSGPro4jMatching[ligneParametreProjectionMatch, 3]
        paramProj4EFonctionEPSG <- toString(paramProj4EFonctionEPSG)
        print(paramProj4EFonctionEPSG)

        # On projette les points GPS en coordonees selon le epsg passe en parametre
        conv <- proj4::project(xy = list(lat, lon),
                        proj = paramProj4EFonctionEPSG,
                        inverse = FALSE,
                        degrees = TRUE,
                        silent = FALSE,
                        ellps.default="sphere")

        XLambert <- conv$x
        YLambert <- conv$y

        # Definition du positionnement des capteurs

        # Pour les premiers points de chaque capteur, on definit son orientation par rapport a la deuxieme valeur
        ##### Definition de l'orientation #####
        arctan <- atan((XLambert[2]-XLambert[1])/(YLambert[2]-YLambert[1]))
        ##### Capteur 0 #####
        XLambertCapt0 <- XLambert[1]+0.337*(cos(arctan)+sin(arctan))
        YLambertCapt0 <- YLambert[1]+0.797*(-cos(arctan)+sin(arctan))
        ##### Capteur 3 #####
        XLambertCapt3 <- XLambert[1]+0.337*(cos(arctan)+sin(arctan))
        YLambertCapt3 <- YLambert[1]+0.15*(-cos(arctan)+sin(arctan))
        ##### Capteur 4 #####
        XLambertCapt4 <- XLambert[1]+0.337*(cos(arctan)+sin(arctan))
        YLambertCapt4 <- YLambert[1]-497*(-cos(arctan)+sin(arctan))
        #### Distance #####
        
        DistancePoints <- sqrt((XLambert[2]-XLambert[1])^2+(YLambert[2]-YLambert[1])^2)
        
        # Boucle pour les points de 2 a n
        for(i in 2:(length(XLambert)-1)) {
                
            ##### Definition de l'orientation #####
            arctan <- atan((XLambert[i]-XLambert[i-1])/(YLambert[i]-YLambert[i-1]))
            ##### Capteur 0 #####
            XLambertCapt0 <- c(XLambertCapt0, XLambert[i]+0.337*(cos(arctan)+sin(arctan)))
            YLambertCapt0 <- c(YLambertCapt0, YLambert[i]+0.797*(-cos(arctan)+sin(arctan)))
            ##### Capteur 3 #####
            XLambertCapt3 <- c(XLambertCapt3, XLambert[i]+0.337*(cos(arctan)+sin(arctan)))
            YLambertCapt3 <- c(YLambertCapt3, YLambert[i]+0.15*(-cos(arctan)+sin(arctan)))
            ##### Capteur 4 #####
            XLambertCapt4 <- c(XLambertCapt4, XLambert[i]+0.337*(cos(arctan)+sin(arctan)))
            YLambertCapt4 <- c(YLambertCapt4, YLambert[i]-0.497*(-cos(arctan)+sin(arctan)))
            #### Distance #####
            DistancePoints <- c(DistancePoints, sqrt((XLambert[i+1]-XLambert[i])^2+(YLambert[i+1]-YLambert[i])^2))
        }

        # Pour les derniers points de chaque capteur, on definit son orientation par rapport a la deuxieme valeur
        ##### Definition de l'orientation #####
        arctan <- atan((XLambert[length(XLambert)]-XLambert[length(XLambert)-1])/(YLambert[length(XLambert)]-YLambert[length(XLambert)-1]))
        ##### Capteur 0 #####
        XLambertCapt0 <- c(XLambertCapt0, XLambert[length(XLambert)]+0.337*(cos(arctan)+sin(arctan)))
        YLambertCapt0 <- c(YLambertCapt0, YLambert[length(XLambert)]+0.797*(-cos(arctan)+sin(arctan)))
        ##### Capteur 3 #####
        XLambertCapt3 <- c(XLambertCapt3, XLambert[length(XLambert)]+0.337*(cos(arctan)+sin(arctan)))
        YLambertCapt3 <- c(YLambertCapt3, YLambert[length(XLambert)]+0.15*(-cos(arctan)+sin(arctan)))
        ##### Capteur 4 #####
        XLambertCapt4 <- c(XLambertCapt4, XLambert[length(XLambert)]+0.337*(cos(arctan)+sin(arctan)))
        YLambertCapt4 <- c(YLambertCapt4, YLambert[length(XLambert)]-0.497*(-cos(arctan)+sin(arctan)))
        #### Distance #####
        DistancePoints <- c(DistancePoints, sqrt((XLambert[length(XLambert)]-XLambert[length(XLambert)-1])^2+(YLambert[length(XLambert)]-YLambert[length(XLambert)-1])^2))
        
        result <- data.table(timestamp = data[,1],
                             xCapt0 = XLambertCapt0,
                             yCapt0 = YLambertCapt0,
                             xCapt3 = XLambertCapt3,
                             yCapt3 = YLambertCapt3,
                             xCapt4 = XLambertCapt4, 
                             yCapt4 = YLambertCapt4,
                             distance = DistancePoints)
        print(head(result), digits = 10)

        print(Sys.time()-tempsExec)
}
