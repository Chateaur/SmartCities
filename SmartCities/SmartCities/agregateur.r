##### Probleme ????? #####
#
# les donnees du fichier gps sont-elles bien :
#TS1-TS2
#TS2-TS3
#TS3-TS4 ...
#
# Dans la fonction lookRange, la for va de 1 a length(lux) ?? Pareil pour le code appelant pseudotimestamp
#
# a <- dataTable[-1,] dans le subsetDataTable et pseudotimestamp??
#
#rbindlist bind bien dans le bon ordre ?
#
# Dans la fonction loop : gps[(i*2)] ??

# loading libraries
for (package in c('data.table', 'bit64','parallel','pbapply','plyr','dplyr','rmarkdown','proj4','rgdal')) {
        if (!require(package, character.only=T, quietly=T)) {
                install.packages(package)
                library(package, character.only=T)
        }
}

#Fonction d'execution principale du programme
start <- function(){

        #TimeStamp du temps d'execution total du script
        t0 <- Sys.time()

        #Chargement du fichier GPS
        gps <- loadGps()

        #Chargement des fichiers de Luxmetres
        lux <- loadLux()
        
        print(head(lux), digits = 15)

        #Cluster
        ta <- Sys.time()
        cat("*Creation des Threads et du cluster\n")
        cl <- makeCluster(detectCores())
        clusterEvalQ(cl,library(bit64))
        clusterEvalQ(cl,library(data.table))
        clusterEvalQ(cl,library(proj4))
        clusterEvalQ(cl,library(rgdal))
        cat("** DONE\n")
        print(Sys.time() - ta)

        #Projection des coordonnees GPS en positions [X:Y]
        taa <- Sys.time()
        cat("* Lancement de la projection des points\n")
        listGps <- subSetDataTable(gps)
        projectedGps <- parLapply(cl,listGps,projection,2154)
        projectedGps <- rbindlist(projectedGps)
        cat("** Done\n")
        print(Sys.time() - taa)

        #Generation des pseudoTimeStamps & lissage des moyennes
        tc <- Sys.time()
        cat("Generation des pseudoTimeStamps & lissage des moyennes\n")
        tempPseudoTimeStampDT <- parLapply(cl,lux,pseudoTimestamp)
        for (i in 1:length(tempPseudoTimeStampDT)) {
                tempPseudoTimeStampDT[[i]][, CAPT := i]
        }
        cat("** DONE\n")
        print(Sys.time() - tc)

        #Definition des constantes pour la calibration des luxmetres
        tb <- Sys.time()
        cat("*Generation des constantes de calibration pour la conversion en lux & conversion\n")
        droite <- parLapply(cl,lux,convertLux,lMax = 3000,uMax = 3.3)

        #Conversion Volts en Lux
        for (i in 1:length(tempPseudoTimeStampDT)) {
                tempPseudoTimeStampDT[[i]] <- cbind(tempPseudoTimeStampDT[[i]],data.frame("LUX" = droite[[i]]$a*tempPseudoTimeStampDT[[i]]$V + droite[[i]]$b))
        }
        cat("** DONE\n")
        print(Sys.time() - tb)

        #Capteurs
        pseudoTimeStampDT <- rbindlist(tempPseudoTimeStampDT)
        keycols <- c("TS", "CAPT")
        pseudoTimeStampDT <- setkeyv(pseudoTimeStampDT, keycols)

        #Script principal
        te <- Sys.time()

        cat("* Agregation et ecriture des donnees\n")

        index <- subSetDataTable(projectedGps)
        parLapply(cl,index,loop,pseudoTimeStampDT)

        stopCluster(cl)
        cat("** DONE\n")
        print(Sys.time() - te)

        #Enregistrement des données dans le fichier
        tf <- Sys.time()
        cat("* Merging des fichiers de donnees\n")
        index <- createIndex()
        for (i in 2:length(index$output)) {
                file.append(index$output[1],index$output[i])
                file.remove(index$output[i])
        }
        cat("** DONE\n")
        print(Sys.time() - tf)

        #Fin
        cat("******************************************\n")
        cat("**  Total Execution time of the script : **\n")
        print(Sys.time() - t0)
}


#Fonction de loop principale
loop <- function(gps,lux){

        ##Definition de cette fonction dans la loop
        #Retourne une matrice contenant les points positionn?s sur une droite reliant deux points avec la valeur luxmetre correspondante
        getDroite <- function(x1,x2,y1,y2,subsetCapt){
                if (nrow(subsetCapt) >= 1) {

                        #Calcul du nombre de points
                        nbPoints <- length(subsetCapt)

                        # Calcul de l'ecart des points sur l'axe des abscisses et des ordonn?es
                        ecartPointsX <- (x2 - x1)/nbPoints
                        ecartPointsY <- (y2 - y1)/nbPoints

                        i <- c(0:nbPoints)

                        xTemp <- i*ecartPointsX + x1
                        yTemp <- i*ecartPointsY + y1

                        res <- data.table (xTemp,yTemp,subsetCapt$LUX)
                }
        }

        #Cherche dans chaque data.table les valeurs de capteurs
        #correspondant au range passe en parametre
        lookRange <- function(gps,lux){

                TS1 <- gps[1]
                TS2 <- gps[9]

                class(TS1) <- "integer64"
                class(TS2) <- "integer64"

                min <- lux[J(TS1), roll=-Inf, which = TRUE, mult = "first", rollends = TRUE]
                max <- lux[J(TS2), roll=Inf, which = TRUE, mult = "last", rollends = TRUE]

                if(min != 1 && max!= 1){

                        temp <- lux[min:max]

                        #Subset des valeurs de capteurs
                        droiteCaptList <- list()
                        for (i in 1:3) {
                                #Calcul des points interm?diaires en leur associant leur valeur de luxm?tre
                                droiteCaptList[[i]] <- getDroite(gps[(i*2)],gps[(8+(i*2))],gps[((i*2)+1)],gps[(8+(i*2))],subset(temp, CAPT == i))
                        }

                        finalToWrite <- rbindlist(droiteCaptList)
                }
        }

        data <- apply(gps,MARGIN = 1,lookRange,lux)
        write.table(rbindlist(data),paste("../data/output_",sample(1:100, 1),".csv",sep = ""), append = FALSE, quote = FALSE, sep = ";",row.names=FALSE)
}

#Fonction qui parcours et applique en parallele les pseudo timestamps aux fichiers de donnees
pseudoTimestamp <- function(lux){

        a <- lux
        a <- a[-1,]
        a <- rbind(a,data.frame(a[length(a$V1)]))
        index <- cbind(lux,a)
        colnames(index) <- c("T1","a","b","CAPT","T2","c","d","CAPT")

        #Generation du pseudoTimestamp
        TempDtPseudoTimestamp <- list()

        TempDtPseudoTimestamp[[1]] <- data.table(lux$V1)
        TempDtPseudoTimestamp[[2]] <- data.table(as.integer64((index$T1 + index$T2)/2))

        #Agregation des timestamps + pseudoTimestamps
        dtPseudoTimestamp <- data.table()
        dtPseudoTimestamp <- rbindlist(TempDtPseudoTimestamp)

        #Lissage des valeurs
        lissageValue <- list()
        parite <- c(0:(length(index$a)-1))*2+1
        lissageValueEven <- data.table((index$a + index$b + index$c)/3)
        lissageValue[[2]] <- cbind(lissageValueEven,parite)

        parite <- c(0:(length(index$a)-1))*2
        lissageValueOdd <- data.table((index$b + index$c + index$d)/3)
        lissageValue[[1]] <- cbind(lissageValueOdd,parite)

        lissage <- rbindlist(lissageValue)
        lissage <- cbind(dtPseudoTimestamp,lissage)
        lissage <- lissage[order(V1),]
        lissage$parite <- NULL
        colnames(lissage) <- c("TS","V")

        lissage
}

##Fonction qui permet de calculer Uo dans un data.table de capteur
##Pour cela, on calcule la moyenne des valeurs pendant la premiere minute
convertLux <- function (captData,lMax,uMax) {
        ##On considere que une minute correspond a un ecart de 60 000 000 en TS
        setkeyv(captData, c("V1"))

        #On repere le premier timeStamp du fichier
        TS1 <- captData[[1,1]]

        #On ajoute une minute ? ce timeStamp
        TS2 <- TS1 + 60000000

        #On subset la data.table pour obtenir la moyenne des valeurs de capteur pendant
        #la premiere minute
        min <- 1
        max <- captData[J(TS2), roll=Inf, which = TRUE, mult = "last", rollends = TRUE]

        temp <- captData[min:max]

        #On obtient notre valeur de calibration
        u0 <- mean(temp$V2)
        a <- lMax/(uMax - u0)
        b <- -a*u0
        retour <- data.frame("u0" = u0,"a" = a,"b" = b)
        retour
}



#Fonction qui permet de substter des data.table afin de les preparer pour des traitements multi-Threads Renvois une liste contenant morceau par morceau le maximum de data par Coeur Processeurs
subSetDataTable <- function(dataTable){
        a <- dataTable[-1,]
        a <- rbind(a,dataTable[nrow(dataTable)])
        index <- cbind(dataTable,a)
        colnames(index)[1] = "V1"

        #Calcul de la longueur de chaque Thread
        lengthThread <- ceiling(length(index$V1)/(detectCores()))

        tempList <- list()
        #Subset et append des data.tables
        for (i in 1:(detectCores())) {
                tempList[[i]] <- subset(index, (V1 >= index[lengthThread*(i-1) +1]$V1) & (V1 < index[lengthThread*i+1]$V1))
        }

        tempList
}

# Fonction qui va permettre d'identifier et d'indexer les fichiers de donnees
createIndex <- function(){

        setwd("../data")
        absolutePath <- getwd()

        tempListFiles <- list.files()
        splitedTempListFiles <- strsplit(tempListFiles,"_")
        trueDataList <- list()

        for(i in 1:length(list.files())){
                if(splitedTempListFiles[[i]][1] == "AllDatas"){
                        if(splitedTempListFiles[[i]][4] == "CalibrationMechanism"){
                                trueDataList$calibration <- c(trueDataList$calibration,paste(absolutePath,"/",tempListFiles[i],sep=""))
                        }
                        else if(splitedTempListFiles[[i]][4] == "temperature"){
                                trueDataList$temperature <- c(trueDataList$temperature,paste(absolutePath,"/",tempListFiles[i],sep=""))
                        }
                        else if(splitedTempListFiles[[i]][4] == "NmeaTimeLatLongSatnbAtlHSpeedVSpeed"){
                                trueDataList$gps <- c(trueDataList$gps,paste(absolutePath,"/",tempListFiles[i],sep=""))
                        }
                        else if(grep("ch",splitedTempListFiles[[i]][4]) > 0) {
                                trueDataList$lux <- c(trueDataList$lux,paste(absolutePath,"/",tempListFiles[[i]],sep=""))
                        }
                }
                else if(splitedTempListFiles[[i]][1] == "output"){
                        trueDataList$output <- c(trueDataList$output,paste(absolutePath,"/",tempListFiles[[i]],sep=""))
                }

        }

        setwd("../smartcities")
        trueDataList
}

# fonction qui renvoie une liste de luxmetres
loadLux <- function(){
        t1 <- Sys.time()
        cat("* Ouverture des fichiers de Luxmetres\n")
        cl <- makeCluster(detectCores())

        index <- createIndex()
        lux <- parLapply(cl,index$lux,fread)

        for (i in 1:length(lux)) {
                lux[[i]][, CAPT := i]
        }

        stopCluster(cl)
        cat("** DONE\n")
        print(Sys.time() - t1)
        lux
}

# fonction qui renvoie la data.table GPS
loadGps <- function(){
        t1 <- Sys.time()
        cat("* Ouverture du fichier GPS\n")

        index <- createIndex()
        gps <- fread(index$gps[1])

        cat("** DONE\n")
        print(Sys.time() - t1)
        gps
}

# Fonction permettant de projetter la latitude et la longitude en coordonnees GPS selon l'EPSG passe en parametre
projection <- function(data,EPSG = 2154) {

    #Recuperation des donnees utiles du fichier GPS
    TS  <- data$V1
    lat <- data$V3
    lon <- data$V4

    # Creation du dataframe contenant les correspondances entre les codes EPSG et les parametre de projection
    EPSGPro4jMatching <- make_EPSG()
    # On remove les NAs
    EPSGPro4jMatching <- EPSGPro4jMatching[complete.cases(EPSGPro4jMatching),]
    # On regarde dans ce data frame a quelle ligne correspond le code epsg passe en parametre de fonction
    ligneParametreProjectionMatch <- match(EPSG, EPSGPro4jMatching$code)

    #On recupere le parametre de proj4 correspondant au code epsg passe en parmaetre
    paramProj4EFonctionEPSG <- as.character(EPSGPro4jMatching[ligneParametreProjectionMatch,3])

    # On projette les points GPS en coordonees selon le epsg passe en parametre
    conv <- proj4::project(xy = list(lat, lon),
                           proj = paramProj4EFonctionEPSG,
                           inverse = FALSE,
                           degrees = TRUE,
                           silent = FALSE,
                           ellps.default="sphere")
    #Conversion en Data Table
    dtConv <- as.data.table(conv)

    XLambert <- dtConv$x
    YLambert <- dtConv$y

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

    # Pour les derniers points de chaque capteur, on definit son orientation en prenant l'angle entre l'avant derniere et la derniere valeur
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

    #Formatage du resultat
    result <- data.table(TS)
    result <- cbind(result,XLambertCapt0,
                    YLambertCapt0,
                    XLambertCapt3,
                    YLambertCapt3,
                    XLambertCapt4,
                    YLambertCapt4,
                    DistancePoints)

    result
}
