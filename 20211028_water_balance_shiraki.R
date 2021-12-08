################
# Water Balance Model Shiraki Plateu
# ML Fischer 11/2021
################

remove(list=ls()) 
library(raster)
library(rgdal)
setwd("C:/Users/Marcus/OneDrive/Desktop/Georgien/20211028_model_new")

############################################################################
############################################################################
# Part 1
# Calculate R_sw:
# Berger, A., Loutre, M.F. and Yin Q. (2010), Total irradiation during any time interval of the year us-ing elliptic integrals, Quaternary Science Reviews, 29, 1968 - 1982, doi:10.1016/j.quascirev.2010.05.007
library(palinsol)
lat = 41.3 * pi / 180
orbit=c(eps=  23.446 * pi/180., ecc= 0.016724, varpi= (102.04 - 180)* pi/180. ) #modern day orbit
Insol_d1d2 (orbit,d1=1,d2=1,lat=41*pi/180,avg=TRUE)# W/m2 


############################################################################
############################################################################
# Part 2 calcute ET
source("ET_model_function.R")


#roughness length forest interval:
#according to bergner et al., 2003
z0_f_u <- 0.5
z0_f_l <- 0.3

#roughness length grassland wet interval:
z0_gw_u <- 0.06
z0_gw_l <- 0.02

#roughness length lake
z0_l <- 0.0002

#roughness length grassland dry interval:
z0_gd_u <- 0.03
z0_gd_l <- 0.008

# surface drag coefficient
# according to Rowntree In; Schmugge und Andre p.9f, Bookhagen 2001, Bergner 2003


#forest cda:
z0       <- z0_f_u
u_a      <- (0.4*1.87)/ log(10/z0) #friction velocity (windspeed =1.87m/s at 10 m, karman konstant = 0.4)
r_a_a    <- log(10/z0)/(u_a * 0.4) #atmospheric resitance
cd_a     <- 0.4*0.4 * (log(r_a_a /(z0))^-2) #surface drag coefficient
cda_f_u   <- cd_a

z0       <- z0_f_l
u_a      <- (0.4*1.87)/ log(10/z0) #friction velocity (windspeed =1.87m/s at 10 m, karman konstant = 0.4)
r_a_a    <- log(10/z0)/(u_a * 0.4) #atmospheric resitance
cd_a     <- 0.4*0.4 * (log(r_a_a /(z0))^-2) #surface drag coefficient
cda_f_l   <- cd_a

#wet grasland cda:
z0       <- z0_gw_u
u_a      <- (0.4*1.87)/ log(10/z0) #friction velocity (windspeed =1.87m/s at 10 m, karman konstant = 0.4)
r_a_a    <- log(10/z0)/(u_a * 0.4) #atmospheric resitance
cd_a     <- 0.4*0.4 * (log(r_a_a /(z0))^-2) #surface drag coefficient
cda_gw_u   <- cd_a

z0       <- z0_gw_l
u_a      <- (0.4*1.87)/ log(10/z0) #friction velocity (windspeed =1.87m/s at 10 m, karman konstant = 0.4)
r_a_a    <- log(10/z0)/(u_a * 0.4) #atmospheric resitance
cd_a     <- 0.4*0.4 * (log(r_a_a /(z0))^-2) #surface drag coefficient
cda_gw_l   <- cd_a

#dry grasland cda:
z0       <- z0_gd_u
u_a      <- (0.4*1.87)/ log(10/z0) #friction velocity (windspeed =1.87m/s at 10 m, karman konstant = 0.4)
r_a_a    <- log(10/z0)/(u_a * 0.4) #atmospheric resitance
cd_a     <- 0.4*0.4 * (log(r_a_a /(z0))^-2) #surface drag coefficient
cda_gd_u   <- cd_a

z0       <- z0_gd_l
u_a      <- (0.4*1.87)/ log(10/z0) #friction velocity (windspeed =1.87m/s at 10 m, karman konstant = 0.4)
r_a_a    <- log(10/z0)/(u_a * 0.4) #atmospheric resitance
cd_a     <- 0.4*0.4 * (log(r_a_a /(z0))^-2) #surface drag coefficient
cda_gd_l   <- cd_a

#lake cda:
z0       <- z0_l
u_a      <- (0.4*1.87)/ log(10/z0) #friction velocity (windspeed =1.87m/s at 10 m, karman konstant = 0.4)
r_a_a    <- log(10/z0)/(u_a * 0.4) #atmospheric resitance
cd_a     <- 0.4*0.4 * (log(r_a_a /(z0))^-2) #surface drag coefficient
cda_l   <- cd_a


######FOREST ET early holocene
ETa_EA_forest_u <- ETa( t_a_ld    = 286.2, 
                 emis_ld   = 0.95, 
                 albedo_ld = 0.1,
                 rh        = 0.63,
                 f_ld      = 0.75,
                 cc        = 0.53,
                 ws        = 1.87,
                 a         = 0.39,
                 b         = 0.38,
                 a_2       = 0.22,
                 b_2       = 2.00,
                 cds       = cda_f_u,
                 p         = 93982,
                 r_swc     = 325)

ETa_EA_forest_l <- ETa( t_a_ld    = 286.2, 
                      emis_ld   = 0.95, 
                      albedo_ld = 0.1,
                      rh        = 0.63,
                      f_ld      = 0.5,
                      cc        = 0.53,
                      ws        = 1.87,
                      a         = 0.39,
                      b         = 0.38,
                      a_2       = 0.22,
                      b_2       = 2.00,
                      cds       = cda_f_l,
                      p         = 93982,
                      r_swc     = 325)

######WET GRASSLAND ET early holocene (iron age?)
ETa_EA_wetgrassland_u <- ETa( t_a_ld    = 286.2, 
                        emis_ld   = 0.9, 
                        albedo_ld = 0.15,
                        rh        = 0.63,
                        f_ld      = 0.3,
                        cc        = 0.53,
                        ws        = 1.87,
                        a         = 0.39,
                        b         = 0.38,
                        a_2       = 0.22,
                        b_2       = 2.00,
                        cds       = cda_gw_u,
                        p         = 93982,
                        r_swc     = 325)

ETa_EA_wetgrassland_l <- ETa( t_a_ld    = 286.2, 
                        emis_ld   = 0.9, 
                        albedo_ld = 0.15,
                        rh        = 0.63,
                        f_ld      = 0.3,
                        cc        = 0.53,
                        ws        = 1.87,
                        a         = 0.39,
                        b         = 0.38,
                        a_2       = 0.22,
                        b_2       = 2.00,
                        cds       = cda_gw_l,
                        p         = 93982,
                        r_swc     = 325)


######DRX GRASSLAND ET ?late bronce?
ETa_EA_drygrassland_u <- ETa( t_a_ld    = 287.5, #1.5°C warmer
                              emis_ld   = 0.9, 
                              albedo_ld = 0.15,
                              rh        = 0.63,
                              f_ld      = 0.3,
                              cc        = 0.53,
                              ws        = 1.87,
                              a         = 0.39,
                              b         = 0.38,
                              a_2       = 0.22,
                              b_2       = 2.00,
                              cds       = cda_gd_u,
                              p         = 93982,
                              r_swc     = 325)

ETa_EA_drygrassland_l <- ETa( t_a_ld    = 287.5, #1.5°C warmer
                              emis_ld   = 0.9, 
                              albedo_ld = 0.15,
                              rh        = 0.63,
                              f_ld      = 0.3,
                              cc        = 0.53,
                              ws        = 1.87,
                              a         = 0.39,
                              b         = 0.38,
                              a_2       = 0.22,
                              b_2       = 2.00,
                              cds       = cda_gd_l,
                              p         = 93982,
                              r_swc     = 325)


######lake
ETa_EA_lake <- ETa( t_a_ld    = 287.5, #1.5°C warmer
                              emis_ld   = 0.90, 
                              albedo_ld = 0.06,
                              rh        = 0.63,
                              f_ld      = 1,
                              cc        = 0.53,
                              ws        = 1.87,
                              a         = 0.39,
                              b         = 0.38,
                              a_2       = 0.22,
                              b_2       = 2.00,
                              cds       = cda_l,
                              p         = 93982,
                              r_swc     = 325)

######DRX GRASSLAND ET modern-day
ETa_EA_drygrassland_MD_u <- ETa( t_a_ld    = 286.2, 
                              emis_ld   = 0.9, 
                              albedo_ld = 0.15,
                              rh        = 0.63,
                              f_ld      = 0.3,
                              cc        = 0.53,
                              ws        = 1.87,
                              a         = 0.39,
                              b         = 0.38,
                              a_2       = 0.22,
                              b_2       = 2.00,
                              cds       = cda_gd_u,
                              p         = 93982,
                              r_swc     = 325)

ETa_EA_drygrassland_MD_l <- ETa( t_a_ld    = 286.2, 
                              emis_ld   = 0.9, 
                              albedo_ld = 0.15,
                              rh        = 0.63,
                              f_ld      = 0.3,
                              cc        = 0.53,
                              ws        = 1.87,
                              a         = 0.39,
                              b         = 0.38,
                              a_2       = 0.22,
                              b_2       = 2.00,
                              cds       = cda_gd_l,
                              p         = 93982,
                              r_swc     = 325)



ETa_EA_drygrassland_MD_780masl_u <- ETa( t_a_ld    = 286.2- (1.5*0.65), 
                                 emis_ld   = 0.9, 
                                 albedo_ld = 0.15,
                                 rh        = 0.63,
                                 f_ld      = 0.3,
                                 cc        = 0.53,
                                 ws        = 1.87,
                                 a         = 0.39,
                                 b         = 0.38,
                                 a_2       = 0.22,
                                 b_2       = 2.00,
                                 cds       = cda_gd_u,
                                 p         = 92323,
                                 r_swc     = 325)

ETa_EA_drygrassland_MD_780masl_l <- ETa( t_a_ld    = 286.2 - (1.5*0.65), 
                                 emis_ld   = 0.9, 
                                 albedo_ld = 0.15,
                                 rh        = 0.63,
                                 f_ld      = 0.3,
                                 cc        = 0.53,
                                 ws        = 1.87,
                                 a         = 0.39,
                                 b         = 0.38,
                                 a_2       = 0.22,
                                 b_2       = 2.00,
                                 cds       = cda_gd_l,
                                 p         = 92323,
                                 r_swc     = 325)

############################################################################
############################################################################
# Part 3 DEM Analysis


# load dem and catchment shapefile
setwd("C:/Users/Marcus/OneDrive/Desktop/Georgien/20211028_model_new/data")
ra <- raster("SRTM_Georgien.tif")
plot(ra)

area.shp <- readOGR(".",layer = "basin_UTM")
plot(area.shp,add=T)
area(area.shp)

#reproject
dem <- projectRaster(ra,crs= area.shp@proj4string)

#clip and masc dem:
area.shp<- spTransform(area.shp,ra@crs@projargs) #transform shapefile 
plot(area.shp,add=T)
dem <- mask(ra,area.shp)
dem <- crop(dem,area.shp)

#show result:
plot(dem)
plot(area.shp,add=T)

#catchment area (m^2)
enszw <- area(area.shp)

#######################
freq<-as.data.frame(freq(dem,digit=1, merge=F))
freq[7,2] <- sum(freq[1:7,2])
freq<-freq[-(1:6),]

# cleaning dataset
freq <- freq[order(freq[,1]),] 

# create a sum array of all frequency heights
# all summed pixels below a height means the area below a certain height
for (i in 1:(length(freq[,1])-1))
{freq[i,3]<-sum(freq[1:i,2])}  

# calcute the area by the resolution of the pixels, km²
freq<-freq[-423,] #delete the last control row
freq[,4]<-freq[,3]*(23^2)/(1000^2) 

# create a sum array of all frequency area-counts
# the volume means that every pixel below a certain value has 1m³ for each column
for (i in 1:(length(freq[,1]))) 
{freq[i,5]<-sum(freq[1:i,3])}   

freq[,6]<-freq[,5]*(23^2)/(1000^3) #volume by height
freq[,7]<-freq[2]*(23^2)/(1000^2)#area per height 
# lake area ratio is the area below height divided throug the catchment size
freq[,8]<- 100*freq[,3]/2173415800 
colnames(freq)<-c("height","pix_by_height","pix_below_height","km²_below_height","pix_vol_below","km³_below","area_by_height","lake area ration")

plot(freq[,1]~freq[,4],type="l",ylab="Height in m",xlab= "Area in km²")
abline(h=700)
abline(h=600)

####





############################################################################
############################################################################
# Part 4 water volume per year

A_basin <- area(area.shp)/1000000

plot(dem)
dem_c<-dem
dem_c[which(dem_c[]<600)]<-1
dem_c[which(dem_c[]>=600)]<-2
A_f600 <- as.numeric(freq(dem_c)[2,2]/(freq(dem_c)[2,2]+freq(dem_c)[1,2])*A_basin)


dem_c<-dem
dem_c[which(dem_c[]<700)]<-1
dem_c[which(dem_c[]>=700)]<-2
A_f700 <- as.numeric(freq(dem_c)[2,2]/(freq(dem_c)[2,2]+freq(dem_c)[1,2])*A_basin)


dem_c<-dem
dem_c[which(dem_c[]>578)]<-1
dem_c[which(dem_c[]<=578)]<-2
A_f578 <- as.numeric(freq(dem_c)[2,2]/(freq(dem_c)[2,2]+freq(dem_c)[1,2])*A_basin)


A_l_u <- 10
A_l_l <- 5

#A_l_max

#?iron age:
V_p700 <- 700/1000000* A_basin #km^3
V_p800 <- 800/1000000* A_basin

V_p700*1000 # hm^3
V_p800*1000

#600m timberline:
V_f600_u <- A_f600 * ETa_EA_forest_u/1000000
V_f600_l <- A_f600 * ETa_EA_forest_l/1000000

V_f600_u*1000
V_f600_l*1000

V_g600_u <- (A_basin-A_f600) * ETa_EA_wetgrassland_u/1000000
V_g600_l <- (A_basin-A_f600) * ETa_EA_wetgrassland_l/1000000

V_g600_u*1000
V_g600_l*1000

#700m timberline:
V_f700_u <- A_f700 * ETa_EA_forest_u/1000000
V_f700_l <- A_f700 * ETa_EA_forest_l/1000000

V_f700_u*1000
V_f700_l*1000

V_g700_u <- (A_basin-A_f700) * ETa_EA_wetgrassland_u/1000000
V_g700_l <- (A_basin-A_f700) * ETa_EA_wetgrassland_l/1000000

V_g700_u*1000
V_g700_l*1000

#?bronce age:
V_p700 <- 700/1000000* A_basin
V_p800 <- 800/1000000* A_basin

V_g_u <- A_basin * ETa_EA_drygrassland_u/1000000
V_g_l <- A_basin * ETa_EA_drygrassland_l/1000000

V_g_u*1000
V_g_l*1000

V_l_u <- A_l_u * ETa_EA_lake/1000000
V_l_l <- A_l_l * ETa_EA_lake/1000000

V_l_u*1000
V_l_l*1000



V_A_f578 <- A_f578 * ETa_EA_lake/1000000
V_A_f578*1000



#modern-day

V_p500 <- 500/1000000* A_basin

V_p500*1000


V_gm_u <- A_basin * ETa_EA_drygrassland_MD_u/1000000
V_gm_l <- A_basin * ETa_EA_drygrassland_MD_l/1000000

V_gm_u*1000
V_gm_l*1000



A_basin *ETa_EA_drygrassland_MD_780masl_u/1000
A_basin *ETa_EA_drygrassland_MD_780masl_l/1000


############################################################################
############################################################################
# Part 5 monthly water balance

get(load("cli.RData"))
cli_h <- cli


cli[2,]<- cli[2,]* 500/sum(cli[2,])

cli_h[2,]<- cli_h[2,]* 700/sum(cli_h[2,])
cli_h<- rbind(cli_h, cli_h[2,]* 800/sum(cli_h[2,]))
cli_h[1,]<- cli_h[1,]+1.5

library(SPEI)

ETp<- as.vector(thornthwaite(cli[1,], lat= 41, na.rm = FALSE))
cli <- rbind(cli,ETp)

ETp<- as.vector(thornthwaite(cli_h[1,], lat= 41, na.rm = FALSE))
cli_h <- rbind(cli_h,ETp)

sum(cli[2,]) #PRE
sum(cli[3,]) #ETp


sum(cli_h[2,]) #PREmin
sum(cli_h[3,]) # PREmax
sum(cli_h[4,]) # ETp


wd <- cli[2,]- cli[3,]

cli_max_storage <- sum(wd[which(wd[]>0)])
cli_neg_storage <- sum(wd[which(wd[]<0)])


res<-rep(0,12)
cli <- rbind(cli,wd,res)

du<- c(9:12,1:9)
cli


cli 

for(i in 1:12){
  cli[5,du[i+1]]<-cli[5,du[i]] + cli[4,du[i+1]]
  cli
}
cli[5,]<-cli[5,] + cli[5,9]





sum(cli_h[4,])
sum(cli_h[3,])

#####
wd <- cli_h[2,]- cli_h[4,]
cli_h_max_storage_1 <- sum(wd[which(wd[]>0)])

wd <- cli_h[3,]- cli_h[4,]
cli_h_max_storage_2 <- sum(wd[which(wd[]>0)])


wd <- cli_h[2,]- cli_h[4,]
cli_h_neg_storage_1 <- sum(wd[which(wd[]<0)])

wd <- cli_h[3,]- cli_h[4,]
cli_h_neg_storage_2 <- sum(wd[which(wd[]<0)])
####

#####
wd <- cli_h[2,]- cli[3,]
cli_h_max_storage_1_e <- sum(wd[which(wd[]>0)])

wd <- cli_h[3,]- cli[3,]
cli_h_max_storage_2_e <- sum(wd[which(wd[]>0)])


wd <- cli_h[2,]- cli[3,]
cli_h_neg_storage_1_e <- sum(wd[which(wd[]<0)])

wd <- cli_h[3,]- cli[3,]
cli_h_neg_storage_2_e <- sum(wd[which(wd[]<0)])





####
Vol <- (cli_h_max_storage_2- 100) / 1000000 * 300
###
ETp<- as.vector(thornthwaite(cli_h[1,], lat= 41, na.rm = FALSE))
cli_h <- rbind(cli_h,ETp)

mon <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
###
#pdf("20200926_Figure_xx_monthly_raw.pdf", width = 6.25, height = 3.5,pointsize=8)

par(mfrow=c(1,3))
par(oma=c(1,2,2,2))


plot(cli_h[2,],type="l",ylim=c(0,200),col="blue"
     ,main="Early to Mid Holocene"
     ,ylab="Water Flux (mm)"
     ,xlab="Month of the Year"
     ,axes=F)
lines(cli_h[3,],type="l",col="blue")
lines(cli[3,],type="l",col="red")

axis(2,las=1)
axis(1,las=1,labels=mon,at=1:12)

mtext(text = paste(round(cli_h_max_storage_1_e),"to" ,round(cli_h_max_storage_2_e),"mm"),
      side = 1,#side 1 = bottom
      line = -15,col="blue",adj=0)

mtext(text = paste(round(cli_h_neg_storage_1_e),"to" ,round(cli_h_neg_storage_2_e),"mm"),
      side = 1,#side 1 = bottom
      line = -15,col="red",adj=1)
m<- 12
polygon(y= c(cli_h[3,1:m], cli[3,m:1]), x= c(1:m, m:1) ,border=NA,col=rgb(0, 0, 0,0.25))

m<- 12
polygon(y= c(cli_h[2,1:m], cli_h[3,m:1]), x= c(1:m, m:1) ,border=NA,col=rgb(0, 0, 0,0.25))

plot(cli_h[2,],type="l",ylim=c(0,200),col="blue"
     ,main="Late Bronze/ Early Iron Age"
     ,ylab="Water Flux (mm)"
     ,xlab="Month of the Year"
     ,axes=F)
lines(cli_h[3,],type="l",col="blue")
lines(cli_h[4,],type="l",col="red")
#lines(cli[3,],type="l",col="green")

axis(2,las=1)
axis(1,las=1,labels=mon,at=1:12)

mtext(text = paste(round(cli_h_max_storage_1),"to" ,round(cli_h_max_storage_2),"mm"),
      side = 1,#side 1 = bottom
      line = -15,col="blue",adj=0)

mtext(text = paste(round(cli_h_neg_storage_1),"to" ,round(cli_h_neg_storage_2),"mm"),
      side = 1,#side 1 = bottom
      line = -15,col="red",adj=1)
m<- 12
polygon(y= c(cli_h[3,1:m], cli_h[4,m:1]), x= c(1:m, m:1) ,border=NA,col=rgb(0, 0, 0,0.25))

m<- 12
polygon(y= c(cli_h[2,1:m], cli_h[3,m:1]), x= c(1:m, m:1) ,border=NA,col=rgb(0, 0, 0,0.25))

plot(cli[2,],type="l",ylim=c(0,200),col="blue"
     ,main="Modern Day"
     ,ylab="Water Flux (mm)"
     ,xlab="Month of the Year"
     ,axes=F)
lines(cli[3,],type="l",col="red")

axis(2,las=1)
axis(1,las=1,labels=mon,at=1:12)

mtext(text = paste(round(cli_max_storage),"mm"),
      side = 1,#side 1 = bottom
      line = -15,col="blue",adj=0)

mtext(text = paste(round(cli_neg_storage),"mm"),
      side = 1,#side 1 = bottom
      line = -15,col="red",adj=1)
m<- 12
polygon(y= c(cli[2,1:m], cli[3,m:1]), x= c(1:m, m:1) ,border=NA,col=rgb(0, 0, 0,0.25))


polygon(y= c(cli[2,1:4], cli[3,4:1]), x= c(1:4, 4:1) ,border=NA,col=rgb(0, 0, 1,0.5))
polygon(y= c(cli[2,10:12], cli[3,12:10]), x= c(10:12, 12:10) ,border=NA,col=rgb(0, 0, 1,0.5))
polygon(y= c(cli[2,4:10], cli[3,10:4]), x= c(4:10, 10:4) ,border=NA,col=rgb(1, 0, 0,0.5))


legend("topleft", legend=c("potential Evapotranspiration", "Precipitation"),
       col=c("red", "blue"), lty=1, bty="n")







