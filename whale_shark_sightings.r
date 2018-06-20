require(lubridate); require(raster); require(RNetCDF); require(HMMoce)
require(fields)

# explore some sst data
sst.dir <- '~/Documents/WHOI/Data/WhaleSharks/Natl/EnvData/sst/'
#dir.create(file.path(sst.dir), recursive = TRUE, showWarnings = FALSE)
setwd('~/Documents/WHOI/Data/WhaleSharks/Natl/')

# load sightings table
dat <- read.table(file='encounterSearchResults_export_jcochrane.csv', sep=',', header=T)
sp.lim <- list(lonmin=-80, lonmax=-60,
               latmin=35, latmax=45)

# use HMMoce::get.bath.data() to download bathymetry
bathy <- raster::raster('~/Documents/WHOI/Data/Swordfish/batch/sword_ctag_bathy_highres.grd')
c200 <- rasterToContour(bathy, levels=-200)

dat$sst <- NA
dat$depth <- NA
dat$dist_to_200m <- NA

for (i in 1:nrow(dat)){
  udates <- as.Date(paste(dat$Year.Identified[i],'-',dat$Month.Identified[i],'-',dat$Day.Identified[i], sep=''), format='%Y-%m-%d')
  
  HMMoce::get.env(udates, filename='whaleshark', type = 'sst', sst.type='ghr', spatLim = sp.lim, save.dir = sst.dir)
  
  nc <- open.nc(paste(sst.dir, 'whaleshark_', udates,'.nc',sep=''))
  lat <- as.numeric(var.get.nc(nc, 'latitude'))
  lon <- as.numeric(var.get.nc(nc, 'longitude'))
  sst <- var.get.nc(nc, 'SST')
  xlims <- c(sp.lim[[1]], sp.lim[[2]]); ylims <- c(sp.lim[[3]], sp.lim[[4]])
  lon.idx <- c(which.min((lon - xlims[1])^2):which.min((lon - xlims[2])^2))
  lat.idx <- c(which.min((lat - ylims[1])^2):which.min((lat - ylims[2])^2))
  
  sst <- flip(raster(t(sst[lon.idx,lat.idx])*(9/5)+32, xmn=min(lon[lon.idx]), xmx=max(lon[lon.idx]),
                     ymn=min(lat[lat.idx]), ymx=max(lat[lat.idx])),2)
  #sst.rnge <- cellStats(sst, 'range')
  sst.rnge <- c(58, 78)
  sst.breaks <- seq(sst.rnge[1], sst.rnge[2], length.out=201)
  sst.mid <- sst.breaks[1:(length(sst.breaks)-1)]
  sst.col <- jet.colors(length(sst.breaks)-1) #[as.vector((dataT))]
  
  bathy.rnge <- c(0, -3000)
  bathy.breaks = seq(bathy.rnge[1], bathy.rnge[2], length.out=201)
  bathy.mid = bathy.breaks[1:(length(bathy.breaks)-1)]
  gry.pal <- colorRampPalette(rev(brewer.pal(9, 'Greys')[1:6]))#(100))
  bathy.col <- gry.pal(length(bathy.breaks)-1)
  
  add.alpha <- function(COLORS, ALPHA){
    if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
    RGB <- col2rgb(COLORS, alpha=TRUE)
    RGB[4,] <- round(RGB[4,]*ALPHA)
    NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
    return(NEW.COLORS)
  }
  sst.col <- add.alpha(sst.col, .6)
  
  png(paste('sst2_',udates,'.png',sep=''),
      width=4000, height=3400, res=400, bg='white')
  nf <- layout(matrix(c(1,2),1,2,byrow=T), widths=c(5,1), heights=c(5))
  plot(bathy,xlim=c(-72,-68.5), ylim=c(38.5,41.5), col=bathy.col)
  image(sst, add=T, col=sst.col, breaks=sst.breaks, zlim=c(55,85), xlim=c(-72,-69), ylim=c(38.5,41.5), xlab='', ylab='')
  image(bathy,xlim=c(-77,-65), ylim=c(35,43))
  #lines(c200, lty=2)
  contour(bathy, lty=2, levels=-100, add=T)
  contour(sst, levels=73.5, add=T, col='black')
  points(-68.1667, 40.5000)
  points(dat$Longitude[i], dat$Latitude[i], pch=16)
  text(dat$Longitude[i], dat$Latitude[i]-.5, paste(udates))
  text(dat$Longitude[i], dat$Latitude[i]-1, paste(round(raster::extract(sst, cbind(dat$Longitude[i], dat$Latitude[i])),2),'F'))
  world(add=T, fill=T, col='grey80',border='grey80')
  
  title(paste('SST from ', udates, sep=''))
  
  image(1, sst.mid, t(as.matrix(sst.mid)), breaks=sst.breaks, col=sst.col, axes=FALSE, xlab="",
        ylab=parse(text=paste('Temperature', "*degree~F", sep="")))
  axis(2)#, ylab=parse(text=paste('Temperature', "*degree~F", sep="")));box();
  box()
  dev.off()
  
  dat$sst[i] <- round(raster::extract(sst, cbind(dat$Longitude[i], dat$Latitude[i])),2)
  dat$depth[i] <- round(raster::extract(bathy, cbind(dat$Longitude[i], dat$Latitude[i])),0)
  
}



dat1 <- dat
coordinates(dat1) <- ~ Longitude + Latitude
proj4string(dat1) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
crs.proj <- CRS('+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs')
dat1 <- spTransform(dat1, crs.proj)
c200 <- spTransform(c200, crs.proj)
dat$dist_to_200m <- c(gDistance(dat1, c200, byid=T)) * 5.3996e-4 # nautical miles
