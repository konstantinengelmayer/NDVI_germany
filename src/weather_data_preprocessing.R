# Get the weather data from one station: 
library(rdwd) # for the dwd download
#TODO: Get data from dwd
link <- selectDWD(id="001420", res="daily", var="kl", per="rescent") # Get the link from the DWD Website Frankfurt
file <- dataDWD(link, read=F) # Creating a file
df <- readDWD(file, varnames = T) # Reading the file with all values and names

head(df) # have a look at the data
tail(df)

# Calculate the Evapotranspiration if necessary:
library(simET)
library(lubridate)
library(solarPos)
library(openxlsx)

unique(is.na(df$SDK.Sonnenscheindauer)) # look if SDK has been measured
#df <- na.omit(df) # delete all NAs

start_date <- ymd("2022-07-18") # Create start and end date
end_date <- ymd("2024-01-18")
date_vector <- seq(start_date, end_date, by = "day") # create date vector
day_of_year <- yday(date_vector) # with the days
day_of_year
df$day_of_year <- day_of_year # put the day to dataframe df

RA <- cal_extraterrestrialRadiation_for_daily(df$day_of_year, 50.110924*(pi/180)) # Calculate the Radiation !different values for every location!
RA
SDK <- df$SDK.Sonnenscheindauer
SDK
Sunset_Angle <- cal_sunsetHourAngle((50.110924*(pi/180)), cal_solarDeclination(day_of_year)) # calculate the sunset !different values for every location!
Sunset_Angle
Daylight <- cal_daylightHours(Sunset_Angle) # calculate daylight hour !different values for every location!
Daylight
Solar_Radiation <- dplyr::mutate(df, 
                    Solar_Radiation=cal_solarRadiation(as=0.25,
                    bs= 0.5,
                    n=SDK.Sonnenscheindauer, 
                    N=Daylight,
                    Ra=RA))


Solar_Radiation <- cal_solarRadiation(as=0.25,
                                      bs= 0.5,
                                      n=SDK, 
                                      N=Daylight,
                                      Ra=RA)
Solar_Radiation

Evapotrans <- dplyr::mutate(df,
                            ET0=cal_ET0_from_PM_for_daily(Latitude=50.110924,
                                                          Altitude = 110,
                                                          J=day_of_year,
                                                          Tmax=TXK.Lufttemperatur_Max,
                                                          Tmin=TNK.Lufttemperatur_Min,
                                                          Rs=Solar_Radiation,
                                                          RHmean=UPM.Relative_Feuchte,
                                                          Wind=FM.Windgeschwindigkeit)) # Calculate the Evapotranspiration !different values for every location!
Evapotranspiration <- Evapotrans$ET0
Evapotranspiration
df$Evapotranspiration <- Evapotranspiration
Wasserbilanz <- Evapotrans$RSK.Niederschlagshoehe-Evapotrans$ET0 # Calculate the "Wasserbilanz"
df$Wasserbilanz <- Wasserbilanz
head(df)

# write.table(df, file="Wetterdaten_Evapotranspiration_1420", sep=",", dec=".", row.names=FALSE)




######################## This is the script from Reudenbach ##########################################
# https://github.com/gisma/meteoGermany/tree/main/src

# Main control script
#'
#' @description controls the main prediction run. 
#'              Download and prepare climate and DEM data 
#'              perform a Kriging with autovariogram.
#'              A two step value correction is performed 
#'              to keep the Interpolation in valid ranges. 
#'              Finally the corrected data is extracted by communities
#'              
#' Sunshine kriging https://rmets.onlinelibrary.wiley.com/doi/pdf/10.1017/S1350482706002362
#' @author Chris Reudenbach creuden@gmail.com

# ---- setup project ----
#devtools::install_github("envima/envimaR")
library(envimaR)
library(rprojroot)
appendProjectDirList = c("data/data_lev0/GhcnDaily",
                         "data/data_lev0/GhcnMonthly")
root_folder = find_rstudio_root_file()
source(file.path(root_folder, "000_setup.R"))

# ---- default arguments ----
crs = raster::crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
epsg=3035
res=500     # resolution of prediction DEM
startDate = "2003-01-01"
endDate = "2021-12-31"
type= "historical"    #c("historical","recent") # recent means rolling the last 500 days
getDEM = FALSE        # download and prepare DEM data (only needed once)
getClimate = FALSE    # download climate data (usually only done once)
minStations = 5       # minimum number of accepted stations

# ---- start processing ----
# ---- prepare the climate and auxiliary data----
source(file.path(envrmt$path_src,"main_prepare_data.R"))

dat_list = sort(as.character(unique(cVar.sf$MESS_DATUM)))[1:length(unique(cVar.sf$MESS_DATUM))]
cVar ="SDK"
for (cVar in c("TXK","TNK","TMK","SDK","PM","UPM")){
  
  matrix_of_sums <- parallel::mclapply( seq_along(dat_list), function(n){
    currentDate = dat_list[n]
    if (as.Date(currentDate) >= as.Date(startDate) & as.Date(currentDate) <= as.Date(endDate)){
      cd= substr(currentDate,1,10)
      if(!file.exists(paste0(envrmt$path_data_lev1,"/",cVar,"/",cd,"_",cVar,".tif"))){
        cVar.sf.day = cVar.sf[as.character(cVar.sf$MESS_DATUM) == as.Date(currentDate),]
        if (cVar == "SDK") {
          dt=suncalc::getSunlightTimes(date = as.Date(currentDate), lat = 51.0, lon = 9.0, tz = "UTC")
          td=dt$sunset-dt$sunrise
          maxDaylight= ceiling(as.numeric(unlist(stringr::str_split(td,"Time difference of "))))
          cVar.sf.day$tmp=NA
          dat = cVar.sf.day %>% mutate(tmp = replace(!!sym(cVar), as.numeric(!!sym(cVar)) == -999, NA))
          dat = cVar.sf.day %>% mutate(tmp = replace(!!sym(cVar), as.numeric(!!sym(cVar)) > maxDaylight, maxDaylight))
        } else if (cVar == "PM") {
          cVar.sf.day$tmp=NA
          dat = cVar.sf.day %>% mutate(tmp = replace(!!sym(cVar), as.numeric(!!sym(cVar)) == -999, NA))
          dat = cVar.sf.day %>% mutate(tmp = replace(!!sym(cVar), as.numeric(!!sym(cVar)) > 1060.6, 1060.6))
          dat = cVar.sf.day %>% mutate(tmp = replace(!!sym(cVar), as.numeric(!!sym(cVar)) < 954.9, 954.9))
        } else if (cVar == "UPM") {
          cVar.sf.day$tmp=NA
          dat = cVar.sf.day %>% mutate(tmp = replace(!!sym(cVar), as.numeric(!!sym(cVar)) == -999, NA))
          dat = cVar.sf.day %>% mutate(tmp = replace(!!sym(cVar), as.numeric(!!sym(cVar)) > 100, 100))
          dat = cVar.sf.day %>% mutate(tmp = replace(!!sym(cVar), as.numeric(!!sym(cVar)) < 0, 0))
        }else if (cVar == "TXK" | cVar == "TNK" | cVar == "TMK") {
          cVar.sf.day$tmp=NA
          dat = cVar.sf.day %>% mutate(tmp = replace(!!sym(cVar), as.numeric(!!sym(cVar)) == -999, NA))
          dat = cVar.sf.day %>% mutate(tmp = replace(!!sym(cVar), as.numeric(!!sym(cVar)) > 42, 42))
          dat = cVar.sf.day %>% mutate(tmp = replace(!!sym(cVar), as.numeric(!!sym(cVar)) < -46.0, -46.0))
        }
        
        dat = dat[,c("Stationshoehe","tmp","geometry")]
        dat$tmp=as.numeric(dat$tmp)
        names(dat) = c("Stationshoehe",cVar,"geometry")
        data <- dat %>% drop_na()
        if (nrow(data)>minStations){
          data <- dplyr::distinct(data, geometry, .keep_all = TRUE)
          data = st_transform(data,st_crs(dem))
          st_crs(dem)=3035
          st_crs(data)=3035
          seed=123
          
          vm.auto = automap::autofitVariogram(formula = as.formula(paste(cVar, "~1")),
                                              input_data = data)
          
          #plot(vm.auto)
          seed=123
          tmax.pred <- krige(formula = as.formula(paste(cVar, "~Stationshoehe")),
                             locations = data,#data[sample.int(nrow(data),min(sample_size,nrow(data))),],
                             newdata = dem,
                             model = vm.auto$var_model,
                             debug.level=-1)
          
          stars::write_stars(tmax.pred,paste0(envrmt$path_data_lev1,"/",cVar,"/",cd,"_",cVar,".tif"),overwrite=TRUE,options="COMPRESS=LZW")
          rm(tmax.pred)
          gc()
        } else {
          stars::write_stars(dem*0-9999,paste0(envrmt$path_data_lev1,"/",cVar,"/",cd,"_",cVar,".tif"),overwrite=TRUE,options="COMPRESS=LZW")
        }
      }
      
    }
  }, mc.cores = 10, mc.allow.recursive = TRUE)
}

# final correction and extraction per community
for (cVar in c("TXK","TNK","TMK","SDK","PM","UPM")){
  source(file.path(root_folder, "src/main_script_calculate_communities.R"))
}
