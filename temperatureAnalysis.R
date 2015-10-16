##############################################################################
####PLOT HISTOGRAM OF DAILY AVERAGE TEMPERATURE FOR ALL DATA
##############################################################################

#Variable of month names for naming files
months <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October",
	"November", "December")

#Generate histograms of the average temperature (mean of max and min) for each day of the year throughout
#the entire dataset

SDs <- vector( ) #track the standard deviation
c <- 1
months30 <- c(9, 4, 6, 11)
for(i in 1:12){
	if(i %in% months30){ #months with 30 days
		lim <- 30
	}else if(i == 2){
		lim <- 29
	}else{
		lim <- 31
	}
	for(j in 1:lim){
		if(i < 10){
			if(j < 10){
				d <- paste(0, i, 0, j, sep = "")
			}else{
				d <- paste(0, i, j, sep = "")
			}
		}else{
			if(j < 10){
				d <- paste(i, 0, j, sep = "")
			}else{
				d <- paste(i, j, sep = "")
			}
		}
		t <- s.df.sep[s.df.sep$MONTH_DAY == d, ]
		a <- ((t$MaxTemp.0.1C / 10) + (t$MinTemp.0.1C / 10)) / 2
		png(filename = paste("...", d, ".png", sep = '')) #Save the graph where "..." is the plot name
		hist(a, breaks = seq(floor(min(a)), ceiling(max(a)), .5), xlab = "Degrees C", main = paste("Average Temperature on", d))
		dev.off( )
		SDs[c] <- sd(a)
		c <- c + 1
	}
}

#Plot the various standard deviations
plot(seq(1, 366, 1), SDs, xlab = "Day of the Year", ylab = "Standard Deviation",
	main = "Standard Deviation of Average Temperature by Day of the Year")

##############################################################################
####PLOT COMPARISON HISTOGRAM OF MONTHLY AVERAGE TEMPERATURE BETWEEN BLOCKS
##############################################################################

MEANS <- vector( )
SDS <- vector( )
co <- 1

for(j in 1:12){
	#Data for month j from the earliest block of size 30
	t1 <- vector( )
	for(i in 1869:1898){
		s <- subsetData(j, i, s.df.sep)
		x <- ((s$MaxTemp.0.1C / 10) + (s$MinTemp.0.1C / 10)) / 2
		t1 <- c(t1, x)
	}

	#Data for month j from the earliest block of size 30
	t2 <- vector( )
	for(i in 1985:2014){
		s <- subsetData(j, i, s.df.sep)
		x <- ((s$MaxTemp.0.1C / 10) + (s$MinTemp.0.1C / 10)) / 2
		t2 <- c(t2, x)
	}

	#Track the mean and standard deviation of each resultant histogram
	MEANS[co] <- mean(t1)
	SDS[co] <- sd(t1)
	co <- co + 1
	MEANS[co] <- mean(t2)
	SDS[co] <- sd(t2)
	co <- co + 1

	#Plot the overlapping histograms
	png(filename = paste("...", months[j], ".png", sep = ''))
	hist(t1, breaks = seq(-30, 35, .5), xlab = "Degrees C", main = paste("Average Temperature in ", months[j], sep = ""),
		col=rgb(1,0,0,0.5))
	hist(t2, breaks = seq(-30, 35, .5), col=rgb(0,0,1,0.5), add = T)
	legend("topright", c("1869-1878", "2005-2014"), col = c("red", "blue"), lwd = 2)
	dev.off( )
}

##############################################################################
####PLOT DIFFERENCE IN OBSERVED AND AVERAGE TEMPERATURE FOR EACH DAY
##############################################################################

#Subset the data into all years between 1981 and 2014 (the most recent 30 year block, and a few additional years of data)
x <- s.df.sep[s.df.sep$YEAR >= 1981, ]
x <- x[x$YEAR <= 2014, ]

for(i in 1:12){
	if(i == 2){
		d <- 28
	}else if(i == 9 || i == 4 || i == 6 || i == 11){
		d <- 30
	}else{
		d <- 31
	}
	for(day in 1:d){
		if(i < 10){
			if(day < 10){
				y <- paste(0, i, 0, day, sep = "")
			}else{
				y <- paste(0, i, day, sep = "")
			}
		}else{
			if(day < 10){
				y <- paste(i, 0, day, sep = "")
			}else{
				y <- paste(i, day, sep = "")
			}
		}
		s <- x[x$MONTH_DAY == y, 7] / 10
		a <- s - mean(s)
		png(filename = paste("C:/Users/Morgan/Desktop/School/Old/Summer15/IndepStudy/DegreesFromAverage/", months[i],
			day, ".png", sep = ''))
		plot(a, xaxt = 'n', xlab = "Year", ylab = "Degrees Difference (C)", main = paste("Degrees away from climatological
			average for this day\n", months[i], day))
		abline(0, 0)
		axis(1, at = 1:34, label = seq(1981, 2014, 1))
		dev.off( )
	}
}




















